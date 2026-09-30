/*
  Copyright (C) 2026 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file LICENSE.  If not see
  <http://www.gnu.org/licenses/>.
*/


#include <aspect/postprocess/melting_rate_instantaneous.h>
#include <aspect/utilities.h>
#include <aspect/simulator.h>
#include <aspect/material_model/interface.h>
#include <aspect/global.h>

#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/mpi.h>
#include <deal.II/fe/fe_values.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>


namespace aspect
{
  namespace Postprocess
  {
    namespace
    {
      /**
       * The names of the chemical components this diagnostic works with. They
       * have to match the compositional field names used by the material model.
       */
      const std::vector<std::string> chemical_component_names = {"dunite", "morb", "cmorb"};

      /**
       * The number of chemical components.
       */
      constexpr unsigned int n_chemical_components = 3;

      /**
       * Number of scalars written per evaluation point after the two
       * coordinates: total Gamma, the component-wise Gamma, the tie-line
       * spread and the three tie-line coordinates.
       */
      constexpr unsigned int n_output_scalars = 1 + n_chemical_components + 1
                                                + n_chemical_components + 1
                                                + 1 + n_chemical_components;
    }


    template <int dim>
    std::pair<std::string,std::string>
    MeltingRateInstantaneous<dim>::execute (TableHandler &statistics)
    {
      AssertThrow(this->include_melt_transport(),
                  ExcMessage("The 'melting rate instantaneous' postprocessor "
                             "requires 'Include melt transport = true'."));

      const Introspection<dim> &introspection = this->introspection();

      // ------------------------------------------------------------------
      // Locate the fields we need.
      // ------------------------------------------------------------------
      AssertThrow(introspection.compositional_name_exists("porosity"),
                  ExcMessage("The 'melting rate instantaneous' postprocessor requires a "
                             "compositional field named 'porosity'."));

      for (const auto &component : chemical_component_names)
        {
          AssertThrow(introspection.compositional_name_exists(component),
                      ExcMessage("The 'melting rate instantaneous' postprocessor requires a "
                                 "compositional field named '" + component + "'."));
          AssertThrow(introspection.compositional_name_exists(component + "_liquid"),
                      ExcMessage("The 'melting rate instantaneous' postprocessor requires a "
                                 "compositional field named '" + component + "_liquid'."));
        }

      const unsigned int porosity_index = introspection.compositional_index_for_name("porosity");

      std::array<unsigned int, n_chemical_components> bulk_index;
      std::array<unsigned int, n_chemical_components> liquid_index;
      for (unsigned int i = 0; i < n_chemical_components; ++i)
        {
          bulk_index[i]   = introspection.compositional_index_for_name(chemical_component_names[i]);
          liquid_index[i] = introspection.compositional_index_for_name(chemical_component_names[i] + "_liquid");
        }

      const unsigned int fluid_pressure_component
        = introspection.variable("fluid pressure").first_component_index;

      const FEValuesExtractors::Scalar fluid_pressure_extractor(fluid_pressure_component);
      const FEValuesExtractors::Vector fluid_velocity_extractor(
        introspection.variable("fluid velocity").first_component_index);

      // ------------------------------------------------------------------
      // Evaluate at the center of every cell. At the center of a Q2 element
      // the finite element gradient is the symmetric, central derivative of
      // the solution; at cell interfaces the one-sided derivatives of the two
      // adjacent cells differ, so evaluating there would make the result
      // ambiguous. Each cell therefore contributes exactly one point.
      // ------------------------------------------------------------------
      const QMidpoint<dim> quadrature;
      const unsigned int n_q_points = quadrature.size();

      FEValues<dim> fe_values(this->get_mapping(),
                              this->get_fe(),
                              quadrature,
                              update_values | update_gradients | update_quadrature_points);

      // The advection system discretizes the bulk concentration equation with
      // a Gauss quadrature of degree
      //   composition_degree + (stokes_velocity_degree+1)/2
      // (see assemble_advection_system() in assembly.cc) and uses the
      // *arithmetic cell average* of div(V) in both the c_bar*div(V) and
      // c_l*div(V) terms -- see the comment "average divergence u over the
      // cell" in melt.cc. Reproduce exactly that discretization so that the
      // perturbation direction is the same cell-averaged rate of change of the
      // bulk composition that the model itself advances.
      const unsigned int advection_quadrature_degree =
        introspection.polynomial_degree.compositional_fields
        + (introspection.polynomial_degree.velocities + 1) / 2;
      const QGauss<dim> advection_quadrature(advection_quadrature_degree);
      const unsigned int n_advection_q_points = advection_quadrature.size();
      FEValues<dim> advection_fe_values(this->get_mapping(),
                                        this->get_fe(),
                                        advection_quadrature,
                                        update_values | update_gradients | update_quadrature_points);

      const unsigned int n_compositional_fields = introspection.n_compositional_fields;

      // Tie-line diagnostic. For the "uniform phase" calibration case the bulk
      // composition is constructed as c_bar = (1-f) S + f L with constant end
      // members S and L, so the melt fraction recovered from each component,
      //   f_i = (c_bar_i - S_i) / (L_i - S_i),
      // must be identical for all i. Their spread therefore measures directly
      // by how much the state violates the tie line, i.e. by how much the
      // equilibrium has been broken. This is only meaningful on that case but
      // harmless (zero) otherwise.
      AssertThrow(tie_line_solid.size() == n_chemical_components
                  && tie_line_liquid.size() == n_chemical_components,
                  ExcMessage("The 'Tie line solid composition' and 'Tie line liquid composition' "
                             "parameters of the 'melting rate instantaneous' postprocessor must "
                             "each contain exactly three values."));

      std::vector<double> local_data;

      const double time_step = this->get_timestep();
      const double tau_seconds = tau_yr * year_in_seconds;

      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            // Cell average of div(V), computed on the same quadrature that the
            // advection assembly uses (see above).
            advection_fe_values.reinit(cell);
            std::vector<double> advection_divergences(n_advection_q_points);
            advection_fe_values[introspection.extractors.velocities]
            .get_function_divergences(this->get_solution(), advection_divergences);
            double divergence_u = 0.0;
            for (unsigned int q = 0; q < n_advection_q_points; ++q)
              divergence_u += advection_divergences[q] / n_advection_q_points;
            // Values and gradients, read from the in-memory solution on the
            // system finite element.
            std::vector<double> temperature(n_q_points);
            std::vector<double> fluid_pressure(n_q_points);
            fe_values[introspection.extractors.temperature]
            .get_function_values(this->get_solution(), temperature);
            fe_values[fluid_pressure_extractor]
            .get_function_values(this->get_solution(), fluid_pressure);

            std::vector<Tensor<1,dim>> velocity(n_q_points);
            std::vector<Tensor<1,dim>> fluid_velocity(n_q_points);
            std::vector<double> velocity_divergence(n_q_points);
            std::vector<double> fluid_velocity_divergence(n_q_points);
            fe_values[introspection.extractors.velocities]
            .get_function_values(this->get_solution(), velocity);
            fe_values[introspection.extractors.velocities]
            .get_function_divergences(this->get_solution(), velocity_divergence);
            fe_values[fluid_velocity_extractor]
            .get_function_values(this->get_solution(), fluid_velocity);
            fe_values[fluid_velocity_extractor]
            .get_function_divergences(this->get_solution(), fluid_velocity_divergence);

            std::vector<std::vector<double>> composition_values(n_compositional_fields,
                                                                std::vector<double>(n_q_points));
            std::vector<std::vector<Tensor<1,dim>>> composition_gradients(n_compositional_fields,
                std::vector<Tensor<1,dim>>(n_q_points));
            for (unsigned int c = 0; c < n_compositional_fields; ++c)
              {
                fe_values[introspection.extractors.compositional_fields[c]]
                .get_function_values(this->get_solution(), composition_values[c]);
                fe_values[introspection.extractors.compositional_fields[c]]
                .get_function_gradients(this->get_solution(), composition_gradients[c]);
              }

            // Rates of change of temperature and fluid pressure. These are not
            // available as a spatial flux divergence (the pressure in particular
            // is a Stokes unknown), so they are evaluated as a finite difference
            // between the current and the previous time step. The very first
            // time step has length zero, so the rates are left at zero there.
            std::vector<double> temperature_rate(n_q_points, 0.0);
            std::vector<double> pressure_rate(n_q_points, 0.0);
            if ((perturb_temperature || perturb_pressure) && time_step > 0.0)
              {
                std::vector<double> old_temperature(n_q_points);
                std::vector<double> old_fluid_pressure(n_q_points);
                fe_values[introspection.extractors.temperature]
                .get_function_values(this->get_old_solution(), old_temperature);
                fe_values[fluid_pressure_extractor]
                .get_function_values(this->get_old_solution(), old_fluid_pressure);

                for (unsigned int q = 0; q < n_q_points; ++q)
                  {
                    if (perturb_temperature)
                      temperature_rate[q] = (temperature[q] - old_temperature[q]) / time_step;
                    if (perturb_pressure)
                      pressure_rate[q] = (fluid_pressure[q] - old_fluid_pressure[q]) / time_step;
                  }
              }

            // Per-point state, rates of change and the liquid transport term.
            std::vector<std::array<double, n_chemical_components>> stored_bulk(n_q_points);
            std::vector<std::array<double, n_chemical_components>> bulk_rate(n_q_points);
            std::vector<std::array<double, n_chemical_components>> liquid_transport(n_q_points);
            std::vector<std::array<double, n_chemical_components>> liquid_transport_cell_averaged(n_q_points);
            // Diagnostics: y-component of grad(c_l) and phi (v_l - V).grad(c_l).
            // For the uniform-phase calibration case c_l is spatially constant,
            // so both should vanish; they are written to the output file.
            std::vector<std::array<double, n_chemical_components>> liquid_gradient_y(n_q_points);
            std::vector<double> q_dot_grad_c_l(n_q_points, 0.0);

            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const double porosity = std::max(0.0, std::min(1.0, composition_values[porosity_index][q]));
                const Tensor<1,dim> porosity_gradient = composition_gradients[porosity_index][q];

                for (unsigned int i = 0; i < n_chemical_components; ++i)
                  {
                    const double bulk   = composition_values[bulk_index[i]][q];
                    const double liquid = composition_values[liquid_index[i]][q];
                    const Tensor<1,dim> bulk_gradient   = composition_gradients[bulk_index[i]][q];
                    const Tensor<1,dim> liquid_gradient = composition_gradients[liquid_index[i]][q];

                    stored_bulk[q][i] = bulk;

                    // Liquid mass per unit volume theta = phi * c_l and its gradient.
                    const double theta = porosity * liquid;
                    const Tensor<1,dim> theta_gradient = porosity_gradient * liquid + porosity * liquid_gradient;

                    // Solid part of the bulk flux: (c_bulk - theta) * V.
                    const double solid_mass = bulk - theta;
                    const Tensor<1,dim> solid_mass_gradient = bulk_gradient - theta_gradient;

                    // Eulerian rate of change of the bulk concentration.
                    if (use_cell_averaged_form)
                      {
                        // Exactly the discrete equation solved by the
                        // advection system for the bulk concentrations
                        // (MeltAdvectionSystem::execute() in melt.cc):
                        //   d(c_bar)/dt = - V.grad(c_bar) - c_bar D + c_l D
                        //                 - phi (u_f - V).grad(c_l)
                        // with D the arithmetic cell average of div(V).
                        bulk_rate[q][i] = -( velocity[q] * bulk_gradient
                                             + bulk * divergence_u )
                                          + liquid * divergence_u
                                          - porosity * (fluid_velocity[q] - velocity[q]) * liquid_gradient;
                      }
                    else
                      {
                        // Pointwise flux divergence form,
                        //   d(c_bar)/dt = - div[ (c_bar - theta) V + theta v_l ].
                        bulk_rate[q][i] = -( solid_mass_gradient * velocity[q]
                                             + solid_mass * velocity_divergence[q]
                                             + theta_gradient * fluid_velocity[q]
                                             + theta * fluid_velocity_divergence[q] );
                      }

                    // Divergence of the liquid mass flux, i.e. the transport term.
                    liquid_transport[q][i] = theta_gradient * fluid_velocity[q]
                                             + theta * fluid_velocity_divergence[q];

                    // Same transport term, but with div(V) and div(q) replaced by
                    // the cell-averaged D and -D, exactly as the advection
                    // system does in the bulk equation. Starting from
                    //   div(theta v_l) = div(theta V) + div(c_l q)
                    // and using div(q) -> -D gives
                    //   div(theta v_l) = V.grad(theta) + theta D
                    //                    + q.grad(c_l) - c_l D.
                    // Summing over i reproduces the total in the porosity form
                    // (V.grad(phi) - (1-phi) D), so sum_i Gamma^i = Gamma holds
                    // identically.
                    liquid_transport_cell_averaged[q][i]
                      = velocity[q] * theta_gradient
                        + theta * divergence_u
                        + porosity * (fluid_velocity[q] - velocity[q]) * liquid_gradient
                        - liquid * divergence_u;

                    // Diagnostics (see the declaration above).
                    liquid_gradient_y[q][i] = liquid_gradient[dim-1];
                    q_dot_grad_c_l[q] += porosity * (fluid_velocity[q] - velocity[q]) * liquid_gradient;
                  }
              }

            // Evaluate the equilibrium map. It depends on the *fluid* pressure,
            // whereas MaterialModelInputs fills 'pressure' from the total
            // pressure variable, so overwrite it.
            MaterialModel::MaterialModelInputs<dim> base_inputs(fe_values, cell, introspection,
                                                                this->get_solution(), true);
            for (unsigned int q = 0; q < n_q_points; ++q)
              base_inputs.pressure[q] = fluid_pressure[q];

            auto evaluate_equilibrium
              = [&](const MaterialModel::MaterialModelInputs<dim> &inputs,
                    std::vector<double> &equilibrium_melt_fraction,
                    std::vector<std::array<double, n_chemical_components>> &equilibrium_liquid_composition)
            {
              MaterialModel::MaterialModelOutputs<dim> outputs(n_q_points, n_compositional_fields);
              this->get_material_model().create_additional_named_outputs(outputs);
              this->get_material_model().evaluate(inputs, outputs);

              const MaterialModel::PrescribedFieldOutputs<dim> *prescribed_fields =
                outputs.template get_additional_output<MaterialModel::PrescribedFieldOutputs<dim>>();
              AssertThrow(prescribed_fields != nullptr,
                          ExcMessage("The material model in use does not provide PrescribedFieldOutputs, "
                                     "which the 'melting rate instantaneous' postprocessor requires to "
                                     "read the equilibrium melt fraction and liquid composition."));

              equilibrium_melt_fraction.resize(n_q_points);
              equilibrium_liquid_composition.resize(n_q_points);
              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  equilibrium_melt_fraction[q] = prescribed_fields->prescribed_field_outputs[q][porosity_index];
                  for (unsigned int i = 0; i < n_chemical_components; ++i)
                    equilibrium_liquid_composition[q][i]
                      = prescribed_fields->prescribed_field_outputs[q][liquid_index[i]];
                }
            };

            // Build the perturbed state x + sign * tau * dx/dt. The bulk
            // composition is renormalized so that it sums to one, exactly as in
            // the external postprocessing implementation.
            auto build_perturbed_inputs = [&](const double sign)
            {
              MaterialModel::MaterialModelInputs<dim> perturbed_inputs = base_inputs;

              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  double sum = 0.0;
                  std::array<double, n_chemical_components> perturbed_bulk;
                  for (unsigned int i = 0; i < n_chemical_components; ++i)
                    {
                      double value = stored_bulk[q][i];
                      if (perturb_composition)
                        value += sign * tau_seconds * bulk_rate[q][i];
                      value = std::max(0.0, value);
                      perturbed_bulk[i] = value;
                      sum += value;
                    }

                  if (sum > 0.0)
                    for (unsigned int i = 0; i < n_chemical_components; ++i)
                      perturbed_inputs.composition[q][bulk_index[i]] = perturbed_bulk[i] / sum;

                  if (perturb_temperature)
                    perturbed_inputs.temperature[q] = temperature[q] + sign * tau_seconds * temperature_rate[q];

                  if (perturb_pressure)
                    perturbed_inputs.pressure[q] = fluid_pressure[q] + sign * tau_seconds * pressure_rate[q];
                }

              return perturbed_inputs;
            };

            std::vector<double> lower_melt_fraction;
            std::vector<double> upper_melt_fraction;
            std::vector<std::array<double, n_chemical_components>> lower_liquid_composition;
            std::vector<std::array<double, n_chemical_components>> upper_liquid_composition;

            double denominator = 1.0;
            if (use_central_difference)
              {
                const MaterialModel::MaterialModelInputs<dim> minus_inputs = build_perturbed_inputs(-1.0);
                const MaterialModel::MaterialModelInputs<dim> plus_inputs  = build_perturbed_inputs(+1.0);

                evaluate_equilibrium(minus_inputs, lower_melt_fraction, lower_liquid_composition);
                evaluate_equilibrium(plus_inputs,  upper_melt_fraction, upper_liquid_composition);
                denominator = 2.0;
              }
            else
              {
                const MaterialModel::MaterialModelInputs<dim> plus_inputs = build_perturbed_inputs(+1.0);

                evaluate_equilibrium(base_inputs, lower_melt_fraction, lower_liquid_composition);
                evaluate_equilibrium(plus_inputs, upper_melt_fraction, upper_liquid_composition);
              }

            // Gamma^i/rho_l = [theta(t+tau) - theta(t)] / (denominator*tau)
            //                 + div(phi c_l^i v_l).
            std::vector<double> gamma(n_chemical_components);
            for (unsigned int q = 0; q < n_q_points; ++q)
              {
                const bool no_melt = (lower_melt_fraction[q] <= melt_fraction_threshold)
                                     && (upper_melt_fraction[q] <= melt_fraction_threshold);

                double total = 0.0;
                double response_total = 0.0;
                for (unsigned int i = 0; i < n_chemical_components; ++i)
                  {
                    const double theta_upper = upper_melt_fraction[q] * upper_liquid_composition[q][i];
                    const double theta_lower = lower_melt_fraction[q] * lower_liquid_composition[q][i];

                    const double response = (theta_upper - theta_lower) / (denominator * tau_seconds);
                    const double transport = (use_cell_averaged_form
                                              ? liquid_transport_cell_averaged[q][i]
                                              : liquid_transport[q][i]);
                    double value = response + transport;

                    if (zero_where_no_melt && no_melt)
                      value = 0.0;

                    gamma[i] = value;
                    total += value;
                    response_total += response;
                  }

                // Same discrete equation, written in the solid (porosity) form
                // that the model's own 'melting_rate' field uses:
                //   Gamma/rho = d(phi)/dt + V.grad(phi) - (1-phi) D.
                const double porosity_q
                  = std::max(0.0, std::min(1.0, composition_values[porosity_index][q]));
                const double cell_averaged_total
                  = response_total
                    + velocity[q] * composition_gradients[porosity_index][q]
                    - (1.0 - porosity_q) * divergence_u;

                // Tie-line coordinates recovered from each component, and their
                // spread (0 if the state lies exactly on the tie line).
                std::array<double, n_chemical_components> tie_line_f;
                double tie_line_min = std::numeric_limits<double>::max();
                double tie_line_max = std::numeric_limits<double>::lowest();
                for (unsigned int i = 0; i < n_chemical_components; ++i)
                  {
                    const double denominator = tie_line_liquid[i] - tie_line_solid[i];
                    tie_line_f[i] = (stored_bulk[q][i] - tie_line_solid[i]) / denominator;
                    tie_line_min = std::min(tie_line_min, tie_line_f[i]);
                    tie_line_max = std::max(tie_line_max, tie_line_f[i]);
                  }

                const Point<dim> &p = fe_values.quadrature_point(q);
                local_data.push_back(p[0]);
                local_data.push_back(p[dim-1]);
                local_data.push_back(total);
                for (unsigned int i = 0; i < n_chemical_components; ++i)
                  local_data.push_back(gamma[i]);
                local_data.push_back(cell_averaged_total);
                for (unsigned int i = 0; i < n_chemical_components; ++i)
                  local_data.push_back(liquid_gradient_y[q][i]);
                local_data.push_back(q_dot_grad_c_l[q]);
                local_data.push_back(tie_line_max - tie_line_min);
                for (unsigned int i = 0; i < n_chemical_components; ++i)
                  local_data.push_back(tie_line_f[i]);
              }
          }

      // ------------------------------------------------------------------
      // Collect everything on rank zero and write it out in double precision.
      // ------------------------------------------------------------------
      const std::vector<std::vector<double>> gathered_data
        = Utilities::MPI::all_gather(this->get_mpi_communicator(), local_data);
      std::vector<double> global_data;
      for (const auto &rank_data : gathered_data)
        global_data.insert(global_data.end(), rank_data.begin(), rank_data.end());

      double local_max = 0.0;
      const unsigned int output_stride = 2 + n_output_scalars;
      for (unsigned int i = 0; i < local_data.size(); i += output_stride)
        local_max = std::max(local_max, std::abs(local_data[i+2]));
      const double global_max = Utilities::MPI::max(local_max, this->get_mpi_communicator());

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) == 0)
        {
          const std::string directory = this->get_output_directory() + "melting_rate_instantaneous";
          std::filesystem::create_directories(directory);

          std::ostringstream filename;
          filename << directory << "/step_"
                   << std::setw(6) << std::setfill('0') << this->get_timestep_number()
                   << ".csv";

          std::ofstream output(filename.str());
          output << std::setprecision(17);
          output << "x,y,melting_rate_instantaneous";
          for (const auto &component : chemical_component_names)
            output << "," << component << "_melting_rate_instantaneous";
          output << ",melting_rate_cell_averaged";
          for (const auto &component : chemical_component_names)
            output << ",d_c_l_dy_" << component;
          output << ",q_dot_grad_c_l";
          output << ",tie_line_f_spread";
          for (const auto &component : chemical_component_names)
            output << ",tie_line_f_" << component;
          output << "\n";

          const unsigned int stride = output_stride;
          for (unsigned int i = 0; i < global_data.size(); i += stride)
            {
              for (unsigned int j = 0; j < stride; ++j)
                output << (j == 0 ? "" : ",") << global_data[i+j];
              output << "\n";
            }
        }

      statistics.add_value ("Maximal melting rate (1/s)", global_max);
      statistics.set_precision ("Maximal melting rate (1/s)", 8);
      statistics.set_scientific ("Maximal melting rate (1/s)", true);

      std::ostringstream description;
      description.precision(4);
      description << "1/s" << " (tau = " << tau_yr << " yr";

      return std::pair<std::string, std::string> ("Maximal melting rate:",
                                                  description.str());
    }



    template <int dim>
    void
    MeltingRateInstantaneous<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Melting rate instantaneous");
        {
          prm.declare_entry("Virtual time step", "1.0",
                            Patterns::Double(0.0),
                            "The virtual time step tau [yr] used to construct the small convective "
                            "perturbation that breaks the instantaneous equilibrium. It is unrelated "
                            "to the time step of the simulation: it only sets the size of the "
                            "finite-difference perturbation of the equilibrium map.");

          prm.declare_entry("Perturb composition", "true",
                            Patterns::Bool(),
                            "Whether the bulk composition is perturbed along its advective rate of "
                            "change, computed from the model's own bulk mass flux.");

          prm.declare_entry("Perturb temperature", "false",
                            Patterns::Bool(),
                            "Whether the temperature is perturbed along its rate of change, "
                            "evaluated as a finite difference between the current and the previous "
                            "time step. Disabled by default: this one-step finite difference is not "
                            "tied to any discrete operator of the model and has no counterpart in "
                            "the liquid transport term, so it can inject structure that is not a "
                            "physical melting rate. Kept for future validation.");

          prm.declare_entry("Perturb pressure", "false",
                            Patterns::Bool(),
                            "Whether the fluid pressure is perturbed along its rate of change, "
                            "evaluated as a finite difference between the current and the previous "
                            "time step. Disabled by default. The fluid pressure is a Stokes unknown, "
                            "so its rate of change is only available as a one-step difference; for "
                            "cases with a pressure-dependent melting curve (A != 0) this term is "
                            "large and concentrated where the porosity gradient is steepest, i.e. at "
                            "the flanks of a porosity wave, and it has no counterpart in the liquid "
                            "transport term. Kept for future validation.");

          prm.declare_entry("Use central difference", "false",
                            Patterns::Bool(),
                            "If true, the equilibrium map is evaluated at both +tau and -tau, which "
                            "is second-order accurate in tau. If false, a forward difference from "
                            "the current state is used, which matches the external postprocessing "
                            "implementation. In both cases two equilibrium solves are required.");

          prm.declare_entry("Use cell averaged form", "true",
                            Patterns::Bool(),
                            "If true, the rate of change of the bulk composition used to perturb the "
                            "equilibrium is evaluated in the same cell-averaged discrete form that "
                            "the advection system solves (arithmetic cell average of div(V) plus the "
                            "c_l*div(V) - phi*(u_f-V).grad(c_l) source of melt.cc), and the total "
                            "melting rate is written as a second column in the porosity form "
                            "d(phi)/dt + V.grad(phi) - (1-phi) D. If false, the pointwise flux "
                            "divergence is used instead.");

          prm.declare_entry("Zero where no melt", "true",
                            Patterns::Bool(),
                            "Whether the melting rate is set to zero at points where the equilibrium "
                            "melt fraction is essentially zero.");

          prm.declare_entry("Melt fraction threshold", "1e-12",
                            Patterns::Double(0.0),
                            "The melt fraction below which the melting rate is set to zero if "
                            "'Zero where no melt' is enabled.");

          prm.declare_entry("Tie line solid composition", "0.70013835799, 0.29933620004, 5.2544197096e-04",
                            Patterns::List(Patterns::Double()),
                            "The three solid end-member concentrations S_i of the tie line used by "
                            "the 'uniform phase' calibration case. Together with the liquid "
                            "end members they define the tie-line-consistency diagnostic "
                            "f_i = (c_bar_i - S_i)/(L_i - S_i) that is written to the output file.");

          prm.declare_entry("Tie line liquid composition", "0.0084701555582, 0.2850051987, 0.70652464574",
                            Patterns::List(Patterns::Double()),
                            "The three liquid end-member concentrations L_i of the tie line used by "
                            "the 'uniform phase' calibration case, see 'Tie line solid composition'.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }



    template <int dim>
    void
    MeltingRateInstantaneous<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Postprocess");
      {
        prm.enter_subsection("Melting rate instantaneous");
        {
          tau_yr                  = prm.get_double("Virtual time step");
          perturb_composition     = prm.get_bool("Perturb composition");
          perturb_temperature     = prm.get_bool("Perturb temperature");
          perturb_pressure        = prm.get_bool("Perturb pressure");
          use_central_difference  = prm.get_bool("Use central difference");
          use_cell_averaged_form      = prm.get_bool("Use cell averaged form");
          zero_where_no_melt      = prm.get_bool("Zero where no melt");
          melt_fraction_threshold = prm.get_double("Melt fraction threshold");

          tie_line_solid  = Utilities::string_to_double(
                              Utilities::split_string_list(prm.get("Tie line solid composition")));
          tie_line_liquid = Utilities::string_to_double(
                              Utilities::split_string_list(prm.get("Tie line liquid composition")));

          AssertThrow(tau_yr > 0.0,
                      ExcMessage("The 'Virtual time step' of the 'melting rate instantaneous' "
                                 "postprocessor must be positive."));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(MeltingRateInstantaneous,
                                  "melting rate instantaneous",
                                  "A postprocessor that computes the melting rate "
                                  "$\\Gamma/\\rho_l$ and the three component-wise values "
                                  "$\\Gamma^i/\\rho_l$ by perturbing the instantaneous thermodynamic "
                                  "equilibrium with a small convective perturbation. It uses the "
                                  "in-memory, double precision solution and the material model's own "
                                  "equilibrium solver, so that the float32 output round trip of the "
                                  "external postprocessing workflow is avoided and no additional "
                                  "advection equation has to be solved. The result is written to "
                                  "'<output directory>/melting_rate_instantaneous/' in double "
                                  "precision. Physical units: 1/s.")
  }
}
