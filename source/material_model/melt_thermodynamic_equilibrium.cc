/*
  Copyright (C) 2015 - 2024 by the authors of the ASPECT code.

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


#include <aspect/material_model/melt_thermodynamic_equilibrium.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/utilities.h>

#include <deal.II/base/signaling_nan.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/numerics/fe_field_function.h>

#  pragma message("Compiling melt_thermodynamic_equilibrium.cc")

// to use pcout for debugging output
#include <aspect/simulator.h>
namespace
{
  const bool local_debug = true;
}

namespace aspect
{
  namespace MaterialModel
  {
    // i m not sure if we need this
    using namespace dealii;

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    reference_darcy_coefficient () const
    {
      // 0.01 = 1% melt
      return reference_permeability * Utilities::fixed_power<3>(0.01) / eta_f;
    }

    template <int dim>
    bool
    MeltThermodynamicEquilibrium<dim>::
    is_compressible () const
    {
      return false;
    }

    // match component index to composition indices
    template <int dim>
    bool
    MeltThermodynamicEquilibrium<dim>::
    match_component_index_to_composition_indices ()
    {
      // if equilibrium calculation is disabled, return true directly
      if (!enable_equilibrium_calculation)
        return true;

      component_index_to_composition_indices.clear();
      for (unsigned int comp_idx = 0; comp_idx < n_components; ++comp_idx)
        {
          std::string solid_name = component_names[comp_idx] + component_name_solid_suffix;
          std::string liquid_name = component_names[comp_idx] + component_name_liquid_suffix;
          bool solid_found = false;
          bool liquid_found = false;
          unsigned int solid_comp_field_idx = numbers::invalid_unsigned_int;
          unsigned int liquid_comp_field_idx = numbers::invalid_unsigned_int;

          // for (unsigned int field_idx = 0; field_idx < this->introspection().n_compositional_fields; ++field_idx)
          //   {
          //     std::string current_name = this->introspection().name_for_compositional_index(field_idx);
          //     if (current_name == solid_name)
          //       {
          //         solid_found = true;
          //         solid_comp_field_idx = field_idx;
          //       }
          //     if (current_name == liquid_name)
          //       {
          //         liquid_found = true;
          //         liquid_comp_field_idx = field_idx;
          //       }
          //   }

          solid_comp_field_idx = this->introspection().compositional_index_for_name(solid_name);
          liquid_comp_field_idx = this->introspection().compositional_index_for_name(liquid_name);
          solid_found = (solid_comp_field_idx != numbers::invalid_unsigned_int);
          liquid_found = (liquid_comp_field_idx != numbers::invalid_unsigned_int);

          if (solid_found && liquid_found)
            {
              component_index_to_composition_indices.push_back(std::make_pair(solid_comp_field_idx, liquid_comp_field_idx));
            }
          else
            {
              return false;
            }
        }
      return true;
    }

    // for convenience, we define functions
    // for melting curves and equilibrium constants

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    temperature_melting (const double pressure,
                         const double temperature_m_0,
                         const double coefficient_A,
                         const double coefficient_B) const
    {
      // T_m = T_m_0 + A * P + B * P^2
      return temperature_m_0
             + coefficient_A * pressure
             + coefficient_B * pressure * pressure;
    }
    
    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    equilibrium_constant (/* const double pressure, */
                          const double temperature,
                          const double latent_heat,
                          const double tuning_parameter,
                          const double _T_m) const
    {
      // K = exp((L / r) * (1 / T - 1 / T_m))
      double exponent = (latent_heat / tuning_parameter)
                        * (1.0 / temperature
                           - 1.0 / _T_m);
      // to avoid overflow of exp function
      if (exponent > 700.0)
        return std::exp(700.0);
      else if (exponent < -700.0)
        return std::exp(-700.0);
      else
      return std::exp(exponent);
    }

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    find_solidus (const std::vector<double>& melting_points,
                  const std::vector<double>& bulk_concentrations,
                  const std::vector<double>& latent_heats,
                  const std::vector<double>& tuning_parameters) const
    {
      bool size_is_matching = (melting_points.size() == bulk_concentrations.size()
                          && melting_points.size() == latent_heats.size()
                          && melting_points.size() == tuning_parameters.size());
      AssertThrow(size_is_matching,
        ExcMessage("The sizes of the input vectors do not match."
                   "while calculating the solidus temperature."));
      auto solidus_equation = [&] (const double temp_solidus)
        {
          std::vector<double> eq_consts(melting_points.size());
          unsigned int comp_idx = 0;
          for (comp_idx = 0; comp_idx < melting_points.size(); ++comp_idx)
            {
              eq_consts[comp_idx] = equilibrium_constant(/* const double pressure, */
                                                  temp_solidus,
                                                  latent_heats[comp_idx],
                                                  tuning_parameters[comp_idx],
                                                  melting_points[comp_idx]);
            }
          double sum = 0.0;
          for (comp_idx = 0; comp_idx < melting_points.size(); ++comp_idx) 
            {
              sum += bulk_concentrations[comp_idx] / eq_consts[comp_idx];
            }
          return sum - 1;
        };
      
      // now we are going to find the root of the solidus_equation
      // we can use the bisection method to find the root
      // we need to guess the upper and lower bounds of the root
      // to make sure we contain the root, make the interval larger
      // on both sides with a shift
      const double shift = 100.0;
      const double lower_bound = *std::min_element(melting_points.begin(), melting_points.end()) - shift;
      const double upper_bound = *std::max_element(melting_points.begin(), melting_points.end()) + shift;
      double solidus = 0.0;
      try
        {
          solidus = bisection(solidus_equation,
                              lower_bound,
                              upper_bound);
        }
      catch(const std::exception& e)
        {
          // AssertThrow(false,
          //   ExcMessage("Error in finding the solidus temperature: "
          //              + std::string(e.what())));
          // return dealii::numbers::signaling_nan<double>();
          AssertThrow(false,ExcMessage("")); // enpty message as a placeholder
        }
      return solidus;

    }

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    find_liquidus (const std::vector<double>& melting_points,
                   const std::vector<double>& bulk_concentrations,
                   const std::vector<double>& latent_heats,
                   const std::vector<double>& tuning_parameters) const
    {
      bool size_is_matching = (melting_points.size() == bulk_concentrations.size()
                          && melting_points.size() == latent_heats.size()
                          && melting_points.size() == tuning_parameters.size());
      AssertThrow(size_is_matching,
        ExcMessage("The sizes of the input vectors do not match."
                   "while calculating the liquidus temperature."));
      auto liquidus_equation = [&] (const double temp_liquidus)
        {
          std::vector<double> eq_consts(melting_points.size());
          unsigned int comp_idx = 0;
          for (comp_idx = 0; comp_idx < melting_points.size(); ++comp_idx)
            {
              eq_consts[comp_idx] = equilibrium_constant(/* const double pressure, */
                                                  temp_liquidus,
                                                  latent_heats[comp_idx],
                                                  tuning_parameters[comp_idx],
                                                  melting_points[comp_idx]);
            }
          double sum = 0.0;
          for (comp_idx = 0; comp_idx < melting_points.size(); ++comp_idx) 
            {
              sum += bulk_concentrations[comp_idx] * eq_consts[comp_idx];
            }
          return sum - 1;
        };
      
      // now we are going to find the root of the liquidus_equation
      // we can use the bisection method to find the root
      // we need to guess the upper and lower bounds of the root
      // to make sure we contain the root, make the interval larger
      // on both sides with a shift
      const double shift = 100.0;
      const double lower_bound = *std::min_element(melting_points.begin(), melting_points.end()) - shift;
      const double upper_bound = *std::max_element(melting_points.begin(), melting_points.end()) + shift;
      double liquidus = 0.0;
      try
        {
          liquidus = bisection(liquidus_equation,
                              lower_bound,
                              upper_bound);
        }
      catch(const std::exception& e)
        {
          // AssertThrow(false,
          //   ExcMessage("Error in finding the liquidus temperature: "
          //              + std::string(e.what())));
          // return dealii::numbers::signaling_nan<double>();
          AssertThrow(false,ExcMessage("")); // enpty message as a placeholder
        }
      return liquidus;
    }

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    solve_eq_melt_fraction (double temperature, // to minus ZERO_CELSIUS_IN_KELVIN
                   const double pressure,
                   std::vector<double> bulk_concentrations) const
    {
      // return 0.001; // temporary return value
      // if (local_debug) 
      //   {
      //     std::string concentration_string = "";
      //     for (unsigned int i=0; i<bulk_concentrations.size(); ++i)
      //       {
      //         concentration_string += std::to_string(bulk_concentrations[i]);
      //         if (i != bulk_concentrations.size()-1)
      //           concentration_string += ", ";
      //       }
      //     this->get_pcout() << "[MM eval] solve_eq_melt_fraction called with T="
      //                       << temperature << " K, P=" << pressure << " Pa "
      //                       << "and concentrations: " << concentration_string << "\n";
      //   }
      
      // check if the bulk concentrations are valid
      // if not, try to fix them and give warning
      for (unsigned int _i = 0; _i < bulk_concentrations.size(); ++_i)
        {
          if (bulk_concentrations[_i] < 0.0)
            {
              // this->get_pcout() << "[MM eval] Warning: negative bulk concentration for component "
              //                   << _i << ": " << bulk_concentrations[_i] << "at P=" << pressure << " Pa"
              //                   << ". Setting it to zero.\n";
              bulk_concentrations[_i] = 0.0;
            }
          if (bulk_concentrations[_i] > 1.0)
            {
              // this->get_pcout() << "[MM eval] Warning: bulk concentration greater than one for component "
              //                   << _i << ": " << bulk_concentrations[_i] << "at P=" << pressure << " Pa"
              //                   << ". Setting it to one.\n";
              bulk_concentrations[_i] = 1.0;
            }
        }


      // Keller and Katz (2016) do thermodynamic equilibrium calculations
      // with temperature in Celsius. We follow their approach here.
      // the temperature field is in Kelvin in ASPECT

      // we've convert temperature into Celsius yet
      // input parameter is required to be in Celsius

      double inner_melt_fraction = 0.0;
      std::vector<double> melting_points;
      melting_points.reserve(bulk_concentrations.size());
      for(unsigned int _i = 0; _i < bulk_concentrations.size(); ++_i)
        {
          melting_points.push_back(temperature_melting(pressure,
                                                       melting_point_0_values[_i],
                                                       melting_curve_coefficient_A_values[_i],
                                                       melting_curve_coefficient_B_values[_i]));
        }
      std::vector<double> equilibrium_constants;
      equilibrium_constants.reserve(bulk_concentrations.size());
      for(unsigned int _i = 0; _i < bulk_concentrations.size(); ++_i)
        {
          equilibrium_constants.push_back(equilibrium_constant(temperature,
                                                                latent_heat_values[_i],
                                                                tuning_parameter_values[_i],
                                                                melting_points[_i]));
        }
      std::vector<double> latent_heats = latent_heat_values;
      std::vector<double> tuning_parameters = tuning_parameter_values;

      // now we are going to describe an equation about the melt fraction "f"
      // and find its root as the melt fraction
      // we gonna difine a lambda function and use the methods in deal.ii
      // it's a pity that i can't find how to use the root finding method in deal.ii
      // so let's just implement the bisection method with fundamental funcitions
      // in std library

      // first we need to define some lambda functions
      auto f_eq_equation = [&](const double f_eq)
      {
        // this is the function whose root we want to find
        double c_solid_sum = 0.0;
        for(unsigned int _i = 0; _i < bulk_concentrations.size(); ++_i)
          {
            c_solid_sum += calculate_concentration_solid(bulk_concentrations[_i],
                                                        f_eq,
                                                        equilibrium_constants[_i]);
          }
        double c_liquid_sum = 0.0;
        for(unsigned int _i = 0; _i < bulk_concentrations.size(); ++_i)
          {
            c_liquid_sum += calculate_concentration_liquid(bulk_concentrations[_i],
                                                          f_eq,
                                                          equilibrium_constants[_i]);
          }
        return c_solid_sum - c_liquid_sum;
      };

      // before we actually solving the equation
      // we'd better check if the temperature is between the solidus and liquidus
      // we need to define the value vector which is need
      // by the find_solidus function and find_liquidus function
      const std::vector<double> _melting_points = melting_points;
      const std::vector<double> _bulk_concentrations = bulk_concentrations;
      const std::vector<double> _latent_heats = latent_heats;
      const std::vector<double> _tuning_parameters = tuning_parameters;
      double solidus = 0.0;
      double liquidus = 0.0;
      try
      {
        solidus = find_solidus(_melting_points,
                                _bulk_concentrations,
                                _latent_heats,
                                _tuning_parameters);
      }
      catch(const std::exception& e)
      {
        std::string concentration_string = "";
        for (unsigned int i=0; i<_bulk_concentrations.size(); ++i)
          {
            concentration_string += std::to_string(_bulk_concentrations[i]);
            if (i != _bulk_concentrations.size()-1)
              concentration_string += ", ";
          }
        AssertThrow(false,
          ExcMessage("Error in finding the solidus temperature: "
                     "at pressure " + std::to_string(pressure) + "\n"
                      + "and with concentrations " + concentration_string
                      + std::string(e.what())));
      }
      try
      {
        liquidus = find_liquidus(_melting_points,
                                  _bulk_concentrations,
                                  _latent_heats,
                                  _tuning_parameters);
      }
      catch(const std::exception& e)
      {
        std::string concentration_string = "";
        for (unsigned int i=0; i<_bulk_concentrations.size(); ++i)
          {
            concentration_string += std::to_string(_bulk_concentrations[i]);
            if (i != _bulk_concentrations.size()-1)
              concentration_string += ", ";
          }
        AssertThrow(false,
          ExcMessage("Error in finding the liquidus temperature: "
                     "at pressure " + std::to_string(pressure) + "\n"
                      + "and with concentrations " + concentration_string
                      + std::string(e.what())));
      }

      // there is another way to check whether the state is between the solidus and liquidus
      // we can calculate the value of f_eq_equation at f_eq = 0.0 and f_eq = 1.0
      // if both values are positive, the state is all solid
      // if both values are negative, the state is all liquid
      // if one value is positive and the other is negative, the state is between the solidus and liquidus
      // we are not going to implement this method yet because it has not been tested
      // we need to check if f_eq_equation is a monotonically increasing function

      // now we can check if the temperature is between the solidus and liquidus
      const bool is_all_solid = (temperature < solidus);
      const bool is_all_liquid = (temperature > liquidus);
      if (is_all_solid)
        {
          // this is a solid state
          return 0.0;
        }
      else if (is_all_liquid)
        {
          // this is a liquid state
          return 1.0;
        }

      // we need to guess the upper and lower bounds of the root
      const double shift = 1e-6;
      const double lower_bound = 0.0 - shift;
      const double upper_bound = 1.0 + shift;
      inner_melt_fraction = 0.0;
      try
      {
        inner_melt_fraction = bisection(f_eq_equation,
                                        lower_bound,
                                        upper_bound);
      }
      catch(const std::exception& e)
      {
        AssertThrow(false,
          ExcMessage("Error in finding the melt fraction: "
                     + std::string(e.what())
                    + "\n Temperature: " + std::to_string(temperature)
                    + ", solidus: " + std::to_string(solidus)
                    + ", liquidus: " + std::to_string(liquidus)));
        return dealii::numbers::signaling_nan<double>();
      }
      
      // we need to check if the root is between 0 and 1
      AssertThrow(inner_melt_fraction >= 0.0 && inner_melt_fraction <= 1.0,
        ExcMessage("The root finder find the melt fraction."
                   "However The melt fraction is not between 0 and 1."));

      return inner_melt_fraction;
    }

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    calculate_concentration_solid (const double c_bulk,
                                    const double f, // melt fraction
                                    const double eq_const) const
    {
      return c_bulk / ((f / eq_const) + (1 - f));
    }
    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    calculate_concentration_liquid (const double c_bulk,
                                   const double f, // melt fraction
                                   const double eq_const) const
    {
      return c_bulk / (f + (1 - f) * eq_const);
    }

    // However, there have to be a function named melt_fractions
    // because the base class has this function declared as a pure virtual function
    // and we have to implement it
    template <int dim>
    void
    MeltThermodynamicEquilibrium<dim>::
    melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                    std::vector<double> &melt_fractions) const
    {
      std::vector<double> old_porosity(in.n_evaluation_points());
      // we want to get the porosity field from the old solution here,
      // because we need a field that is not updated in the nonlinear iterations
      if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid
          && this->get_timestep_number() > 0 && !this->get_parameters().use_operator_splitting)
        {
          // Prepare the field function

          Functions::FEFieldFunction<dim, LinearAlgebra::BlockVector>

          fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          fe_value.set_active_cell(in.current_cell);
          fe_value.value_list(in.position,
                              old_porosity,
                              this->introspection().component_indices.compositional_fields[porosity_idx]);
        }
      
      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          if (this->get_parameters().use_operator_splitting)
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              old_porosity[q] = in.composition[q][porosity_idx];
            }
 
          if (this->include_melt_transport())
            {
              // we need to pair the solid and liquid fields
              std::vector<unsigned int> solid_indices = {};
              std::vector<unsigned int> liquid_indices = {};
              for (unsigned int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                {
                  solid_indices.push_back(component_index_to_composition_indices[comp_idx].first);
                  liquid_indices.push_back(component_index_to_composition_indices[comp_idx].second);
                }
              // check if the length of solid_indices and liquid_indices equals to n_components
              AssertThrow(solid_indices.size() == n_components && liquid_indices.size() == n_components,
                          ExcMessage("The number of solid and liquid components do not match the expected number of components."));
              
              // we temporarily don't distinguish between melt fraction (mass) and porosity (volume)
              const double old_melt_fraction = old_porosity[q];
              std::vector<double> bulk_concentrations(n_components);
              for (unsigned int _i = 0; _i < n_components; ++_i)
                {
                  bulk_concentrations[_i] = (1 - old_melt_fraction) * in.composition[q][solid_indices[_i]]
                                            + old_melt_fraction * in.composition[q][liquid_indices[_i]];
                }
              melt_fractions[q] = this->solve_eq_melt_fraction(in.temperature[q],
                                                               std::max(0.0, in.pressure[q]),
                                                               bulk_concentrations);
            }
          else
            {
              melt_fractions[q] = 0.0;
            }
        }
    }

    template <int dim>
    void
    MeltThermodynamicEquilibrium<dim>::
    evaluate(const typename Interface<dim>::MaterialModelInputs &in, typename Interface<dim>::MaterialModelOutputs &out) const
    {
      // if (local_debug)
      // {
      //   const auto props = in.requested_properties;
      //   const bool on_cell = (in.current_cell.state() == IteratorState::valid);
      //   const bool need_sr = !in.strain_rate.empty(); // was compute_strain_rate=true?
      //   this->get_pcout() << "[MM eval] tstep=" << this->get_timestep_number()
      //               << " on_cell=" << on_cell
      //               << " need_sr=" << need_sr
      //               << " props=("
      //               << props
      //               << ")"
      //               << std::endl;
      // }
      const unsigned int num_compositional_fields = this->introspection().n_compositional_fields;
      std::vector<double> old_porosity(in.n_evaluation_points());
      std::vector<std::vector<double>> old_fields(in.n_evaluation_points(),
                                                  std::vector<double>(num_compositional_fields));

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();
      PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>();

      // we want to get the porosity field from the old solution here,
      // because we need a field that is not updated in the nonlinear iterations
      if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid
          && this->get_timestep_number() >= 0 && !this->get_parameters().use_operator_splitting)
        {
          // Prepare the field function

          Functions::FEFieldFunction<dim, LinearAlgebra::BlockVector>

          fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          fe_value.set_active_cell(in.current_cell);
          // // NOTE: if the extraction sentence for porosity is put outside the loop,
          // // the porosity values will be failed to be extracted and keep as zero.
          // // The reason is still unknown.
          // fe_value.value_list(in.position,
          //                     old_porosity,
          //                     this->introspection().component_indices.compositional_fields[porosity_idx]);
          for(unsigned int c=0; c<this->introspection().n_compositional_fields; ++c)
            {
              std::vector<double> temp_field(in.n_evaluation_points());
              // if (c == porosity_idx) continue;
              // fe_value.value_list(in.position,
              //                     temp_field,
              //                     this->introspection().component_indices.compositional_fields[c]);
              for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
                {
                  if (c == porosity_idx)
                    old_porosity[i] = temp_field[i];
                  old_fields[i][c] = temp_field[i];
                }
            }

          // // temporary fix: use current_linearization_point as old_field
          // old_fields = in.composition;
          // for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          // {
          //   const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
          //   old_porosity[i] = in.composition[i][porosity_idx];
          // }
          
        }
      else if (this->get_parameters().use_operator_splitting)
        for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
          {
            const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
            old_porosity[i] = in.composition[i][porosity_idx];
          }

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          // calculate density first, we need it for the reaction term
          // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
          double temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                      * thermal_expansivity;
          else
            temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

          // keep density constant for now
          out.densities[i] = reference_rho_s;

          out.viscosities[i] = eta_0;
          // By default, no melting or freezing --> set all reactions to zero
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            {
              out.reaction_terms[i][c] = 0.0;

              if (this->get_parameters().use_operator_splitting && reaction_rate_out != nullptr)
                reaction_rate_out->reaction_rates[i][c] = 0.0;
            }

          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));

              // calculate viscosity based on local melt
              out.viscosities[i] *= std::exp(- alpha_phi * porosity);

              if (include_melting_and_freezing && (in.requests_property(MaterialProperties::reaction_terms) ||
                                                   in.requests_property(MaterialProperties::reaction_rates) ||
                                                   in.requests_property(MaterialProperties::viscosity)))
                {
                  Assert(this->get_timestep_number()<=1 || std::isfinite(in.strain_rate[i].norm()),
                         ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                    "not filled by the caller."));
                  const double trace_strain_rate =
                    (this->get_timestep_number() > 1) ?
                    (trace(in.strain_rate[i]))
                    :
                    numbers::signaling_nan<double>();

                  // Get all chemical field indices and group them
                  // TODO: we'd better check if the matching is successful during setup
                  // or there will be segmentation fault while accessing invalid indices here
                  std::vector<unsigned int> solid_indices = {};
                  std::vector<unsigned int> liquid_indices = {};
                  if (enable_equilibrium_calculation)
                    {
                      for (unsigned int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                        {
                          solid_indices.push_back(component_index_to_composition_indices[comp_idx].first);
                          liquid_indices.push_back(component_index_to_composition_indices[comp_idx].second);
                        }
                      // if (local_debug)
                      //   {
                      //     std::string solid_idx_str = "";
                      //     std::string liquid_idx_str = "";
                      //     for (unsigned int _i = 0; _i < solid_indices.size(); ++_i)
                      //       {
                      //         solid_idx_str += std::to_string(solid_indices[_i]);
                      //         liquid_idx_str += std::to_string(liquid_indices[_i]);
                      //         if (_i != solid_indices.size() - 1)
                      //           {
                      //             solid_idx_str += ", ";
                      //             liquid_idx_str += ", ";
                      //           }
                      //       }
                      //     this->get_pcout() << "[MM eval] solid indices: " << solid_idx_str << "\n";
                      //     this->get_pcout() << "[MM eval] liquid indices: " << liquid_idx_str << "\n";
                      //   }
                      // check if the length of solid_indices and liquid_indices equals to n_components
                      AssertThrow(solid_indices.size() == n_components && liquid_indices.size() == n_components,
                                  ExcMessage("The number of solid and liquid components do not match the expected number of components."));
                    }

                  // we temporarily don't distinguish between melt fraction (mass) and porosity (volume)
                  const double old_melt_fraction = old_porosity[i];
                  std::vector<double> bulk_concentrations(n_components);
                  if (enable_equilibrium_calculation)
                    {
                    for (unsigned int _i = 0; _i < n_components; ++_i)
                      {
                        bulk_concentrations[_i] = (1 - old_melt_fraction) * old_fields[i][solid_indices[_i]]
                                                  + old_melt_fraction * old_fields[i][liquid_indices[_i]];
                        // if (local_debug)
                        //   {
                        //     this->get_pcout() << "[MM eval] bulk concentration of component "
                        //                       << _i << ": " << bulk_concentrations[_i] << "\n";
                        //   }
                      }
                    }
                  // Keller and Katz (2016) do thermodynamic equilibrium calculations
                  // with temperature in Celsius. We follow their approach here.
                  const double temperature_for_equilibrium_calculation = in.temperature[i] - ZERO_CELSIUS_IN_KELVIN;
                  const double eq_melt_fraction = (enable_equilibrium_calculation ? 
                                                   this->solve_eq_melt_fraction(temperature_for_equilibrium_calculation,
                                                                                std::max(0.0, in.pressure[i]),
                                                                                bulk_concentrations) : 
                                                   old_melt_fraction);
                  double porosity_change = 0.0;
                  if (prescribed_field_out != nullptr)
                    {
                      prescribed_field_out->prescribed_field_outputs[i][0] = std::min(background_porosity + eq_melt_fraction, 1.0);

                      // // check initial equilibrium at timestep 0
                      // if (this->get_timestep_number() == 0 && local_debug)
                      //   {
                      //     if (std::abs(eq_melt_fraction - old_porosity[i]) > 1e-8)
                      //       {
                      //         this->get_pcout() << "[MM eval] Initial porosity field at point "
                      //                           << in.position[i]
                      //                           << " does not match the equilibrium melt fraction. "
                      //                           << "Initial porosity: " << old_porosity[i]
                      //                           << ", equilibrium melt fraction: " << eq_melt_fraction
                      //                           << "\n";
                      //       }
                      //   }
                    }
                  porosity_change = std::min(background_porosity + eq_melt_fraction, 1.0) - old_porosity[i];
                  // do not allow negative porosity
                  if (old_porosity[i] + porosity_change < 0.0)
                    porosity_change = - old_porosity[i];
                  // do not allow porosity larger than 1.0
                  if (old_porosity[i] + porosity_change > 1.0)
                    porosity_change = 1.0 - old_porosity[i];
                  // now we calculate the equilibrium concentration of each component
                  std::vector<double> melting_points(n_components);
                  std::vector<double> eq_consts(n_components);
                  std::vector<double> c_solid_eq_values(n_components);
                  std::vector<double> c_liquid_eq_values(n_components);
                  if (enable_equilibrium_calculation)
                    {
                      for (unsigned int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                        {
                          melting_points[comp_idx] = temperature_melting(std::max(0.0, in.pressure[i]),
                                                                    melting_point_0_values[comp_idx],
                                                                    melting_curve_coefficient_A_values[comp_idx],
                                                                    melting_curve_coefficient_B_values[comp_idx]);
                          eq_consts[comp_idx] = equilibrium_constant(temperature_for_equilibrium_calculation,
                                                              latent_heat_values[comp_idx],
                                                              tuning_parameter_values[comp_idx],
                                                              melting_points[comp_idx]);
                          c_solid_eq_values[comp_idx] = calculate_concentration_solid(bulk_concentrations[comp_idx],
                                                                                eq_melt_fraction,
                                                                                eq_consts[comp_idx]);
                          c_liquid_eq_values[comp_idx] = calculate_concentration_liquid(bulk_concentrations[comp_idx],
                                                                                  eq_melt_fraction,
                                                                                  eq_consts[comp_idx]);
                        }
                    }
                  // now we can calculate the reaction rate of each component
                  for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                    {
                      // temporary implementation without reaction rate output

                      if (c == porosity_idx && this->get_timestep_number() > 0)
                        // this term will be treated as melting rate, it should multiply density and divide timestep
                        // we decide to calculate reaction term of porosity field
                        // whether equilibrium calculation is enabled or disabled
                        // for we've already deal with the disabled case above for porosity
                        out.reaction_terms[i][c] = porosity_change * out.densities[i] / this->get_timestep();
                      else if (this->get_timestep_number() > 0 && enable_equilibrium_calculation)
                      {
                        // if equilibrium calculation is disabled
                        // the method "match_component_index_to_composition_indices" would fail
                        // and solid_indices and liquid_indices would be meaningless
                        double c_eq = 0.0;
                        for (unsigned int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                          {
                            if (c == solid_indices[comp_idx])
                              {
                                c_eq = c_solid_eq_values[comp_idx];
                                break;
                              }
                            else if (c == liquid_indices[comp_idx])
                              {
                                c_eq = c_liquid_eq_values[comp_idx];
                                break;
                              }
                          }
                        if (prescribed_field_out != nullptr)
                          {
                            prescribed_field_out->prescribed_field_outputs[i][c] = c_eq;
                          }
                        out.reaction_terms[i][c] = (c_eq - old_fields[i][c]);
                      }
                      else
                      {
                        out.reaction_terms[i][c] = 0.0;
                        if (prescribed_field_out != nullptr)
                          {
                            prescribed_field_out->prescribed_field_outputs[i][c] = old_fields[i][c];
                          }
                      }

                      // fill reaction rate outputs if the model uses operator splitting
                      if (this->get_parameters().use_operator_splitting)
                        {
                          AssertThrow(this->get_melt_handler().melt_parameters.melt_without_porosity_advection_field == false,
                                      ExcMessage("Operator splitting with prescribed porosity field is not supported."));
                          if (reaction_rate_out != nullptr)
                            {
                              if (c == porosity_idx && this->get_timestep_number() > 0)
                                reaction_rate_out->reaction_rates[i][c] = porosity_change / melting_time_scale;
                              else if (this->get_timestep_number() > 0 && enable_equilibrium_calculation)
                              {
                                double c_eq = 0.0;
                                for (unsigned int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                                  {
                                    if (c == solid_indices[comp_idx])
                                      {
                                        c_eq = c_solid_eq_values[comp_idx];
                                        break;
                                      }
                                    else if (c == liquid_indices[comp_idx])
                                      {
                                        c_eq = c_liquid_eq_values[comp_idx];
                                        break;
                                      }
                                  }
                                reaction_rate_out->reaction_rates[i][c] = (c_eq - in.composition[i][c]) / melting_time_scale;
                              }
                              else
                                reaction_rate_out->reaction_rates[i][c] = 0.0;
                            }
                          out.reaction_terms[i][c] = 0.0;
                        }
                    }

                  if (prescribed_field_out != nullptr)
                    {
                      // temporary output c_bulk_an for debugging
                      // prescribed_field_out->prescribed_field_outputs[i][5] = bulk_concentrations[0]; // c_bulk_ol
                      // prescribed_field_out->prescribed_field_outputs[i][6] = bulk_concentrations[1]; // c_bulk_an
                      // prescribed_field_out->prescribed_field_outputs[i][5] = old_fields[i][liquid_indices[1]];
                      // prescribed_field_out->prescribed_field_outputs[i][6] = in.composition[i][liquid_indices[1]];
                    }

                }  
            }

          out.entropy_derivative_pressure[i]    = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = 0.0;

          double visc_temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
          {
            const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
            visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
          }
          else if (thermal_viscosity_exponent != 0.0)
          {
            const double delta_temp = in.temperature[i]-reference_T;
            visc_temperature_dependence = std::max(std::min(std::exp(-thermal_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
          }
          out.viscosities[i] *= visc_temperature_dependence;
        }    

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();

      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              // this is old porosity
              double porosity = std::max(in.composition[i][porosity_idx],0.0);

              // seems unnecessary to fill melt_out for porosity here
              // for we've fill it as the prescribed field.
              melt_out->porosities[i] = 0.0;

              melt_out->fluid_viscosities[i] = eta_f;
              melt_out->permeabilities[i] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<2>(1.0-porosity);
              melt_out->fluid_density_gradients[i] = Tensor<1,dim>();

              // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
              double temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                          * thermal_expansivity;
              else
                temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
              melt_out->fluid_densities[i] = reference_rho_f * temperature_dependence
                                             * std::exp(melt_compressibility * (in.pressure[i] - this->get_surface_pressure()));

              melt_out->compaction_viscosities[i] = xi_0 * std::exp(- alpha_phi * porosity);

              double visc_temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                {
                  const double delta_temp = in.temperature[i]-this->get_adiabatic_conditions().temperature(in.position[i]);
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/this->get_adiabatic_conditions().temperature(in.position[i])),1e4),1e-4);
                }
              else if (thermal_viscosity_exponent != 0.0)
                {
                  const double delta_temp = in.temperature[i]-reference_T;
                  visc_temperature_dependence = std::max(std::min(std::exp(-thermal_bulk_viscosity_exponent*delta_temp/reference_T),1e4),1e-4);
                }
              melt_out->compaction_viscosities[i] *= visc_temperature_dependence;
            }
        }
    }

    // the new declare_parameters function
    template <int dim>
    void
    MeltThermodynamicEquilibrium<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt thermodynamic equilibrium");
        {
          prm.declare_entry ("Reference solid density", "3000.",
                             Patterns::Double (0.),
                             "Reference density of the solid $\\rho_{s,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference melt density", "2500.",
                             Patterns::Double (0.),
                             "Reference density of the melt/fluid$\\rho_{f,0}$. "
                             "Units: \\si{\\kilogram\\per\\meter\\cubed}.");
          prm.declare_entry ("Reference temperature", "293.",
                             Patterns::Double (0.),
                             "The reference temperature $T_0$. The reference temperature is used "
                             "in both the density and viscosity formulas. Units: \\si{\\kelvin}.");
          prm.declare_entry ("Reference shear viscosity", "5e20",
                             Patterns::Double (0.),
                             "The value of the constant viscosity $\\eta_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference bulk viscosity", "1e22",
                             Patterns::Double (0.),
                             "The value of the constant bulk viscosity $\\xi_0$ of the solid matrix. "
                             "This viscosity may be modified by both temperature and porosity "
                             "dependencies. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference melt viscosity", "10.",
                             Patterns::Double (0.),
                             "The value of the constant melt viscosity $\\eta_f$. Units: \\si{\\pascal\\second}.");
          prm.declare_entry ("Reference specific heat", "1250.",
                              Patterns::Double (0.),
                              "The value of the specific heat $C_p$. "
                              "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}.");
          prm.declare_entry ("Reference permeability", "1e-8",
                                Patterns::Double(),
                                "Reference permeability of the solid host rock."
                                "Units: \\si{\\meter\\squared}.");

          prm.declare_entry ("Thermal viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of the shear viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Thermal bulk viscosity exponent", "0.0",
                             Patterns::Double (0.),
                             "The temperature dependence of the bulk viscosity. Dimensionless exponent. "
                             "See the general documentation "
                             "of this model for a formula that states the dependence of the "
                             "viscosity on this factor, which is called $\\beta$ there.");
          prm.declare_entry ("Exponential melt weakening factor", "27.",
                              Patterns::Double (0.),
                              "The porosity dependence of the viscosity. Units: dimensionless.");
          prm.declare_entry ("Thermal expansion coefficient", "2e-5",
                                Patterns::Double (0.),
                                "The value of the thermal expansion coefficient $\\beta$. "
                                "Units: \\si{\\per\\kelvin}.");
          prm.declare_entry ("Solid compressibility", "0.0",
                                  Patterns::Double (0.),
                                  "The value of the compressibility of the solid matrix. "
                                  "Units: \\si{\\per\\pascal}.");
          prm.declare_entry ("Melt compressibility", "0.0",
                                  Patterns::Double (0.),
                                  "The value of the compressibility of the melt. "
                                  "Units: \\si{\\per\\pascal}.");
          
          prm.declare_entry ("Enable equilibrium calculation", "true",
                             Patterns::Bool (),
                             "Whether to enable the calculation of the equilibrium melt fraction "
                             "based on the local temperature, pressure, and composition. If false, "
                             "then the model degrades to a simple melt transport model "
                             "without thermodynamic equilibrium calculation and melting/freezing source term.");

          prm.declare_entry ("Equilibrium solving method", "bisection",
                             Patterns::Selection ("bisection|newton"),
                             "The method used to solve the equation for the melt fraction. "
                             "Either bisection or newton.");

          prm.declare_entry ("Thermal conductivity", "4.7",
                             Patterns::Double (0.),
                             "The value of the thermal conductivity $k$. "
                             "Units: \\si{\\watt\\per\\meter\\per\\kelvin}.");
          prm.declare_entry ("Include melting and freezing", "true",
                             Patterns::Bool (),
                             "Whether to include melting and freezing (according to a simplified "
                             "linear melting approximation in the model (if true), or not (if "
                             "false).");
          prm.declare_entry ("Melting time scale for operator splitting", "1e3",
                             Patterns::Double (0.),
                             "In case the operator splitting scheme is used, the porosity field can not "
                             "be set to a new equilibrium melt fraction instantly, but the model has to "
                             "provide a melting time scale instead. This time scale defines how fast melting "
                             "happens, or more specifically, the parameter defines the time after which "
                             "the deviation of the porosity from the equilibrium melt fraction will be "
                             "reduced to a fraction of $1/e$. So if the melting time scale is small compared "
                             "to the time step size, the reaction will be so fast that the porosity is very "
                             "close to the equilibrium melt fraction after reactions are computed. Conversely, "
                             "if the melting time scale is large compared to the time step size, almost no "
                             "melting and freezing will occur."
                             "\n\n"
                             "Also note that the melting time scale has to be larger than or equal to the reaction "
                             "time step used in the operator splitting scheme, otherwise reactions can not be "
                             "computed. If the model does not use operator splitting, this parameter is not used. "
                             "Units: yr or s, depending on the ``Use years "
                             "in output instead of seconds'' parameter.");

          // Chemical component name list
          // The user need to enter chemical component names manually in the material model section of the parameter file
          // while they are actually defined in the compositional fields section already.
          // One chemical component corresponds to two compositional fields: one for solid and one for liquid
          prm.declare_entry ("Number of chemical components", "0",
                             Patterns::Integer (0),
                             "Number of chemical components used in the melt model. "
                             "Each chemical component corresponds to two compositional fields: "
                             "one for the solid part and one for the liquid part. "
                             "For example, if there are two chemical components, "
                             "then there should be four compositional fields defined in the "
                             "'Compositional fields' section of the parameter file: "
                             "'component1_solid', 'component1_liquid', 'component2_solid', 'component2_liquid'.");

          prm.declare_entry ("Chemical component names", "",
                             Patterns::List (Patterns::Anything()),
                             "List of names of chemical components used in the melt model. "
                             "Each name corresponds to two compositional fields: "
                             "one for the solid part and one for the liquid part. "
                             "The order of the names has to correspond to the order of the "
                             "compositional fields defined in the 'Compositional fields' section of the parameter file. "
                             "For example, if the compositional fields are defined as "
                             "'peridotite_solid', 'peridotite_liquid', 'basalt_solid', 'basalt_liquid', "
                             "then the chemical component names should be defined as 'peridotite, basalt'.");
          
          // Phase suffixes for solid and liquid compositional fields
          prm.declare_entry ("Phase suffix for solid compositional fields", "_solid",
                              Patterns::Anything(),
                              "Suffix used to identify solid compositional fields."
                              "For example, if the solid compositional fields are named '_solid', "
                              "then the suffix should be defined as '_solid'."
                            );
          prm.declare_entry ("Phase suffix for liquid compositional fields", "_liquid",
                              Patterns::Anything(),
                              "Suffix used to identify liquid compositional fields."
                              "For example, if the liquid compositional fields are named '_liquid', "
                              "then the suffix should be defined as '_liquid'."
                            );

          prm.declare_entry ("Melting point for each component at surface", // name
                             "", // default value
                             Patterns::List(Patterns::Double(0)), // pattern
                             "As we follow Keller and Katz (2016), the melting point here is in Celsius. "
                             "List of melting points at surface for each component. "
                             "T_m = T_m0 + A * P + B * P^2. "
                             "Units: \\si{\\kelvin}" // description
                            );
          prm.declare_entry ("Melting curve coefficient A for each component", // name
                             "", // default value
                             Patterns::List(Patterns::Double()), // pattern
                             "List of melting curve coefficients A for each component. "
                             "T_m = T_m0 + A * P + B * P^2. "
                             "Units: \\si{\\kelvin\\per\\pascal}" // description
                            );
          prm.declare_entry ("Melting curve coefficient B for each component", // name
                             "", // default value
                             Patterns::List(Patterns::Double()), // pattern
                             "List of melting curve coefficients B for each component. "
                             "T_m = T_m0 + A * P + B * P^2. "
                             "Units: \\si{\\kelvin\\per\\pascal\\square}" // description
                            ); 
          prm.declare_entry ("Latent heat for each component", // name
                             "", // default value
                             Patterns::List(Patterns::Double(0)), // pattern
                             "List of latent heats for each component. "
                             "If a unit of component melts, it consumes latent heat L. "
                             "Units: \\si{\\joule\\per\\kilogram}" // description
                            );
          prm.declare_entry ("Tuning parameter for each component", // name
                             "", // default value
                             Patterns::List(Patterns::Double(0)), // pattern
                             "List of tuning parameters for each component. "
                             "The tuning parameter controls the equilibrium constant K. "
                             "Units: \\si{\\joule\\per\\kelvin\\per\\kilogram}" // description
                            );   
                            
          prm.declare_entry ("Background porosity", 
                             "0.0",
                             Patterns::Double (0.0),
                             "A small background porosity added to avoid zero permeability.");    

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltThermodynamicEquilibrium<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Material model");
      {
        prm.enter_subsection("Melt thermodynamic equilibrium");
        {
          reference_rho_s                   = prm.get_double ("Reference solid density");
          reference_rho_f                   = prm.get_double ("Reference melt density");
          reference_T                       = prm.get_double ("Reference temperature");
          eta_0                             = prm.get_double ("Reference shear viscosity");
          xi_0                              = prm.get_double ("Reference bulk viscosity");
          eta_f                             = prm.get_double ("Reference melt viscosity");
          reference_permeability            = prm.get_double ("Reference permeability");
          reference_specific_heat           = prm.get_double ("Reference specific heat");

          thermal_viscosity_exponent        = prm.get_double ("Thermal viscosity exponent");
          thermal_bulk_viscosity_exponent   = prm.get_double ("Thermal bulk viscosity exponent");
          thermal_expansivity               = prm.get_double ("Thermal expansion coefficient");
          alpha_phi                         = prm.get_double ("Exponential melt weakening factor");
          compressibility                   = prm.get_double ("Solid compressibility");
          melt_compressibility              = prm.get_double ("Melt compressibility");

          thermal_conductivity              = prm.get_double ("Thermal conductivity");
          include_melting_and_freezing      = prm.get_bool ("Include melting and freezing");
          melting_time_scale                = prm.get_double ("Melting time scale for operator splitting");

          enable_equilibrium_calculation    = prm.get_bool ("Enable equilibrium calculation");

          equilibrium_solving_method = prm.get ("Equilibrium solving method");

          AssertThrow(!(this->get_melt_handler().melt_parameters.melt_without_porosity_advection_field 
                        && this->get_parameters().use_operator_splitting),
                      ExcMessage("Error: Material model Melt thermodynamic equilibrium "
                                 "can not use the option 'Melt without porosity advection field' together with operator splitting."));

          // parse chemical component number
          n_components = prm.get_integer ("Number of chemical components");

          // parse chemical component names
          component_names = Utilities::split_string_list(prm.get ("Chemical component names"));
          AssertThrow(component_names.size() == n_components,
                      ExcMessage("Error: The number of chemical component names provided (" +
                                 Utilities::to_string(component_names.size()) + ") does not match the number of chemical components defined (" +
                                 Utilities::to_string(n_components) + ")."));

          // parse phase suffixes
          component_name_solid_suffix = prm.get ("Phase suffix for solid compositional fields");
          component_name_liquid_suffix = prm.get ("Phase suffix for liquid compositional fields");

          melting_point_0_values = 
            Utilities::possibly_extend_from_1_to_N (
              Utilities::string_to_double(
                Utilities::split_string_list(prm.get ("Melting point for each component at surface"))
              ),
              n_components, // we only need to write a number once for each component
              "Melting point for each component at surface"
            );
          melting_curve_coefficient_A_values = 
            Utilities::possibly_extend_from_1_to_N (
              Utilities::string_to_double(
                Utilities::split_string_list(prm.get ("Melting curve coefficient A for each component"))
              ),
              n_components, // we only need to write a number once for each component
              "Melting curve coefficient A for each component"
            );
          melting_curve_coefficient_B_values = 
            Utilities::possibly_extend_from_1_to_N (
              Utilities::string_to_double(
                Utilities::split_string_list(prm.get ("Melting curve coefficient B for each component"))
              ),
              n_components, // we only need to write a number once for each component
              "Melting curve coefficient B for each component"
            );
          latent_heat_values = 
            Utilities::possibly_extend_from_1_to_N (
              Utilities::string_to_double(
                Utilities::split_string_list(prm.get ("Latent heat for each component"))
              ),
              n_components, // we only need to write a number once for each component
              "Latent heat for each component"
            );
          tuning_parameter_values = 
            Utilities::possibly_extend_from_1_to_N (
              Utilities::string_to_double(
                Utilities::split_string_list(prm.get ("Tuning parameter for each component"))
              ),
              n_components, // we only need to write a number once for each component
              "Tuning parameter for each component"
            );
          
          background_porosity = prm.get_double("Background porosity");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt thermodynamic equilibrium with Thermal viscosity exponent can not have reference_T=0."));

          if (this->convert_output_to_years() == true)
            melting_time_scale *= year_in_seconds;

          if (this->get_parameters().use_operator_splitting)
            {
              AssertThrow(melting_time_scale >= this->get_parameters().reaction_time_step,
                          ExcMessage("The reaction time step " + Utilities::to_string(this->get_parameters().reaction_time_step)
                                     + " in the operator splitting scheme is too large to compute melting rates! "
                                     "You have to choose it in such a way that it is smaller than the 'Melting time scale for "
                                     "operator splitting' chosen in the material model, which is currently "
                                     + Utilities::to_string(melting_time_scale) + "."));
              AssertThrow(melting_time_scale > 0,
                          ExcMessage("The Melting time scale for operator splitting must be larger than 0!"));
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Melt thermodynamic equilibrium with melt transport only "
                                     "works if there is a compositional field called porosity."));
            }

          if (this->include_melt_transport())
            {
              AssertThrow(this->introspection().compositional_name_exists("porosity"),
                          ExcMessage("Material model Melt thermodynamic equilibrium with melt transport only "
                                     "works if there is a compositional field called porosity."));
              if (include_melting_and_freezing)
                {
                  // check if components matches compositional fields
                  // if equilibrium calculation is disabled, this check is always true
                  const bool match_result = match_component_index_to_composition_indices();
                  AssertThrow(match_result,
                              ExcMessage("Error: Material model Melt thermodynamic equilibrium could not match the chemical components to compositional fields. "
                                         "Please check if the compositional fields for each chemical component are defined correctly according to the naming convention."));

                  if (enable_equilibrium_calculation)
                  {
                    // check if the number of melting parameters matches n_components
                    AssertThrow(melting_point_0_values.size() == n_components,
                                ExcMessage("Error: Material model Melt thermodynamic equilibrium needs the same number of melting points at surface as the number of compositional fields participating in melting. "
                                            "There are " + Utilities::to_string(n_components) + " compositional fields participating in melting, but you provided "
                                            + Utilities::to_string(melting_point_0_values.size()) + " melting points at surface."));
                    AssertThrow(melting_curve_coefficient_A_values.size() == n_components,
                                ExcMessage("Error: Material model Melt thermodynamic equilibrium needs the same number of melting curve coefficients A as the number of compositional fields participating in melting. "
                                            "There are " + Utilities::to_string(n_components) + " compositional fields participating in melting, but you provided "
                                            + Utilities::to_string(melting_curve_coefficient_A_values.size()) + " melting curve coefficients A."));
                    AssertThrow(melting_curve_coefficient_B_values.size() == n_components,
                                ExcMessage("Error: Material model Melt thermodynamic equilibrium needs the same number of melting curve coefficients B as the number of compositional fields participating in melting. "
                                            "There are " + Utilities::to_string(n_components) + " compositional fields participating in melting, but you provided "
                                            + Utilities::to_string(melting_curve_coefficient_B_values.size()) + " melting curve coefficients B."));
                    AssertThrow(latent_heat_values.size() == n_components,
                                ExcMessage("Error: Material model Melt thermodynamic equilibrium needs the same number of latent heats as the number of compositional fields participating in melting. "
                                            "There are " + Utilities::to_string(n_components) + " compositional fields participating in melting, but you provided "
                                            + Utilities::to_string(latent_heat_values.size()) + " latent heats."));
                    AssertThrow(tuning_parameter_values.size() == n_components,
                                ExcMessage("Error: Material model Melt thermodynamic equilibrium needs the same number of tuning parameters as the number of compositional fields participating in melting. "
                                            "There are " + Utilities::to_string(n_components) + " compositional fields participating in melting, but you provided "
                                            + Utilities::to_string(tuning_parameter_values.size()) + " tuning parameters."));
                  }
                }
            }
          
          if (equilibrium_solving_method != "bisection")
            {
              AssertThrow(false, ExcMessage(
                "Error: Material model Melt thermodynamic equilibrium currently only supports bisection method to solve equations."
                "please set the parameter 'Equilibrium solving method' to 'bisection'."));
            }
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    MeltThermodynamicEquilibrium<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (this->get_melt_handler().melt_parameters.melt_without_porosity_advection_field
          && out.template get_additional_output<MaterialModel::PrescribedFieldOutputs<dim>>() == nullptr)
        {
          // AssertThrow(false, ExcMessage("Successfully run in to create_additional_named_outputs method."))
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            // we create prescribed field outputs for all compositional fields
            // but we only use the first one for porosity. [q][0] is the porosity at point q.
            std::make_unique<MaterialModel::PrescribedFieldOutputs<dim>> (n_points, this->n_compositional_fields())
          );
        }
      if (this->get_parameters().use_operator_splitting
          && out.template get_additional_output<ReactionRateOutputs<dim>>() == nullptr)
        {
          const unsigned int n_points = out.n_evaluation_points();
          out.additional_outputs.push_back(
            std::make_unique<MaterialModel::ReactionRateOutputs<dim>> (n_points, this->n_compositional_fields()));
        }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    ASPECT_REGISTER_MATERIAL_MODEL(MeltThermodynamicEquilibrium,
                                   "melt thermodynamic equilibrium",
                                   "A material model that implements ...")
  }
}
