/*
  Copyright (C) 2011 - 2023 by the authors of the ASPECT code.

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

#include <aspect/initial_composition/thermodynamic_equilibrium.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/geometry_model/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/material_model/interface.h>
#include <aspect/melt.h>
#include <aspect/simulator_signals.h>
#include <aspect/utilities.h>

#include <aspect/structured_data.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/signaling_nan.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iterator>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>

namespace aspect
{
  namespace InitialComposition
  {
    namespace
    {
      // The pressure threshold / melting curve machinery of the material model
      // is not needed here; everything thermodynamic goes through the material
      // model. These helpers only deal with the solitary wave profile.
      // pretty() prints doubles with six decimals, which is useless for
      // the residuals we report in error messages.
      std::string
      pretty (const double value)
      {
        std::ostringstream stream;
        stream.precision(12);
        stream << value;
        return stream.str();
      }



      double
      implicit_function (const double relative_porosity,
                         const double amplitude)
      {
        if (relative_porosity < 1.0 || relative_porosity > amplitude)
          return std::numeric_limits<double>::infinity();

        if (relative_porosity == 1.0)
          return std::numeric_limits<double>::infinity();

        const double a1 = std::sqrt(amplitude - 1.0);
        const double a2 = std::sqrt(amplitude - relative_porosity);
        const double a3 = std::sqrt(amplitude + 0.5);

        return std::abs(a3 * (2.0 * a2
                              - (1.0 / a1) * std::log((a1 - a2) / (a1 + a2))));
      }
    }



    template <int dim>
    ThermodynamicEquilibrium<dim>::ThermodynamicEquilibrium ()
      :
      specification (Specification::melt_fraction),
      mixing_family (MixingFamily::tie_line),
      background_model (FieldModel::uniform),
      target_model (FieldModel::background),
      solid_model (FieldModel::uniform),
      liquid_model (FieldModel::uniform),
      uniform_target (0.0),
      family_component (0),
      consistency_tolerance (1e-8),
      inconsistent_behaviour (InconsistentBehaviour::error),
      use_solitary_wave (false),
      peak_position (0.0),
      peak_porosity (0.0),
      truncation_threshold (1e-4),
      derive_compaction_length (true),
      prescribed_compaction_length (0.0),
      solving_tolerance (1e-14),
      report_initial_equilibrium (true),
      dump_generated_state (false),
      porosity_index (numbers::invalid_unsigned_int),
      notified_about_reference_profile (false),
      derived_parameters_are_ready (false),
      solitary_wave_is_possible (false),
      peak_background_melt_fraction (0.0),
      relative_amplitude (0.0),
      compaction_length (0.0),
      window_radius (0.0),
      max_equilibrium_residual (0.0),
      max_phase_residual (0.0),
      max_composition_sum_error (0.0),
      n_evaluated_points (0),
      n_solved_points (0),
      n_kept_background_points (0)
    {}



    template <int dim>
    void
    ThermodynamicEquilibrium<dim>::initialize ()
    {
      // Which chemical components does the material model have?
      components = this->get_material_model().get_equilibrium_components();

      AssertThrow(!components.empty(),
                  ExcMessage("The initial composition model 'thermodynamic equilibrium' requires a "
                             "material model that computes a solid-liquid thermodynamic equilibrium, "
                             "but the current material model does not provide any equilibrium "
                             "components (get_equilibrium_components() returned an empty list)."));

      const unsigned int n_components = components.size();
      bulk_indices.resize(n_components);
      solid_indices.resize(n_components);
      liquid_indices.resize(n_components);

      for (unsigned int c = 0; c < n_components; ++c)
        {
          bulk_indices[c] = this->introspection().compositional_index_for_name(components[c].name);

          try
            {
              solid_indices[c] = this->introspection().compositional_index_for_name(components[c].name
                                                                                    + components[c].solid_field_suffix);
            }
          catch (...)
            {
              solid_indices[c] = numbers::invalid_unsigned_int;
            }

          liquid_indices[c] = this->introspection().compositional_index_for_name(components[c].name
                                                                                 + components[c].liquid_field_suffix);
        }

      porosity_index = this->introspection().compositional_index_for_name("porosity");

      if (specification == Specification::phase_compositions)
        {
          if (solid_model == FieldModel::uniform)
            AssertThrow(uniform_solid.size() == n_components,
                        ExcMessage("'Uniform solid composition' has "
                                   + pretty(uniform_solid.size())
                                   + " entries, but the material model has "
                                   + pretty(n_components) + " chemical components."));

          if (liquid_model == FieldModel::uniform)
            AssertThrow(uniform_liquid.size() == n_components,
                        ExcMessage("'Uniform liquid composition' has "
                                   + pretty(uniform_liquid.size())
                                   + " entries, but the material model has "
                                   + pretty(n_components) + " chemical components."));
        }

      // Keep the initial temperature model alive: the simulator releases its
      // own pointer after the first time step, but this plugin keeps being
      // called through the boundary composition model afterwards.
      initial_temperature_manager = this->get_initial_temperature_manager_pointer();

      // The objects themselves were created in parse_parameters(); here we
      // only finish their initialization.
      if (background_profile != nullptr)
        background_profile->initialize(this->get_mpi_communicator());

      if (target_profile != nullptr)
        target_profile->initialize(this->get_mpi_communicator());

      if (solid_profile != nullptr)
        solid_profile->initialize(this->get_mpi_communicator());

      if (liquid_profile != nullptr)
        liquid_profile->initialize(this->get_mpi_communicator());

      // Report once the initial state has been set (the initial composition
      // manager is still alive at that point).
      if (report_initial_equilibrium)
        this->get_signals().post_set_initial_state.connect(
          [this] (const SimulatorAccess<dim> &)
        {
          this->report();
        });
    }



    template <int dim>
    std::vector<double>
    ThermodynamicEquilibrium<dim>::
    background_bulk_composition (const Point<dim> &position) const
    {
      const unsigned int n_components = components.size();
      std::vector<double> composition(n_components, 0.0);

      switch (background_model)
        {
          case FieldModel::uniform:
            composition = uniform_background;
            break;

          case FieldModel::depth_profile:
            {
              const double depth = this->get_geometry_model().depth(position);
              for (unsigned int c = 0; c < n_components; ++c)
                composition[c] = background_profile->get_data_component(Point<1>(depth),
                                                                        background_profile->get_column_index_from_name(components[c].name));
              break;
            }

          default:
            AssertThrow(false, ExcNotImplemented());
        }

      // Normalize, so that the family members always sum to one.
      const double sum = std::accumulate(composition.begin(), composition.end(), 0.0);
      AssertThrow(sum > 0.0,
                  ExcMessage("The background bulk composition sums to zero at position ("
                             + pretty(position[0]) + ", "
                             + pretty(position[1]) + ")."));
      for (double &value : composition)
        value /= sum;

      return composition;
    }



    template <int dim>
    double
    ThermodynamicEquilibrium<dim>::
    base_target_melt_fraction (const Point<dim> &position,
                               const double      background_melt_fraction) const
    {
      switch (target_model)
        {
          case FieldModel::background:
            return background_melt_fraction;

          case FieldModel::uniform:
            return uniform_target;

          case FieldModel::depth_profile:
            {
              const double depth = this->get_geometry_model().depth(position);
              return target_profile->get_data_component(Point<1>(depth), 0);
            }

          case FieldModel::function:
            return target_function->value(position, 0);

          default:
            AssertThrow(false, ExcNotImplemented());
            return 0.0;
        }
    }



    template <int dim>
    void
    ThermodynamicEquilibrium<dim>::compute_derived_parameters () const
    {
      derived_parameters_are_ready = true;
      solitary_wave_is_possible = false;
      peak_background_melt_fraction = 0.0;
      relative_amplitude = 0.0;
      compaction_length = 0.0;
      window_radius = 0.0;

      if (!use_solitary_wave)
        return;

      // Before the reference profile exists we cannot evaluate anything.
      if (!this->get_adiabatic_conditions().is_initialized())
        {
          derived_parameters_are_ready = false;
          return;
        }

      // The peak state: use the centre of the domain in the horizontal
      // directions and the prescribed vertical position.
      Point<dim> peak_point = this->get_geometry_model().representative_point(0.0);
      peak_point[dim-1] = peak_position;

      const double pressure = this->get_adiabatic_conditions().pressure(peak_point);
      const double temperature = initial_temperature_manager->initial_temperature(peak_point);
      const std::vector<double> background = background_bulk_composition(peak_point);

      double f0 = 0.0;
      std::vector<double> c_solid, c_liquid;
      evaluate_equilibrium(pressure, temperature, background, f0, c_solid, c_liquid);

      AssertThrow(f0 > 0.0,
                  ExcMessage("The solitary wave cannot be generated: the background composition is "
                             "below the solidus at the requested peak position ("
                             + pretty(peak_position) + "). Move the peak into the "
                             "melting region or enrich the background composition."));

      peak_background_melt_fraction = f0;
      relative_amplitude = peak_porosity / f0;

      AssertThrow(relative_amplitude > 1.0,
                  ExcMessage("The solitary wave cannot be generated: the prescribed peak porosity ("
                             + pretty(peak_porosity) + ") is not larger than the "
                             "background melt fraction at the peak ("
                             + pretty(f0) + ")."));

      // The compaction length, from the material model itself.
      if (derive_compaction_length)
        {
          const unsigned int n_fields = this->n_compositional_fields();

          MaterialModel::MaterialModelInputs<dim> in(1, n_fields);
          MaterialModel::MaterialModelOutputs<dim> out(1, n_fields);
          this->get_material_model().create_additional_named_outputs(out);

          if (out.template get_additional_output<MaterialModel::MeltOutputs<dim>>() == nullptr)
            out.additional_outputs.push_back(
              std::make_unique<MaterialModel::MeltOutputs<dim>>(1, n_fields));

          in.position[0] = peak_point;
          in.temperature[0] = temperature;
          in.pressure[0] = pressure;
          in.pressure_gradient[0] = Tensor<1,dim>();
          in.velocity[0] = Tensor<1,dim>();
          in.strain_rate[0] = SymmetricTensor<2,dim>();
          in.composition[0].assign(n_fields, 0.0);

          for (unsigned int c = 0; c < components.size(); ++c)
            in.composition[0][bulk_indices[c]] = background[c];
          in.composition[0][porosity_index] = f0;

          in.requested_properties = in.requested_properties
                                    | MaterialModel::MaterialProperties::viscosity;

          this->get_material_model().evaluate(in, out);

          MaterialModel::MeltOutputs<dim> *melt_out =
            out.template get_additional_output<MaterialModel::MeltOutputs<dim>>();

          AssertThrow(melt_out != nullptr,
                      ExcMessage("The material model did not create MeltOutputs, so the "
                                 "compaction length cannot be derived. Set "
                                 "'Derive compaction length = false' and prescribe it instead."));

          const double permeability = melt_out->permeabilities[0];
          const double compaction_viscosity = melt_out->compaction_viscosities[0];
          const double fluid_viscosity = melt_out->fluid_viscosities[0];
          const double shear_viscosity = out.viscosities[0];

          compaction_length = std::sqrt(permeability
                                        * (compaction_viscosity + 4.0 / 3.0 * shear_viscosity)
                                        / fluid_viscosity);
        }
      else
        compaction_length = prescribed_compaction_length;

      AssertThrow(compaction_length > 0.0 && std::isfinite(compaction_length),
                  ExcMessage("The compaction length used for the solitary wave is not positive "
                             "and finite. Check the material model's permeability and viscosities, "
                             "or prescribe the compaction length explicitly."));

      // The window radius follows from the truncation threshold.
      window_radius = compaction_length * implicit_function(1.0 + truncation_threshold,
                                                            relative_amplitude);

      solitary_wave_is_possible = true;

      if (this->get_pcout().is_active())
        {
          this->get_pcout() << "[TEQ initial condition] solitary wave parameters:" << std::endl
                            << "    peak position (vertical coordinate) = " << peak_position << " m" << std::endl
                            << "    background melt fraction at peak   = " << f0 << std::endl
                            << "    relative amplitude A               = " << relative_amplitude << std::endl
                            << "    compaction length delta            = " << compaction_length << " m" << std::endl
                            << "    window radius (threshold "
                            << truncation_threshold << ")      = " << window_radius << " m" << std::endl;
        }
    }



    template <int dim>
    double
    ThermodynamicEquilibrium<dim>::
    solitary_wave_implicit_function (const double relative_porosity) const
    {
      return implicit_function(relative_porosity, relative_amplitude);
    }



    template <int dim>
    double
    ThermodynamicEquilibrium<dim>::
    solitary_wave_profile (const double scaled_distance) const
    {
      // Solve implicit_function(phi_r) = scaled_distance for phi_r, with
      // phi_r in [1, A]. The function is +infinity at phi_r = 1 and 0 at
      // phi_r = A, so the interval is always a valid bracket.
      double lower = 1.0;
      double upper = relative_amplitude;

      for (unsigned int iteration = 0; iteration < 200; ++iteration)
        {
          const double middle = 0.5 * (lower + upper);

          if (solitary_wave_implicit_function(middle) > scaled_distance)
            lower = middle;
          else
            upper = middle;

          if (upper - lower < 1e-15 * std::max(1.0, upper))
            break;
        }

      return 0.5 * (lower + upper);
    }



    template <int dim>
    double
    ThermodynamicEquilibrium<dim>::
    target_melt_fraction (const Point<dim> &position,
                          const double      background_melt_fraction) const
    {
      const double base = base_target_melt_fraction(position, background_melt_fraction);

      if (!use_solitary_wave)
        return base;

      if (!derived_parameters_are_ready)
        compute_derived_parameters();

      if (!solitary_wave_is_possible)
        return base;

      const double distance = std::abs(position[dim-1] - peak_position);

      if (distance > window_radius)
        return base;

      // The wave multiplies the background melt fraction by the dimensionless
      // profile, so that phi* = f_bg * phi_r and the target vanishes wherever
      // the background is below the solidus.
      return background_melt_fraction * solitary_wave_profile(distance / compaction_length);
    }



    template <int dim>
    void
    ThermodynamicEquilibrium<dim>::
    evaluate_equilibrium (const double               pressure,
                          const double               temperature,
                          const std::vector<double> &bulk_composition,
                          double                    &melt_fraction,
                          std::vector<double>       &solid_composition,
                          std::vector<double>       &liquid_composition) const
    {
      const bool implemented =
        this->get_material_model().evaluate_equilibrium_state(pressure,
                                                              temperature,
                                                              bulk_composition,
                                                              melt_fraction,
                                                              solid_composition,
                                                              liquid_composition,
                                                              solving_tolerance);

      AssertThrow(implemented,
                  ExcMessage("The material model does not implement "
                             "evaluate_equilibrium_state(), which is required by the initial "
                             "composition model 'thermodynamic equilibrium'."));

      AssertThrow(std::isfinite(melt_fraction)
                  && melt_fraction >= 0.0 && melt_fraction <= 1.0,
                  ExcMessage("The material model returned an invalid melt fraction ("
                             + pretty(melt_fraction) + ") for P = "
                             + pretty(pressure) + " Pa, T = "
                             + pretty(temperature) + " K."));
    }



    template <int dim>
    typename ThermodynamicEquilibrium<dim>::EquilibriumState
    ThermodynamicEquilibrium<dim>::compute_state (const Point<dim> &position) const
    {
      EquilibriumState state;

      const unsigned int n_components = components.size();
      const std::vector<double> background = background_bulk_composition(position);

      // During the construction of the adiabatic reference profile the material
      // model -- and therefore this plugin -- is evaluated before the reference
      // pressure exists. Return the background composition in that case; it is
      // the natural reference composition, and none of the derived quantities
      // can be computed yet.
      if (!this->get_adiabatic_conditions().is_initialized())
        {
          if (!notified_about_reference_profile)
            {
              notified_about_reference_profile = true;
              this->get_pcout() << "[TEQ initial condition] the adiabatic reference profile is not "
                                   "available yet (it is being built by evaluating the material "
                                   "model right now); returning the background bulk composition for "
                                   "that pass." << std::endl;
            }

          state.bulk_composition = background;
          state.solid_composition.assign(n_components, 0.0);
          state.liquid_composition.assign(n_components, 0.0);
          return state;
        }

      state.pressure = this->get_adiabatic_conditions().pressure(position);
      state.temperature = initial_temperature_manager->initial_temperature(position);

      double f_background = 0.0;
      std::vector<double> c_solid_background, c_liquid_background;
      evaluate_equilibrium(state.pressure, state.temperature, background,
                           f_background, c_solid_background, c_liquid_background);

      const bool prescribe_phase_compositions =
        (specification == Specification::phase_compositions);

      const double phi_target = (specification == Specification::bulk_composition)
                                ? 0.0
                                : target_melt_fraction(position, f_background);

      // Finalize: check that the composition is valid and evaluate the
      // equilibrium there.
      const auto finalize = [&] (const std::vector<double> &composition) -> void
      {
        const double total = std::accumulate(composition.begin(), composition.end(), 0.0);
        max_composition_sum_error = std::max(max_composition_sum_error, std::abs(total - 1.0));

        AssertThrow(std::abs(total - 1.0) < 1e-10,
                    ExcMessage("The generated bulk composition does not sum to one (sum = "
                               + pretty(total) + ") at position ("
                               + pretty(position[0]) + ", "
                               + pretty(position[1]) + "). For "
                               "'Specification = phase compositions' this means that the "
                               "prescribed solid and liquid compositions do not satisfy the "
                               "phase mass balance."));

        state.bulk_composition = composition;
        evaluate_equilibrium(state.pressure, state.temperature, state.bulk_composition,
                             state.melt_fraction, state.solid_composition, state.liquid_composition);
      };

      // The two anchors of the one-dimensional family, for the mixing families
      // that are defined through the background state.
      const auto family_anchors = [&] (std::vector<double> &anchor_a,
                                       std::vector<double> &anchor_b) -> void
      {
        anchor_a.resize(n_components);
        anchor_b.resize(n_components);

        switch (mixing_family)
          {
            case MixingFamily::tie_line:
              anchor_a = c_solid_background;
              anchor_b = c_liquid_background;
              break;

            case MixingFamily::background_to_liquid:
              anchor_a = background;
              anchor_b = c_liquid_background;
              break;

            case MixingFamily::background_to_solid:
              anchor_a = background;
              anchor_b = c_solid_background;
              break;

            case MixingFamily::background_to_component:
              anchor_a = background;
              anchor_b.assign(n_components, 0.0);
              anchor_b[family_component] = 1.0;
              break;

            default:
              AssertThrow(false, ExcNotImplemented());
          }
      };

      // Solve for the mixing fraction s along
      //   composition(s) = (1-s) * anchor_a + s * anchor_b
      // such that the equilibrium melt fraction equals target. Along the tie
      // line the mixing fraction equals the melt fraction exactly, so no
      // iteration is required there.
      const auto solve_along_family =
        [&] (const std::vector<double> &anchor_a,
             const std::vector<double> &anchor_b,
             const double               target,
             const bool                 allow_analytic,
             std::vector<double>       &composition) -> void
      {
        composition.resize(n_components);

        if (allow_analytic && mixing_family == MixingFamily::tie_line)
          {
            for (unsigned int c = 0; c < n_components; ++c)
              composition[c] = (1.0 - target) * anchor_a[c] + target * anchor_b[c];
            return;
          }

        const auto melt_fraction_of_mixing_fraction = [&] (const double s) -> double
        {
          std::vector<double> trial(n_components);
          for (unsigned int c = 0; c < n_components; ++c)
            trial[c] = (1.0 - s) * anchor_a[c] + s * anchor_b[c];

          double f = 0.0;
          std::vector<double> cs, cl;
          evaluate_equilibrium(state.pressure, state.temperature, trial, f, cs, cl);
          return f;
        };

        const double f_at_zero = melt_fraction_of_mixing_fraction(0.0);
        const double f_at_one = melt_fraction_of_mixing_fraction(1.0);

        AssertThrow(target >= f_at_zero && target <= f_at_one,
                    ExcMessage("The target melt fraction (" + pretty(target)
                               + ") cannot be reached along the chosen mixing family: the family "
                               "spans [" + pretty(f_at_zero) + ", "
                               + pretty(f_at_one) + "] there. This happened at position ("
                               + pretty(position[0]) + ", "
                               + pretty(position[1]) + ")."));

        double lower = 0.0;
        double upper = 1.0;
        for (unsigned int iteration = 0; iteration < 100; ++iteration)
          {
            const double middle = 0.5 * (lower + upper);

            if (melt_fraction_of_mixing_fraction(middle) < target)
              lower = middle;
            else
              upper = middle;

            if (upper - lower < 1e-16)
              break;
          }

        const double s = 0.5 * (lower + upper);
        for (unsigned int c = 0; c < n_components; ++c)
          composition[c] = (1.0 - s) * anchor_a[c] + s * anchor_b[c];
      };

      if (specification == Specification::bulk_composition || phi_target <= 0.0)
        {
          // Nothing to invert: keep the background composition.
          ++n_kept_background_points;
          finalize(background);
          ++n_evaluated_points;
          return state;
        }

      AssertThrow(phi_target <= 1.0,
                  ExcMessage("The prescribed target melt fraction ("
                             + pretty(phi_target)
                             + ") is not smaller than one."));

      if (prescribe_phase_compositions)
        {
          // The bulk composition follows from the phase mass balance.
          const std::vector<double> prescribed_solid = solid_phase_composition(position);
          const std::vector<double> prescribed_liquid = liquid_phase_composition(position);

          std::vector<double> composition(n_components);
          for (unsigned int c = 0; c < n_components; ++c)
            composition[c] = (1.0 - phi_target) * prescribed_solid[c]
                             + phi_target * prescribed_liquid[c];

          finalize(composition);
        }
      else
        {
          if (mixing_family == MixingFamily::tie_line && f_background <= 0.0)
            AssertThrow(false,
                        ExcMessage("The mixing family 'tie line' is not defined where the "
                                   "background is below the solidus: the hypothetical liquid "
                                   "composition does not sum to one there. Either keep the target "
                                   "melt fraction at zero in that region, or use "
                                   "'background to component' instead. This happened at position ("
                                   + pretty(position[0]) + ", "
                                   + pretty(position[1]) + ")."));

          std::vector<double> anchor_a, anchor_b;
          family_anchors(anchor_a, anchor_b);

          std::vector<double> composition;
          solve_along_family(anchor_a, anchor_b, phi_target, true, composition);
          finalize(composition);
        }

      ++n_solved_points;

      // For 'Specification = phase compositions' check whether the
      // prescribed phase compositions are reproduced by the material model. The
      // liquid composition is only compared where there actually is melt: below
      // the solidus it is not defined and does not sum to one.
      if (prescribe_phase_compositions)
        {
          const std::vector<double> prescribed_solid = solid_phase_composition(position);
          const std::vector<double> prescribed_liquid = liquid_phase_composition(position);

          double residual = std::abs(state.melt_fraction - phi_target);
          for (unsigned int c = 0; c < n_components; ++c)
            residual = std::max(residual,
                                std::abs(state.solid_composition[c] - prescribed_solid[c]));

          if (state.melt_fraction > 0.0)
            for (unsigned int c = 0; c < n_components; ++c)
              residual = std::max(residual,
                                  std::abs(state.liquid_composition[c] - prescribed_liquid[c]));

          max_phase_residual = std::max(max_phase_residual, residual);

          if (residual > consistency_tolerance)
            {
              AssertThrow(inconsistent_behaviour != InconsistentBehaviour::error,
                          ExcMessage("The prescribed phase compositions and melt fraction are not "
                                     "consistent with the material model's equilibrium: the largest "
                                     "difference is " + pretty(residual) + " at position ("
                                     + pretty(position[0]) + ", "
                                     + pretty(position[1]) + ") (P = "
                                     + pretty(state.pressure) + " Pa, T = "
                                     + pretty(state.temperature) + " K). Set "
                                     "'If inconsistent' to 'keep melt fraction' or 'keep liquid' "
                                     "to fall back to a well-posed construction, or adjust the "
                                     "prescribed compositions."));

              // A single well-posed fallback.
              std::vector<double> anchor_a, anchor_b;
              if (inconsistent_behaviour == InconsistentBehaviour::keep_melt_fraction)
                {
                  family_anchors(anchor_a, anchor_b);
                  std::vector<double> composition;
                  solve_along_family(anchor_a, anchor_b, phi_target, true, composition);
                  finalize(composition);
                }
              else
                {
                  // Use the prescribed liquid composition as the enriched end member.
                  anchor_a = background;
                  anchor_b = prescribed_liquid;
                  std::vector<double> composition;
                  solve_along_family(anchor_a, anchor_b, phi_target, false, composition);
                  finalize(composition);
                }
            }
        }

      // Record the final residual, i.e. after any fallback has been applied.
      if (specification != Specification::bulk_composition)
        max_equilibrium_residual = std::max(max_equilibrium_residual,
                                            std::abs(state.melt_fraction - phi_target));

      ++n_evaluated_points;

      return state;
    }



    template <int dim>
    std::vector<double>
    ThermodynamicEquilibrium<dim>::
    solid_phase_composition (const Point<dim> &position) const
    {
      const unsigned int n_components = components.size();
      std::vector<double> composition(n_components, 0.0);

      if (solid_model == FieldModel::uniform)
        composition = uniform_solid;
      else if (solid_model == FieldModel::depth_profile)
        {
          const double depth = this->get_geometry_model().depth(position);
          for (unsigned int c = 0; c < n_components; ++c)
            composition[c] = solid_profile->get_data_component(Point<1>(depth),
                                                               solid_profile->get_column_index_from_name(components[c].name));
        }
      else
        AssertThrow(false, ExcNotImplemented());

      return composition;
    }



    template <int dim>
    std::vector<double>
    ThermodynamicEquilibrium<dim>::
    liquid_phase_composition (const Point<dim> &position) const
    {
      const unsigned int n_components = components.size();
      std::vector<double> composition(n_components, 0.0);

      if (liquid_model == FieldModel::uniform)
        composition = uniform_liquid;
      else if (liquid_model == FieldModel::depth_profile)
        {
          const double depth = this->get_geometry_model().depth(position);
          for (unsigned int c = 0; c < n_components; ++c)
            composition[c] = liquid_profile->get_data_component(Point<1>(depth),
                                                                liquid_profile->get_column_index_from_name(components[c].name));
        }
      else
        AssertThrow(false, ExcNotImplemented());

      return composition;
    }



    template <int dim>
    double
    ThermodynamicEquilibrium<dim>::
    value_for_field (const EquilibriumState &state,
                     const unsigned int      n_comp) const
    {
      for (unsigned int c = 0; c < components.size(); ++c)
        {
          if (n_comp == bulk_indices[c])
            return state.bulk_composition[c];

          if (n_comp == liquid_indices[c])
            return state.liquid_composition[c];

          if (solid_indices[c] != numbers::invalid_unsigned_int
              && n_comp == solid_indices[c])
            return state.solid_composition[c];
        }

      if (n_comp == porosity_index)
        return state.melt_fraction;

      // Any other field (melting rate, debug fields, ...) is filled either by
      // the material model or not at all; zero is a safe initial value.
      return 0.0;
    }



    template <int dim>
    double
    ThermodynamicEquilibrium<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int n_comp) const
    {
      // While the adiabatic reference profile is being built the equilibrium
      // cannot be evaluated yet (see compute_state()). Those provisional states
      // must not enter the cache: the points at which the profile is built are
      // not mesh nodes, and some of them coincide with mesh nodes, which would
      // then keep the provisional (empty) state forever.
      if (!this->get_adiabatic_conditions().is_initialized())
        return value_for_field(compute_state(position), n_comp);

      std::array<double, dim> key;
      for (unsigned int d = 0; d < dim; ++d)
        key[d] = position[d];

      auto entry = state_cache.find(key);
      if (entry == state_cache.end())
        entry = state_cache.emplace(key, compute_state(position)).first;

      return value_for_field(entry->second, n_comp);
    }



    template <int dim>
    void
    ThermodynamicEquilibrium<dim>::report () const
    {
      const unsigned int n_points = n_evaluated_points;
      const unsigned int n_solved = n_solved_points;
      const unsigned int n_background = n_kept_background_points;

      const double global_points = Utilities::MPI::sum(n_points, this->get_mpi_communicator());
      const double global_solved = Utilities::MPI::sum(n_solved, this->get_mpi_communicator());
      const double global_background = Utilities::MPI::sum(n_background, this->get_mpi_communicator());
      const double global_residual = Utilities::MPI::max(max_equilibrium_residual, this->get_mpi_communicator());
      const double global_sum_error = Utilities::MPI::max(max_composition_sum_error, this->get_mpi_communicator());
      const double global_phase_residual = Utilities::MPI::max(max_phase_residual, this->get_mpi_communicator());

      this->get_pcout() << "[TEQ initial condition] generated thermodynamic equilibrium initial state:" << std::endl
                        << "    specification = "
                        << (specification == Specification::bulk_composition ? "bulk composition" :
                            specification == Specification::melt_fraction ? "melt fraction" :
                            "phase compositions")
                        << ", mixing family = "
                        << (mixing_family == MixingFamily::tie_line ? "tie line" :
                            mixing_family == MixingFamily::background_to_liquid ? "background to liquid" :
                            mixing_family == MixingFamily::background_to_solid ? "background to solid" :
                            "background to component")
                        << ", pressure = adiabatic" << std::endl
                        << "    evaluated points             : " << global_points << std::endl
                        << "    points solved for a target   : " << global_solved << std::endl
                        << "    points kept at the background: " << global_background << std::endl
                        << "    max |Phi(c_bar) - phi*|      : " << global_residual << std::endl
                        << "    max |sum(c_bar) - 1|         : " << global_sum_error << std::endl;

      if (specification == Specification::phase_compositions)
        this->get_pcout() << "    max |phase compositions - prescribed| : "
                          << global_phase_residual << std::endl;

      if (dump_generated_state)
        dump_generated_state_to_file();
    }

    template <int dim>
    void
    ThermodynamicEquilibrium<dim>::dump_generated_state_to_file () const
    {
      const unsigned int n_fields = this->n_compositional_fields();
      const unsigned int stride = dim + n_fields;

      // Pack this process' nodes into a flat array and collect everything on
      // the first process.
      std::vector<double> local;
      local.reserve(state_cache.size() * stride);
      for (const auto &entry : state_cache)
        {
          for (unsigned int d = 0; d < dim; ++d)
            local.push_back(entry.first[d]);
          for (unsigned int c = 0; c < n_fields; ++c)
            local.push_back(value_for_field(entry.second, c));
        }

      const std::vector<std::vector<double>> gathered =
        Utilities::MPI::all_gather(this->get_mpi_communicator(), local);

      if (Utilities::MPI::this_mpi_process(this->get_mpi_communicator()) != 0)
        return;

      // The quadrature points of neighbouring cells differ in the last few
      // bits, so the coordinates have to be rounded before they can be used as
      // grid indices; otherwise the grid would contain near-duplicate rows.
      const double tolerance = 1e-9 * this->get_geometry_model().maximal_depth();
      const auto quantize = [tolerance] (const double value) -> double
      {
        return std::round(value / tolerance) * tolerance;
      };

      std::map<std::array<double, dim>, std::vector<double>> values;
      std::array<std::set<double>, dim> coordinates;

      for (const auto &chunk : gathered)
        for (unsigned int i = 0; i + stride <= chunk.size(); i += stride)
          {
            std::array<double, dim> point;
            for (unsigned int d = 0; d < dim; ++d)
              {
                point[d] = quantize(chunk[i + d]);
                coordinates[d].insert(point[d]);
              }
            values[point] = std::vector<double>(chunk.begin() + i + dim,
                                                chunk.begin() + i + stride);
          }

      std::array<unsigned int, dim> n_points;
      std::array<unsigned int, dim> strides;
      unsigned int total = 1;
      for (unsigned int d = 0; d < dim; ++d)
        {
          n_points[d] = coordinates[d].size();
          total *= n_points[d];
        }

      // The first coordinate (x) varies fastest, as in the ascii data files.
      strides[0] = 1;
      for (unsigned int d = 1; d < dim; ++d)
        strides[d] = strides[d - 1] * n_points[d - 1];

      std::ofstream file (output_file_name);
      AssertThrow(file, ExcMessage("Could not open '" + output_file_name + "' for writing."));

      file << "# Thermodynamic equilibrium initial condition generated by ASPECT\n";
      file << "# POINTS:";
      for (unsigned int d = 0; d < dim; ++d)
        file << " " << n_points[d];
      file << "\n";

      file << "# Columns:";
      const char *coordinate_names[3] = {"x", "y", "z"};
      for (unsigned int d = 0; d < dim; ++d)
        file << " " << coordinate_names[d];
      for (unsigned int c = 0; c < n_fields; ++c)
        file << " " << this->introspection().name_for_compositional_index(c);
      file << "\n";

      file.precision(17);
      for (unsigned int i = 0; i < total; ++i)
        {
          std::array<double, dim> point;
          for (unsigned int d = 0; d < dim; ++d)
            {
              const unsigned int index = (i / strides[d]) % n_points[d];
              point[d] = *std::next(coordinates[d].begin(), index);
            }

          const auto entry = values.find(point);
          if (entry == values.end())
            continue;

          for (unsigned int d = 0; d < dim; ++d)
            file << point[d] << " ";
          for (unsigned int c = 0; c < n_fields; ++c)
            file << entry->second[c] << (c + 1 == n_fields ? "\n" : " ");
        }

      this->get_pcout() << "[TEQ initial condition] wrote the generated state to "
                        << output_file_name << std::endl;
    }



    template <int dim>
    void
    ThermodynamicEquilibrium<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial composition model");
      {
        prm.enter_subsection ("Thermodynamic equilibrium");
        {
          prm.declare_entry ("Specification", "melt fraction",
                             Patterns::Selection ("bulk composition|melt fraction|phase compositions"),
                             "Select which quantity is prescribed. "
                             "'bulk composition': the background bulk composition given below is "
                             "used directly, and the melt fraction and the phase compositions "
                             "follow from the material model's equilibrium. "
                             "'melt fraction': a target melt fraction is prescribed and the bulk "
                             "composition is found by moving along the chosen one-dimensional "
                             "family in composition space until the equilibrium melt fraction "
                             "matches the target. "
                             "'phase compositions': the solid composition, the liquid "
                             "composition and the melt fraction are prescribed, and the bulk "
                             "composition follows from the phase mass balance. The result is "
                             "then checked against the material model's equilibrium.");

          prm.declare_entry ("Mixing family", "tie line",
                             Patterns::Selection ("tie line|background to liquid|background to solid|background to component"),
                             "The one-dimensional family in composition space along which the bulk "
                             "composition is moved when a target melt fraction is prescribed. "
                             "'tie line': between the solid and the liquid composition in "
                             "equilibrium with the background; along this family the mixing "
                             "fraction equals the melt fraction exactly. It is only defined where "
                             "the background is above the solidus. "
                             "'background to liquid' and 'background to solid': between the "
                             "background and the respective phase composition. "
                             "'background to component': between the background and a pure "
                             "component; use 'Family component index' to select it. This is the "
                             "family to use if melt is to be created below the solidus.");

          prm.declare_entry ("Family component index", "0",
                             Patterns::Integer (0),
                             "The index of the chemical component used by the mixing family "
                             "'background to component'.");

          prm.declare_entry ("Solving tolerance", "1e-14",
                             Patterns::Double (0.0),
                             "The tolerance to which the equilibrium is solved while the initial "
                             "condition is generated. This is independent of the material model "
                             "parameter 'Equilibrium solving tolerance' so that the initial "
                             "condition can be generated more accurately than the time stepping "
                             "without changing the latter.");

          prm.declare_entry ("Pressure source", "adiabatic",
                             Patterns::Selection ("adiabatic"),
                             "The pressure at which the equilibrium is evaluated while the initial "
                             "condition is generated. 'adiabatic' is the reference lithostatic "
                             "profile AdiabaticConditions::pressure(), which is exactly the "
                             "pressure that is written into the fluid pressure field at time "
                             "zero, so the initial condition and the model agree by construction.");

          prm.declare_entry ("Background composition model", "uniform",
                             Patterns::Selection ("uniform|depth profile"),
                             "How the background bulk composition is prescribed. "
                             "'uniform': a constant composition. "
                             "'depth profile': a one-column ASCII file, with one named column per "
                             "chemical component.");

          prm.declare_entry ("Uniform bulk composition", "",
                             Patterns::List (Patterns::Double (0.0)),
                             "The constant background bulk composition, one value per chemical "
                             "component, in the order given by the material model.");

          prm.declare_entry ("Target melt fraction model", "background",
                             Patterns::Selection ("background|uniform|depth profile|function"),
                             "How the target melt fraction is prescribed. "
                             "'background': the equilibrium melt fraction of the background "
                             "composition. "
                             "'uniform', 'depth profile' and 'function': a constant, a "
                             "one-column ASCII profile, or an expression.");

          prm.declare_entry ("Uniform target melt fraction", "0.0",
                             Patterns::Double (0.0, 1.0),
                             "The constant target melt fraction.");

          prm.declare_entry ("Add solitary wave", "false",
                             Patterns::Bool (),
                             "Whether to superimpose a solitary-wave-shaped melt fraction "
                             "anomaly on the background equilibrium. The wave is defined through "
                             "its dimensionless amplitude profile phi_r, which tends to one far "
                             "from the peak and equals the relative amplitude A at the peak, so "
                             "that the target melt fraction is f_bg * phi_r.");

          prm.declare_entry ("Peak position", "0.0",
                             Patterns::Double (),
                             "The position of the solitary wave peak, given as the coordinate "
                             "along the vertical direction of the model (y for a 2D box, i.e. "
                             "increasing upward). The wave is a horizontal band centred on this "
                             "position.");

          prm.declare_entry ("Peak porosity", "0.0",
                             Patterns::Double (0.0, 1.0),
                             "The absolute melt fraction at the peak of the solitary wave.");

          prm.declare_entry ("Truncation threshold", "1e-4",
                             Patterns::Double (0.0),
                             "The relative amplitude below which the solitary wave is truncated: "
                             "the wave is applied where phi_r - 1 is larger than this value. This "
                             "determines the width of the window around the peak, and therefore "
                             "the resolution that is required to represent the wave.");

          prm.declare_entry ("Derive compaction length", "true",
                             Patterns::Bool (),
                             "Whether the compaction length that sets the width of the solitary "
                             "wave is derived from the material model's permeability and "
                             "viscosities, evaluated at the background state of the peak. This is "
                             "the recommended choice: the wave then follows the model's own "
                             "constitutive law. The derived value is printed during the "
                             "initialization.");

          prm.declare_entry ("Compaction length", "0.0",
                             Patterns::Double (0.0),
                             "The compaction length used for the solitary wave, in m. Only used "
                             "if 'Derive compaction length' is false.");

          prm.declare_entry ("Report initial equilibrium", "true",
                             Patterns::Bool (),
                             "Whether to print a summary of the generated initial condition, "
                             "including the residual of the inverse solve, once the initial state "
                             "has been set.");

          Utilities::AsciiDataBase<dim>::declare_parameters (prm,
                                                             "$ASPECT_SOURCE_DIR/data/initial-composition/",
                                                             "background_bulk_composition.txt",
                                                             "Background depth profile");

          Utilities::AsciiDataBase<dim>::declare_parameters (prm,
                                                             "$ASPECT_SOURCE_DIR/data/initial-composition/",
                                                             "target_melt_fraction.txt",
                                                             "Target depth profile");

          prm.declare_entry ("Consistency check tolerance", "1e-8",
                             Patterns::Double (0.0),
                             "Only used for 'Specification = phase compositions'. The largest "
                             "difference between the prescribed melt fraction and phase "
                             "compositions and the values that the material model returns for the "
                             "generated bulk composition must not exceed this tolerance, otherwise "
                             "the behaviour selected by 'If inconsistent' is applied. The liquid "
                             "composition is only compared where the melt fraction is positive, "
                             "because below the solidus it is not defined.");

          prm.declare_entry ("If inconsistent", "error",
                             Patterns::Selection ("error|keep melt fraction|keep liquid"),
                             "Only used for 'Specification = phase compositions'. What to do if "
                             "the prescribed melt fraction and phase compositions are not "
                             "consistent with the material model's equilibrium. "
                             "'error': abort with a message that reports the residual. "
                             "'keep melt fraction': ignore the prescribed phase compositions and "
                             "find the bulk composition along the selected mixing family that "
                             "reproduces the prescribed melt fraction (this is the behaviour of "
                             "'Specification = melt fraction'). "
                             "'keep liquid': keep the prescribed liquid composition as the "
                             "enriched end member of the mixing family and solve for the melt "
                             "fraction along that family.");

          prm.declare_entry ("Solid phase composition model", "uniform",
                             Patterns::Selection ("uniform|depth profile"),
                             "Only used for 'Specification = phase compositions'. How the solid "
                             "phase composition is prescribed: a constant list of values or a "
                             "one-column ASCII file with one named column per chemical component.");

          prm.declare_entry ("Uniform solid composition", "",
                             Patterns::List (Patterns::Double (0.0)),
                             "The constant solid phase composition, one value per chemical "
                             "component, in the order given by the material model.");

          prm.declare_entry ("Liquid phase composition model", "uniform",
                             Patterns::Selection ("uniform|depth profile"),
                             "Only used for 'Specification = phase compositions'. How the liquid "
                             "phase composition is prescribed.");

          prm.declare_entry ("Uniform liquid composition", "",
                             Patterns::List (Patterns::Double (0.0)),
                             "The constant liquid phase composition, one value per chemical "
                             "component, in the order given by the material model.");

          prm.declare_entry ("Dump generated state to file", "false",
                             Patterns::Bool (),
                             "Whether to write the generated state to an ASCII file once the "
                             "initial condition has been set. The file uses the same layout as "
                             "the 'ascii data' initial composition model reads, so it can be "
                             "compared pointwise with an externally generated initial condition, "
                             "and read back in as a regression test.");

          prm.declare_entry ("Output file name", "generated_initial_composition.txt",
                             Patterns::FileName (),
                             "The name of the file written when 'Dump generated state to file' is "
                             "true. The file is written by the first MPI process only.");

          Utilities::AsciiDataBase<dim>::declare_parameters (prm,
                                                             "$ASPECT_SOURCE_DIR/data/initial-composition/",
                                                             "solid_phase_composition.txt",
                                                             "Solid phase depth profile");

          Utilities::AsciiDataBase<dim>::declare_parameters (prm,
                                                             "$ASPECT_SOURCE_DIR/data/initial-composition/",
                                                             "liquid_phase_composition.txt",
                                                             "Liquid phase depth profile");

          prm.enter_subsection ("Target function");
          {
            Functions::ParsedFunction<dim>::declare_parameters (prm, 1);
          }
          prm.leave_subsection ();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }



    template <int dim>
    void
    ThermodynamicEquilibrium<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial composition model");
      {
        const std::vector<std::string> model_names =
          Utilities::split_string_list (prm.get ("List of model names"));

        AssertThrow(model_names.size() == 1
                    && model_names[0] == "thermodynamic equilibrium",
                    ExcMessage("The initial composition model 'thermodynamic equilibrium' must be "
                               "the only entry in 'List of model names', because it generates all "
                               "compositional fields that take part in the equilibrium."));

        prm.enter_subsection ("Thermodynamic equilibrium");
        {
          const std::string specification_string = prm.get ("Specification");
          if (specification_string == "bulk composition")
            specification = Specification::bulk_composition;
          else if (specification_string == "melt fraction")
            specification = Specification::melt_fraction;
          else if (specification_string == "phase compositions")
            specification = Specification::phase_compositions;
          else
            AssertThrow(false, ExcNotImplemented());

          const std::string family_string = prm.get ("Mixing family");
          if (family_string == "tie line")
            mixing_family = MixingFamily::tie_line;
          else if (family_string == "background to liquid")
            mixing_family = MixingFamily::background_to_liquid;
          else if (family_string == "background to solid")
            mixing_family = MixingFamily::background_to_solid;
          else if (family_string == "background to component")
            mixing_family = MixingFamily::background_to_component;
          else
            AssertThrow(false, ExcNotImplemented());

          family_component = prm.get_integer ("Family component index");
          solving_tolerance = prm.get_double ("Solving tolerance");

          AssertThrow(prm.get ("Pressure source") == "adiabatic",
                      ExcNotImplemented());

          const std::string background_model_string = prm.get ("Background composition model");
          if (background_model_string == "uniform")
            background_model = FieldModel::uniform;
          else if (background_model_string == "depth profile")
            background_model = FieldModel::depth_profile;
          else
            AssertThrow(false, ExcNotImplemented());

          uniform_background = Utilities::string_to_double(
                                 Utilities::split_string_list (prm.get ("Uniform bulk composition")));

          const std::string target_model_string = prm.get ("Target melt fraction model");
          if (target_model_string == "background")
            target_model = FieldModel::background;
          else if (target_model_string == "uniform")
            target_model = FieldModel::uniform;
          else if (target_model_string == "depth profile")
            target_model = FieldModel::depth_profile;
          else if (target_model_string == "function")
            target_model = FieldModel::function;
          else
            AssertThrow(false, ExcNotImplemented());

          uniform_target = prm.get_double ("Uniform target melt fraction");

          const std::string solid_model_string = prm.get ("Solid phase composition model");
          if (solid_model_string == "uniform")
            solid_model = FieldModel::uniform;
          else if (solid_model_string == "depth profile")
            solid_model = FieldModel::depth_profile;
          else
            AssertThrow(false, ExcNotImplemented());

          const std::string liquid_model_string = prm.get ("Liquid phase composition model");
          if (liquid_model_string == "uniform")
            liquid_model = FieldModel::uniform;
          else if (liquid_model_string == "depth profile")
            liquid_model = FieldModel::depth_profile;
          else
            AssertThrow(false, ExcNotImplemented());

          if (prm.get ("Uniform solid composition") != "")
            uniform_solid = Utilities::string_to_double(
                              Utilities::split_string_list (prm.get ("Uniform solid composition")));

          if (prm.get ("Uniform liquid composition") != "")
            uniform_liquid = Utilities::string_to_double(
                               Utilities::split_string_list (prm.get ("Uniform liquid composition")));

          consistency_tolerance = prm.get_double ("Consistency check tolerance");

          const std::string inconsistent_string = prm.get ("If inconsistent");
          if (inconsistent_string == "error")
            inconsistent_behaviour = InconsistentBehaviour::error;
          else if (inconsistent_string == "keep melt fraction")
            inconsistent_behaviour = InconsistentBehaviour::keep_melt_fraction;
          else if (inconsistent_string == "keep liquid")
            inconsistent_behaviour = InconsistentBehaviour::keep_liquid;
          else
            AssertThrow(false, ExcNotImplemented());

          dump_generated_state = prm.get_bool ("Dump generated state to file");
          output_file_name = prm.get ("Output file name");

          use_solitary_wave = prm.get_bool ("Add solitary wave");
          peak_position = prm.get_double ("Peak position");
          peak_porosity = prm.get_double ("Peak porosity");
          truncation_threshold = prm.get_double ("Truncation threshold");
          derive_compaction_length = prm.get_bool ("Derive compaction length");
          prescribed_compaction_length = prm.get_double ("Compaction length");
          report_initial_equilibrium = prm.get_bool ("Report initial equilibrium");

          AssertThrow(use_solitary_wave == false || target_model == FieldModel::background,
                      ExcMessage("The solitary wave multiplies the background equilibrium melt "
                                 "fraction, so 'Target melt fraction model' has to be "
                                 "'background' when 'Add solitary wave' is true."));

          AssertThrow(use_solitary_wave == false || derive_compaction_length
                      || prescribed_compaction_length > 0.0,
                      ExcMessage("'Derive compaction length' is false, so a positive "
                                 "'Compaction length' has to be prescribed."));

          if (background_model == FieldModel::depth_profile)
            {
              background_profile = std::make_unique<Utilities::AsciiDataProfile<dim>>();
              background_profile->parse_parameters (prm, "Background depth profile");
            }

          if (target_model == FieldModel::depth_profile)
            {
              target_profile = std::make_unique<Utilities::AsciiDataProfile<dim>>();
              target_profile->parse_parameters (prm, "Target depth profile");
            }
          else if (target_model == FieldModel::function)
            {
              target_function = std::make_unique<Functions::ParsedFunction<dim>>(1);
              prm.enter_subsection ("Target function");
              target_function->parse_parameters (prm);
              prm.leave_subsection ();
            }

          if (solid_model == FieldModel::depth_profile)
            {
              solid_profile = std::make_unique<Utilities::AsciiDataProfile<dim>>();
              solid_profile->parse_parameters (prm, "Solid phase depth profile");
            }

          if (liquid_model == FieldModel::depth_profile)
            {
              liquid_profile = std::make_unique<Utilities::AsciiDataProfile<dim>>();
              liquid_profile->parse_parameters (prm, "Liquid phase depth profile");
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }
  }
}


// explicit instantiation
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(ThermodynamicEquilibrium,
                                              "thermodynamic equilibrium",
                                              "Generate an initial state that is in thermodynamic "
                                              "equilibrium by inverting the material model's own "
                                              "equilibrium map. Only the bulk composition is a "
                                              "free field; the melt fraction and the phase "
                                              "compositions are derived from it (they are "
                                              "prescribed fields that are recomputed at every "
                                              "time step anyway). A target melt fraction and a "
                                              "solitary-wave-shaped anomaly can be prescribed, "
                                              "in which case the bulk composition that realizes "
                                              "them is computed.")
  }
}
