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

#ifndef _aspect_initial_composition_thermodynamic_equilibrium_h
#define _aspect_initial_composition_thermodynamic_equilibrium_h

#include <aspect/initial_composition/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

#include <array>
#include <map>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * Initial compositional fields that are consistent with the thermodynamic
     * equilibrium of the material model.
     *
     * In a melt model only the bulk composition is a real degree of freedom:
     * the melt fraction (porosity), the solid and the liquid concentrations
     * are prescribed fields that are recomputed from the equilibrium map at
     * every time step.  This plugin therefore only has to produce the bulk
     * composition, and it does so by inverting the material model's own
     * equilibrium map at the grid nodes, using exactly the same thermodynamics
     * as the time stepping.  This removes the parameter-rounding, ASCII-file
     * interpolation and float32 inconsistencies of an externally generated
     * initial condition.
     *
     * Two ways of specifying the state are implemented:
     *
     * - `bulk composition`: the bulk composition is given directly (uniform,
     *   depth profile or function) and the melt fraction and phase
     *   compositions follow from the equilibrium.
     * - `melt fraction`: a target melt fraction is given and the bulk
     *   composition is found by moving along a one-dimensional family in
     *   composition space until the equilibrium melt fraction matches the
     *   target.  For the `tie line` family this can be done analytically.
     *
     * A solitary-wave-like target melt fraction can be superimposed on the
     * background equilibrium; its relative amplitude, its compaction length
     * and its truncation window are all derived from the material model
     * rather than being user parameters (see the documentation of the
     * individual parameters).
     *
     * @ingroup InitialCompositions
     */
    template <int dim>
    class ThermodynamicEquilibrium : public Interface<dim>,
      public aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        ThermodynamicEquilibrium ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize () override;

        /**
         * Return the initial composition of compositional field @p n_comp at
         * position @p position.
         */
        double
        initial_composition (const Point<dim> &position,
                             const unsigned int n_comp) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Which quantity the user specifies.
         */
        enum class Specification
        {
          bulk_composition,
          melt_fraction,
          phase_compositions
        };

        /**
         * What to do if the prescribed phase compositions and melt fraction are
         * not consistent with the material model's equilibrium.
         */
        enum class InconsistentBehaviour
        {
          error,
          keep_melt_fraction,
          keep_liquid
        };

        /**
         * The one-dimensional family in composition space along which the bulk
         * composition is moved to reach the target melt fraction.
         */
        enum class MixingFamily
        {
          tie_line,
          background_to_liquid,
          background_to_solid,
          background_to_component
        };

        /**
         * How a scalar field (bulk composition or target melt fraction) is
         * prescribed.
         */
        enum class FieldModel
        {
          uniform,
          depth_profile,
          function,
          background
        };

        /**
         * The full thermodynamic state at one point.
         */
        struct EquilibriumState
        {
          double pressure = 0.0;
          double temperature = 0.0;
          std::vector<double> bulk_composition;
          double melt_fraction = 0.0;
          std::vector<double> solid_composition;
          std::vector<double> liquid_composition;
        };

        Specification specification;
        MixingFamily  mixing_family;

        FieldModel background_model;
        FieldModel target_model;
        FieldModel solid_model;
        FieldModel liquid_model;

        std::vector<double> uniform_background;
        std::vector<double> uniform_solid;
        std::vector<double> uniform_liquid;
        double uniform_target;
        unsigned int family_component;

        double consistency_tolerance;
        InconsistentBehaviour inconsistent_behaviour;

        bool   use_solitary_wave;
        double peak_position;
        double peak_porosity;
        double truncation_threshold;
        bool   derive_compaction_length;
        double prescribed_compaction_length;

        double solving_tolerance;
        bool   report_initial_equilibrium;
        bool   dump_generated_state;
        std::string output_file_name;

        /**
         * The chemical components of the material model and the indices of the
         * corresponding compositional fields.
         */
        std::vector<MaterialModel::EquilibriumComponent> components;
        std::vector<unsigned int> bulk_indices;
        std::vector<unsigned int> solid_indices;
        std::vector<unsigned int> liquid_indices;
        unsigned int porosity_index;

        /**
         * A shared pointer to the manager of the initial temperature model.
         *
         * The simulator releases its own pointer after the first time step,
         * but the boundary composition model keeps calling this plugin through
         * the initial composition manager, so we have to keep the temperature
         * initial condition alive ourselves.
         */
        std::shared_ptr<const InitialTemperature::Manager<dim>> initial_temperature_manager;

        /**
         * Models for the background bulk composition.
         */
        std::unique_ptr<Utilities::AsciiDataProfile<dim>> background_profile;
        std::unique_ptr<Functions::ParsedFunction<dim>>    background_function;

        /**
         * Models for the target melt fraction.
         */
        std::unique_ptr<Utilities::AsciiDataProfile<dim>> target_profile;
        std::unique_ptr<Functions::ParsedFunction<dim>>    target_function;

        /**
         * Models for the prescribed phase compositions (only used for
         * 'Specification = phase compositions').
         */
        std::unique_ptr<Utilities::AsciiDataProfile<dim>> solid_profile;
        std::unique_ptr<Utilities::AsciiDataProfile<dim>> liquid_profile;

        /**
         * Quantities that are derived from the material model once, on the
         * first evaluation for which the adiabatic reference profile exists.
         */
        mutable bool notified_about_reference_profile;
        mutable bool derived_parameters_are_ready;
        mutable bool solitary_wave_is_possible;
        mutable double peak_background_melt_fraction;
        mutable double relative_amplitude;
        mutable double compaction_length;
        mutable double window_radius;

        /**
         * A cache of the computed states.
         *
         * The initial composition manager asks this plugin for every
         * compositional field in turn and traverses the whole mesh for each of
         * them, so without a cache the same state would be computed once per
         * compositional field (eleven times in the melt models). The cache also
         * holds the full node set, which is what `Dump generated state to file'
         * writes out.
         */
        mutable std::map<std::array<double, dim>, EquilibriumState> state_cache;

        /**
         * Diagnostics accumulated while the initial condition is evaluated.
         */
        mutable double max_equilibrium_residual;
        mutable double max_phase_residual;
        mutable double max_composition_sum_error;
        mutable unsigned int n_evaluated_points;
        mutable unsigned int n_solved_points;
        mutable unsigned int n_kept_background_points;

        /**
         * Evaluate the equilibrium of the material model at one point, and
         * abort if the model does not implement it.
         */
        void
        evaluate_equilibrium (const double               pressure,
                              const double               temperature,
                              const std::vector<double> &bulk_composition,
                              double                    &melt_fraction,
                              std::vector<double>       &solid_composition,
                              std::vector<double>       &liquid_composition) const;

        /**
         * The background bulk composition at a point (not normalized).
         */
        std::vector<double>
        background_bulk_composition (const Point<dim> &position) const;

        /**
         * The base target melt fraction (before the solitary wave is applied).
         */
        double
        base_target_melt_fraction (const Point<dim> &position,
                                   const double      background_melt_fraction) const;

        /**
         * The target melt fraction at a point, including the solitary wave if
         * it is enabled.
         */
        double
        target_melt_fraction (const Point<dim> &position,
                              const double      background_melt_fraction) const;

        /**
         * Compute the derived solitary-wave parameters from the state at the
         * peak and the material model's permeability and viscosities.
         */
        void
        compute_derived_parameters () const;

        /**
         * The dimensionless solitary-wave profile phi_r(x/delta), with
         * phi_r -> 1 for large arguments and phi_r(0) = A.
         */
        double
        solitary_wave_profile (const double scaled_distance) const;

        /**
         * The dimensionless profile that corresponds to a given scaled
         * distance, i.e. the inverse of the implicit function.
         */
        double
        solitary_wave_implicit_function (const double relative_porosity) const;

        /**
         * The prescribed solid phase composition at a point.
         */
        std::vector<double>
        solid_phase_composition (const Point<dim> &position) const;

        /**
         * The prescribed liquid phase composition at a point.
         */
        std::vector<double>
        liquid_phase_composition (const Point<dim> &position) const;

        /**
         * The value that the plugin returns for compositional field @p n_comp.
         */
        double
        value_for_field (const EquilibriumState &state,
                         const unsigned int      n_comp) const;

        /**
         * The full state at a point, including the inverse solve for the bulk
         * composition if a target melt fraction is prescribed.
         */
        EquilibriumState
        compute_state (const Point<dim> &position) const;

        /**
         * Print a summary of the generated initial condition.
         */
        void
        report () const;

        /**
         * Write the generated state to an ASCII file in the format that the
         * 'ascii data' initial composition model reads, so that the internal
         * and the external way of generating an initial condition can be
         * compared pointwise.
         */
        void
        dump_generated_state_to_file () const;
    };
  }
}


#endif
