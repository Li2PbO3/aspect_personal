/*
  Copyright (C) 2015 - 2023 by the authors of the ASPECT code.

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

#ifndef _aspect_material_model_melt_thermodynamic_equilibrium_h
#define _aspect_material_model_melt_thermodynamic_equilibrium_h

#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/postprocess/melt_statistics.h>
#include <aspect/melt.h>

#  pragma message("Compiling melt_thermodynamic_equilibrium.h")


# ifdef ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS
// New implementation, advecting bulk concentration fields.

#  pragma message("ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS is ON, using the new implementation advecting bulk concentration fields.")

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

    /**
     * Additional material model inputs carrying compositional fields sampled
     * from old_solution on the current cell.
     */
    template <int dim>
    class OldFieldStateInputs : public AdditionalMaterialInputs<dim>
    {
      public:
        OldFieldStateInputs (const unsigned int n_points,
                             const unsigned int n_comp)
          :
          old_porosity(n_points, numbers::signaling_nan<double>()),
          old_old_porosity(n_points, numbers::signaling_nan<double>()), // idle
          old_fields(n_points, std::vector<double>(n_comp, numbers::signaling_nan<double>())),
          old_old_fields(n_points, std::vector<double>(n_comp, numbers::signaling_nan<double>())) // idle
        {}

        std::vector<double> old_porosity;
        std::vector<double> old_old_porosity;
        std::vector<std::vector<double>> old_fields;
        std::vector<std::vector<double>> old_old_fields;

        void
        fill (const LinearAlgebra::BlockVector &solution,
              const FEValuesBase<dim>          &fe_values,
              const Introspection<dim>         &introspection) override
        {
          const unsigned int n_q_points = fe_values.n_quadrature_points;
          const unsigned int n_comp = introspection.n_compositional_fields;

          old_porosity.resize(n_q_points);
          old_fields.assign(n_q_points, std::vector<double>(n_comp));
          old_old_porosity.resize(n_q_points); // idle
          old_old_fields.assign(n_q_points, std::vector<double>(n_comp)); // idle

          const unsigned int porosity_idx = introspection.compositional_index_for_name("porosity");

          for (unsigned int c = 0; c < n_comp; ++c)
            {
              std::vector<double> temp_field(n_q_points);
              fe_values[introspection.extractors.compositional_fields[c]]
              .get_function_values(solution, temp_field);

              for (unsigned int q = 0; q < n_q_points; ++q)
                {
                  old_fields[q][c] = temp_field[q];
                  if (c == porosity_idx)
                    old_porosity[q] = temp_field[q];
                }
            }
        }
    };

    /**
     * Additional material model inputs carrying pointwise compositional
     * gradients and velocity divergence sampled from the provided solution
     * vector.
     */
    template <int dim>
    class CompositionFieldGradientsInputs : public AdditionalMaterialInputs<dim>
    {
      public:
        CompositionFieldGradientsInputs (const unsigned int n_comp)
          :
          composition_gradients(1, std::vector<Tensor<1,dim>>(n_comp)),
          velocity_divergence(1, numbers::signaling_nan<double>())
        {}

        std::vector<std::vector<Tensor<1,dim>>> composition_gradients;
        std::vector<double> velocity_divergence;

        // [P3-6] Explicit marker telling the caller whether the gradients below were
        // actually sampled from the solution. The size of composition_gradients is
        // always n_q_points (because of the assign() in fill()), so a size check
        // cannot distinguish "filled with real values" from "left at zero because
        // this FEValues object was not set up with update_gradients".
        bool gradients_available = false;

        void
        fill (const LinearAlgebra::BlockVector &solution,
              const FEValuesBase<dim>          &fe_values,
              const Introspection<dim>         &introspection) override
        {
          const unsigned int n_q_points = fe_values.n_quadrature_points;
          const unsigned int n_comp = introspection.n_compositional_fields;
          composition_gradients.assign(n_q_points, std::vector<Tensor<1,dim>>(n_comp));
          velocity_divergence.assign(n_q_points, numbers::signaling_nan<double>());
          gradients_available = false;

          if (((fe_values.get_update_flags() & update_gradients) == update_default)
              || n_q_points == 0)
            return;

          fe_values[introspection.extractors.velocities]
          .get_function_divergences(solution, velocity_divergence);

          for (unsigned int c = 0; c < n_comp; ++c)
            {
              std::vector<Tensor<1,dim>> grad_values(n_q_points);
              fe_values[introspection.extractors.compositional_fields[c]]
              .get_function_gradients(solution, grad_values);

              for (unsigned int q = 0; q < n_q_points; ++q)
                composition_gradients[q][c] = grad_values[q];
            }

          gradients_available = true;
        }
    };

    /**
     * Additional material model inputs carrying fluid pressure sampled
     * from the provided solution vector.
     */
    template <int dim>
    class FluidPressureInputs : public AdditionalMaterialInputs<dim>
    {
      public:
        FluidPressureInputs (const unsigned int n_points)
          :
          fluid_pressure(n_points, numbers::signaling_nan<double>())
        {}

        std::vector<double> fluid_pressure;

        void
        fill (const LinearAlgebra::BlockVector &solution,
              const FEValuesBase<dim>          &fe_values,
              const Introspection<dim>         &introspection) override
        {
          const unsigned int n_q_points = fe_values.n_quadrature_points;
          fluid_pressure.resize(n_q_points, numbers::signaling_nan<double>());

          if (!introspection.variable_exists("fluid pressure"))
            return;

          const FEValuesExtractors::Scalar ex_p_f = introspection.variable("fluid pressure").extractor_scalar();
          fe_values[ex_p_f].get_function_values(solution, fluid_pressure);
        }
    };

    /**
     * A material model that implements a simple formulation of the
     * material parameters required for the modeling of melt transport
     * in a global model, including a source term for the porosity according
     * a simplified linear melting model.
     *
     * The model is considered incompressible, following the definition
     * described in Interface::is_compressible.
     *
     * @ingroup MaterialModels
     */
    template <int dim>
    class MeltThermodynamicEquilibrium : public MaterialModel::MeltInterface<dim>,
      public MaterialModel::MeltFractionModel<dim>,
      public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Return whether the model is compressible or not.  Incompressibility
         * does not necessarily imply that the density is constant; rather, it
         * may still depend on temperature or pressure. In the current
         * context, compressibility means whether we should solve the continuity
         * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
         * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).
         */
        bool is_compressible () const override;

        void evaluate(const typename Interface<dim>::MaterialModelInputs &in,
                      typename Interface<dim>::MaterialModelOutputs &out) const override;

        /**
         * Compute the equilibrium melt fractions for the given input conditions.
         * @p in and @p melt_fractions need to have the same size.
         *
         * @param in Object that contains the current conditions.
         * @param melt_fractions Vector of doubles that is filled with the
         * equilibrium melt fraction for each given input conditions.
         */
        void melt_fractions (const MaterialModel::MaterialModelInputs<dim> &in,
                             std::vector<double> &melt_fractions) const override;

        /**
         * Return the pressure at which the thermodynamic equilibrium should be
         * evaluated for evaluation point @p q of @p in.
         *
         * Which of the available pressure fields is used is selected by the
         * run-time parameter "Pressure for thermodynamic equilibrium":
         *
         * - 'fluid pressure': the melt (fluid) pressure p_f.  This is the
         *   physically consistent two-phase choice, but in a model with melt
         *   transport p_f carries the dynamic pressure and the compaction
         *   pressure, so the thermodynamic equilibrium inherits the
         *   compaction-wave structure.
         * - 'adiabatic pressure': the reference lithostatic profile
         *   AdiabaticConditions::pressure(), integrated downward from the
         *   surface pressure using the laterally averaged reference density.
         *   It is a smooth, laterally uniform function of depth only, so the
         *   equilibrium state becomes a well-defined function of
         *   (depth, T, c) and is free of any compaction-wave imprint.
         * - 'solid pressure': the (total) solid pressure p, i.e. the pressure
         *   variable of the Stokes system.
         *
         * Whenever the adiabatic profile is not available yet (it is built
         * during the first time step by evaluating the material model itself),
         * this function falls back to the pressure passed in @p in.
         */
        double
        equilibrium_pressure (const MaterialModel::MaterialModelInputs<dim> &in,
                              const unsigned int                         q,
                              const FluidPressureInputs<dim>            *fluid_pressure_input) const;

        /**
         * Return the chemical components for which this model computes a
         * solid-liquid thermodynamic equilibrium. See the documentation of the
         * base class.
         */
        std::vector<MaterialModel::EquilibriumComponent>
        get_equilibrium_components () const override;

        /**
         * Evaluate the equilibrium state for a single point. See the
         * documentation of the base class. This is a thin wrapper around
         * solve_eq_melt_fraction() and calculate_concentration_solid/liquid(),
         * i.e. it uses exactly the same thermodynamics as the time stepping.
         *
         * The temperature is expected in K and is converted to degree Celsius
         * internally, as everywhere else in this model.
         */
        bool
        evaluate_equilibrium_state (const double               pressure,
                                    const double               temperature,
                                    const std::vector<double> &bulk_composition,
                                    double                    &melt_fraction,
                                    std::vector<double>       &solid_composition,
                                    std::vector<double>       &liquid_composition,
                                    const double               tolerance = -1.0) const override;

        /**
         * @name Reference quantities
         * @{
         */
        double reference_darcy_coefficient () const override;

        // getter methods for component latent heating models
        unsigned int get_n_chemical_components () const;
        const std::vector<std::string> & get_chemical_component_names () const;
        const std::vector<double> & get_latent_heat_values () const;

        /**
         * @}
         */

        /**
         * @name Functions used in dealing with run-time parameters
         * @{
         */
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
        
        // this material model needs to fill old solution fields.
        // so we need to override this function to fill the additional material model inputs
        void
        fill_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &input,
                      const LinearAlgebra::BlockVector        &solution,
                      const FEValuesBase<dim>                 &fe_values,
                      const Introspection<dim>                &introspection) const override;

        /**
         * @}
         */

        void
        create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const override;

      private:

        const double ZERO_CELSIUS_IN_KELVIN = 273.15;

        /**
         * we need some new variable now.
         */
        
        // reference quantities
        double reference_rho_s;
        double reference_rho_f;
        double reference_T;
        double eta_0;
        double xi_0;
        double eta_f;
        double reference_permeability;
        double reference_specific_heat;
        
        // quantities that control the quantities 
        // changing with p, T, and phi
        double thermal_viscosity_exponent;
        double thermal_bulk_viscosity_exponent;
        double viscosity_activation_energy;
        double viscosity_activation_volume;
        double thermal_expansivity;
        double alpha_phi;
        double compressibility;
        double melt_compressibility;
        
        // quantities that do not change
        double thermal_conductivity;       
        bool include_melting_and_freezing;
        double melting_time_scale;

        // switch for equilibrium calculation
        // if it's false, then the model degrades to a simple melt transport model
        // without thermodynamic equilibrium calculation and melting/freezing source term
        bool enable_equilibrium_calculation;
        // [P3-1] default-initialized; the corresponding parameter is now parsed as well
        bool enable_chemical_reaction_rate = false;

        bool fill_debug_fields;
        bool fill_prescribed_melting_rate_field;
        
        // select a method to solve the equilibrium 
        std::string equilibrium_solving_method;

        /**
         * Selects which pressure field is used as the pressure of the
         * thermodynamic equilibrium calculation.
         */
        enum class EquilibriumPressure
        {
          fluid_pressure,
          adiabatic_pressure,
          solid_pressure
        };

        /**
         * The pressure field used for the thermodynamic equilibrium
         * calculation.  See equilibrium_pressure().
         */
        EquilibriumPressure equilibrium_pressure_choice = EquilibriumPressure::adiabatic_pressure;

        // Convergence tolerance of the bisection used for the equilibrium
        // calculation. It is a tolerance on the *residual of the equilibrium
        // equation*, not on the melt fraction itself; the corresponding error
        // in the liquid concentrations is roughly tolerance / |dF/df|, which
        // for these parameters is of the order 1e-10 for the default 1e-10.
        // Tightening it reduces the spurious spatial variation of c_l (and
        // hence the q.grad(c_l) term of the melting-rate estimator).
        double equilibrium_tolerance;

        // about chemical component name list
        unsigned int n_components;

        std::vector<std::string> component_names;
        std::string component_name_solid_suffix;
        std::string component_name_liquid_suffix;
        
        // we need a method to match the component index to composition index
        struct ComponentIndicesTriplet
        {
          unsigned int bulk_index;
          unsigned int solid_index;
          unsigned int liquid_index;
        };
        std::vector<ComponentIndicesTriplet> component_indices_list;

        virtual
        bool
        match_component_index_to_composition_index ();

        // quantities about melting process

        // firstly the melting lines for virtual compositions
        // T_m(P) = T_m_0 + A * P + B * P^2

        // T_m_0
        std::vector<double> melting_point_0_values; // in K
        // A
        std::vector<double> melting_curve_coefficient_A_values;
        // B
        std::vector<double> melting_curve_coefficient_B_values;
        // pressure threshold
        std::vector<double> melting_curve_pressure_thresholds;

        // secondly the quantities in equilibrium constant
        // K = c_sol / c_liq
        //   = exp((L / r) * (1 / T - 1 / T_m))

        // L
        std::vector<double> latent_heat_values;
        // r
        std::vector<double> tuning_parameter_values;
        // The physical meaning of this tuning parameter is not clear enough,
        // and it is roughly the "effective gas constant" of a certain virtual component.

        // background porosity to avoid zero permeability
        // TODO: maybe useless...
        double background_porosity;
        
        virtual
        double
        temperature_melting (const double pressure,
                             const double temperature_m_0,
                             const double coefficient_A,
                             const double coefficient_B,
                             const double pressure_threshold) const;

        virtual
        double
        equilibrium_constant (/* const double pressure, */
                              const double temperature,
                              const double latent_heat,
                              const double tuning_parameter,
                              const double _T_m) const;
        // we already use the pressure while calculating the melting point
        // so we don't need to pass it again
        
        // seems like it's going to be helpful to declare a pair of functions
        // to calculate the solidus and liquidus temperature
        // [P1-2] find_solidus()/find_liquidus() were removed: the solid/liquid
        // branch is now decided from the endpoint signs of the melt-fraction
        // equation, which is algebraically equivalent and much cheaper.

        virtual
        double
        solve_eq_melt_fraction (const double temperature,
                                const double pressure,
                                const std::vector<double> &bulk_concentrations,
                                const double tolerance = -1.0,
                                // [P2-2] optional outputs so callers can reuse the
                                // melting points / equilibrium constants instead of
                                // recomputing them with extra exp() evaluations
                                std::vector<double> *melting_points_out = nullptr,
                                std::vector<double> *equilibrium_constants_out = nullptr) const;

        virtual
        double
        calculate_concentration_solid (const double c_bulk,
                                       const double f, // melt fraction
                                       const double eq_const // equilibrium constant
                                       ) const;
        virtual
        double
        calculate_concentration_liquid (const double c_bulk,
                                       const double f, // melt fraction
                                       const double eq_const // equilibrium constant
                                       ) const;

        double
        pressure_temperature_viscosity_factor (const double temperature,
                       const double pressure) const;

        // i decide to define my own solvers to find roots
        // of equations, because i have no idea how to use the root finding method in deal.ii

        // this is a bisection method
        double
        bisection (const std::function<double(const double)> &f,
              const double lower_bound,
              const double upper_bound,
              const unsigned int max_iter = 1000,
              const double tolerance = 1e-10) const
        {
          double a = lower_bound;
          double b = upper_bound;

          // [P1-1] Evaluate f at the two bounds once and cache the value at the
          // current midpoint. The previous version recomputed f(a), f(b) once each
          // and f(c) three times per iteration (6 evaluations per loop), which the
          // compiler cannot eliminate because f is a std::function.
          const double fa0 = f(a);
          const double fb0 = f(b);
          AssertThrow(fa0 * fb0 <= 0,
                      ExcMessage("Melt thermodynamic equilibrium (the enabled material model): "
                                 "the function must have opposite signs at the bounds. "
                                 "a=" + std::to_string(a) + " b=" + std::to_string(b)));
          double fa = fa0;
          double fb = fb0;
          double c = 0.0;
          double fc = 0.0;
          for (unsigned int i = 0; i < max_iter; ++i)
            {
              c = 0.5 * (a + b);
              fc = f(c);
              if (fc == 0.0)
                break;
              else if (fc * fa < 0)
                {
                  b = c;
                  fb = fc;
                }
              else
                {
                  a = c;
                  fa = fc;
                }
              if (std::fabs(fc) < tolerance)
                break;
            }
          (void)fb;
          return c;
        }


    }; // class MeltThermodynamicEquilibrium

    /**
     * A class derived from AdditionalMaterialOutputs 
     * that contains the outputs related to component exchange between solid and liquid.
     * This output is for heating models calculating the heat source term 
     * related to component exchange, so called component latent heat.
     * @ingroup MaterialModels
     */
    template <int dim>
    class ComponentPhaseExchangeOutputs : public AdditionalMaterialOutputs<dim>
    {
      public:
        ComponentPhaseExchangeOutputs (const unsigned int n_points)
          :
          effective_latent_heat(n_points, numbers::signaling_nan<double>()),
          partial_phi_partial_T(n_points, numbers::signaling_nan<double>()),
          melting_rate(n_points, numbers::signaling_nan<double>())
        {}

        /**
         * L_eff (effective latent heat) is the latent heat associated with the component exchange between solid and liquid.
         * Latent heat for different components may be different, but the temperature equation can't deal with 
         * multiple latent heats. So we need a effective latent heat.
         */
        std::vector<double> effective_latent_heat; // [n_points]
        /**
         * We also need partial phi partial T.
         */
        std::vector<double> partial_phi_partial_T; // [n_points]
        /**
         * if we choose to apply the phase exchange heating in a source term form,
         * we will use this melting_rate value.
         */
        // here the melting_rate is the changing rate of melt fraction or porosity.
        // this means the melting_rate is within the unit of 1/s or 1/yr
        std::vector<double> melting_rate; // [n_points]
    };

  } // namespace MaterialModel
  
} // namespace aspect

# endif // ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS

#endif
