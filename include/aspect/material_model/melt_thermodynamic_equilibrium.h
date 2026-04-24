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
        bool enable_chemical_reaction_rate;

        bool fill_debug_fields;
        
        // select a method to solve the equilibrium 
        std::string equilibrium_solving_method;

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
                             const double coefficient_B) const;

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
        virtual
        double
        find_solidus (const std::vector<double>& melting_points,
                      const std::vector<double>& bulk_concentrations,
                      const std::vector<double>& latent_heats,
                      const std::vector<double>& tuning_parameters) const;
        virtual
        double
        find_liquidus (const std::vector<double>& melting_points,
                       const std::vector<double>& bulk_concentrations,
                       const std::vector<double>& latent_heats,
                       const std::vector<double>& tuning_parameters) const;

        virtual
        double
        solve_eq_melt_fraction (const double temperature, 
                                const double pressure,
                                std::vector<double> bulk_concentrations) const;

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

          // // Check if the initial bounds are valid
          // if (f(a) * f(b) >= 0)
          // {
          //   // throw std::invalid_argument("The function must have opposite signs at the bounds.");
          //   AssertThrow(false,
          //               ExcMessage("The function must have opposite signs at the bounds."));
          //   return 0.5 * (a + b);
          // }

          double c = 0.0;
          for (unsigned int i = 0; i < max_iter; ++i)
          {
            // Perform checks during the iteration
            if (f(a) * f(b) > 0)
            {
              AssertThrow(false,
                    ExcMessage(
                      "Melt thermodynamic equilibrium (the enabled material model):"
                      "The new bounds do not make the function have opposite signs."
                    ));
              return 0.5 * (a + b);
            }
            c = (a + b) / 2.0;
            if (f(c) == 0.0)
            break;
            else if (f(c) * f(a) < 0)
            b = c;
            else
            a = c;
            if (fabs(f(c)) < tolerance)
            break;
          }
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
          partial_phi_partial_T(n_points, numbers::signaling_nan<double>())
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
    };

  } // namespace MaterialModel
  
} // namespace aspect

# endif // end of the new implementation


# ifndef ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS
// The original implementation, advecting concentration fields in phases.
#  pragma message("ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS is OFF, using the original implementation advecting concentration fields in phases.")

namespace aspect
{
  namespace MaterialModel
  {
    using namespace dealii;

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
         * @name Reference quantities
         * @{
         */
        double reference_darcy_coefficient () const override;


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

        // TODO: we'd better define a "n_components" variable 
        // as the number of virtual components
        // it should equals to a half of n_chemical_composition_fields()
        
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
        bool enable_chemical_reaction_rate;

        bool fill_debug_fields;
        
        // select a method to solve the equilibrium 
        std::string equilibrium_solving_method;

        // about chemical component name list
        unsigned int n_components;

        std::vector<std::string> component_names;

        std::string component_name_solid_suffix;
        std::string component_name_liquid_suffix;

        // we need a method to match the component index to composition indices
        std::vector<std::pair<unsigned int, unsigned int>> component_index_to_composition_indices;

        virtual
        bool
        match_component_index_to_composition_indices ();

        // quantities about melting process

        // firstly the melting lines for virtual compositions
        // T_m(P) = T_m_0 + A * P + B * P^2

        // T_m_0
        std::vector<double> melting_point_0_values; // in K
        // A
        std::vector<double> melting_curve_coefficient_A_values;
        // B
        std::vector<double> melting_curve_coefficient_B_values;

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
        double background_porosity;        
        
        virtual
        double
        temperature_melting (const double pressure,
                             const double temperature_m_0,
                             const double coefficient_A,
                             const double coefficient_B) const;

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
        virtual
        double
        find_solidus (const std::vector<double>& melting_points,
                      const std::vector<double>& bulk_concentrations,
                      const std::vector<double>& latent_heats,
                      const std::vector<double>& tuning_parameters) const;
        virtual
        double
        find_liquidus (const std::vector<double>& melting_points,
                       const std::vector<double>& bulk_concentrations,
                       const std::vector<double>& latent_heats,
                       const std::vector<double>& tuning_parameters) const;

        virtual
        double
        solve_eq_melt_fraction (const double temperature, 
                                const double pressure,
                                std::vector<double> bulk_concentrations) const;

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

          // // Check if the initial bounds are valid
          // if (f(a) * f(b) >= 0)
          // {
          //   // throw std::invalid_argument("The function must have opposite signs at the bounds.");
          //   AssertThrow(false,
          //               ExcMessage("The function must have opposite signs at the bounds."));
          //   return 0.5 * (a + b);
          // }

          double c = 0.0;
          for (unsigned int i = 0; i < max_iter; ++i)
          {
            // Perform checks during the iteration
            if (f(a) * f(b) > 0)
            {
              AssertThrow(false,
                    ExcMessage(
                      "Melt thermodynamic equilibrium (the enabled material model):"
                      "The new bounds do not make the function have opposite signs."
                    ));
              return 0.5 * (a + b);
            }
            c = (a + b) / 2.0;
            if (f(c) == 0.0)
            break;
            else if (f(c) * f(a) < 0)
            b = c;
            else
            a = c;
            if (fabs(f(c)) < tolerance)
            break;
          }
          return c;
        } 


    };

  } // namespace MaterialModel

} // namespace aspect

# endif // end of the original implementation

#endif
