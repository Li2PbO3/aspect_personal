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
# include <aspect/simulator.h>
namespace
{
  const bool local_debug = true;
}

// TODO: this->introspection().compositional_index_for_name(component_name)
// would meet exceptions if the compositional field with the given name is not found. 
// For now we catch the exception and set the index to invalid_unsigned_int 
// to avoid the program from crashing.
// However there is a method compositional_name_exists() under struct Introspection.
// We'd better use this method to check the existence of the compositional field
// before calling compositional_index_for_name() to avoid the exception.

# ifdef ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS
// New implementation, advecting bulk concentration fields.

#  pragma message("ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS is ON, using the new implementation advecting bulk concentration fields.")

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

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    pressure_temperature_viscosity_factor (const double temperature,
                                           const double pressure) const
    {
      const double universal_gas_constant = 8.31446261815324;
      const double clamped_temperature = std::max(temperature, 1.0);
      const double clamped_reference_temperature = std::max(reference_T, 1.0);
      const double local_pressure = std::max(0.0, pressure);
      const double reference_pressure = std::max(0.0, this->get_surface_pressure());

      const double log_viscosity_reference =
        (viscosity_activation_energy + viscosity_activation_volume * reference_pressure)
        / (universal_gas_constant * clamped_reference_temperature);
      const double log_viscosity_local =
        (viscosity_activation_energy + viscosity_activation_volume * local_pressure)
        / (universal_gas_constant * clamped_temperature);

      const double bounded_log_factor =
        std::max(std::min(log_viscosity_local - log_viscosity_reference, 50.0), -50.0);

      return std::exp(bounded_log_factor);
    }

    // match component index to composition index
    template <int dim>
    bool
    MeltThermodynamicEquilibrium<dim>::
    match_component_index_to_composition_index ()
    {
      // if equilibrium calculation is disabled, return true directly
      if (!enable_equilibrium_calculation)
        return true;

      component_indices_list.clear();
      for (unsigned int component_idx = 0; component_idx < n_components; ++component_idx)
        {
          const std::string component_name = component_names[component_idx];

          const std::string bulk_field_name = component_name;
          bool bulk_index_found = false;
          unsigned int bulk_field_idx = numbers::invalid_unsigned_int;
          try
          {
            bulk_field_idx = this->introspection().compositional_index_for_name(component_name);
          }
          catch (const std::exception &e)
          {
            AssertThrow(false,
                        ExcMessage("Melt thermodynamic equilibrium (the enabled material model):"
                                   "Error while trying to find the compositional field for the bulk concentration of component '" + component_name + "'. "
                                   "Please make sure there is a compositional field with the name '" + component_name + "' defined in the 'Compositional fields' section of the parameter file."));
            return false;
          }          
          bulk_index_found = (bulk_field_idx != numbers::invalid_unsigned_int);

          const std::string solid_field_name = component_name + component_name_solid_suffix;
          bool solid_index_found = false;
          unsigned int solid_field_idx = numbers::invalid_unsigned_int;
          try
          {
            solid_field_idx = this->introspection().compositional_index_for_name(solid_field_name);
          }
          catch (const std::exception &e)
          {
            solid_index_found = false;
            solid_field_idx = numbers::invalid_unsigned_int;
            // AssertThrow(false,
            //             ExcMessage("Melt thermodynamic equilibrium (the enabled material model):"
            //                        "Error while trying to find the compositional field for the solid concentration of component '" + component_name + "'. "
            //                        "Please make sure there is a compositional field with the name '" + solid_field_name + "' defined in the 'Compositional fields' section of the parameter file."));
            // return false;
          }
          solid_index_found = (solid_field_idx != numbers::invalid_unsigned_int);

          const std::string liquid_field_name = component_name + component_name_liquid_suffix;
          bool liquid_index_found = false;
          unsigned int liquid_field_idx = numbers::invalid_unsigned_int;
          try
          {
            liquid_field_idx = this->introspection().compositional_index_for_name(liquid_field_name);
          }
          catch (const std::exception &e)
          {
            liquid_index_found = false;
            liquid_field_idx = numbers::invalid_unsigned_int;
            AssertThrow(false,
                        ExcMessage("Melt thermodynamic equilibrium (the enabled material model):"
                                   "Error while trying to find the compositional field for the liquid concentration of component '" + component_name + "'. "
                                   "Please make sure there is a compositional field with the name '" + liquid_field_name + "' defined in the 'Compositional fields' section of the parameter file."));
            return false;
          }
          liquid_index_found = (liquid_field_idx != numbers::invalid_unsigned_int);

          if (bulk_index_found)
            {
              ;
            }
          else
            {
              // this->get_pcout() << "Component '" << component_name 
              //                   << "': bulk_found=" << bulk_found 
              //                   << ", liquid_found=" << liquid_found << std::endl;              
              return false;
            }
          ComponentIndicesTriplet indices;
          indices.bulk_index = bulk_field_idx;
          indices.solid_index = solid_field_idx;
          indices.liquid_index = liquid_field_idx;
          component_indices_list.push_back(indices);
        }
      return true;
    }

    template <int dim>
    unsigned int MeltThermodynamicEquilibrium<dim>::get_n_chemical_components() const
    {
      return n_components;
    }

    template <int dim>
    const std::vector<std::string> & MeltThermodynamicEquilibrium<dim>::get_chemical_component_names() const
    {
      return component_names;
    }

    template <int dim>
    const std::vector<double> & MeltThermodynamicEquilibrium<dim>::get_latent_heat_values() const
    {
      return latent_heat_values;
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
      // threshold at 6 GPa
      const double pressure_threshold = 6.0e9; // 6 GPa
      double T_m = 0.0;
      if (pressure <= pressure_threshold)
        {
          T_m = temperature_m_0 + coefficient_A * pressure + coefficient_B * pressure * pressure;
        }
      else
        {
          const double T_threshold = temperature_m_0 + coefficient_A * pressure_threshold + coefficient_B * pressure_threshold * pressure_threshold;
          const double dT_dP_threshold = coefficient_A + 2.0 * coefficient_B * pressure_threshold;
          // extrapolate linearly with the slope at the threshold pressure
          T_m = T_threshold + dT_dP_threshold * (pressure - pressure_threshold);
        }
      return T_m;
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
      const double temperature_epsilon = 0.0; // 1e-4; // 0.0001 degree Celsius
      const bool is_all_solid = (temperature - solidus < temperature_epsilon);
      const bool is_all_liquid = (temperature - liquidus > -temperature_epsilon);
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
      
      // // we need to check if the root is between 0 and 1
      // AssertThrow(inner_melt_fraction >= 0.0 && inner_melt_fraction <= 1.0,
      //   ExcMessage("The root finder find the melt fraction. "
      //              "at T = " + std::to_string(temperature) + " deg C, and P = " + std::to_string(pressure) + " Pa. "
      //              "In this condition, the solidus is " + std::to_string(solidus) + " deg C, and the liquidus is " + std::to_string(liquidus) + " deg C. "
      //               + "However The melt fraction is not between 0 and 1: " + std::to_string(inner_melt_fraction)));

      // we need to check if the root is between 0 and 1
      // However we don't throw exceptions but just give warnings and set the melt fraction to 0 or 1
      if (inner_melt_fraction < 0.0)
        {
          // this->get_pcout() << "[MM eval] Warning: the melt fraction is smaller than 0: " << inner_melt_fraction << " at T = " << temperature << " deg C, and P = " << pressure << " Pa. "
          //                   << "In this condition, the solidus is " << solidus << " deg C, and the liquidus is " << liquidus << " deg C. "
          //                   << "Setting the melt fraction to 0.\n";
          inner_melt_fraction = 0.0;
        }
      if (inner_melt_fraction > 1.0)
        {
          // this->get_pcout() << "[MM eval] Warning: the melt fraction is greater than 1: " << inner_melt_fraction << " at T = " << temperature << " deg C, and P = " << pressure << " Pa. "
          //                   << "In this condition, the solidus is " << solidus << " deg C, and the liquidus is " << liquidus << " deg C. "
          //                   << "Setting the melt fraction to 1.\n";
          inner_melt_fraction = 1.0;
        }

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
      // // For we have bulk concentration fields advected,
      // // we don't need to extract the porosity from the old solution
      // std::vector<double> old_porosity(in.n_evaluation_points());
      // // we want to get the porosity field from the old solution here,
      // // because we need a field that is not updated in the nonlinear iterations
      // if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid
      //     && this->get_timestep_number() > 0 && !this->get_parameters().use_operator_splitting)
      //   {
      //     // Prepare the field function

      //     Functions::FEFieldFunction<dim, LinearAlgebra::BlockVector>

      //     fe_value(this->get_dof_handler(), this->get_old_solution(), this->get_mapping());

      //     const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

      //     fe_value.set_active_cell(in.current_cell);
      //     fe_value.value_list(in.position,
      //                         old_porosity,
      //                         this->introspection().component_indices.compositional_fields[porosity_idx]);
      //   }

      // we need to match compositional fields for components
      std::vector<unsigned int> compositional_field_indices = {};
      for (unsigned int component_idx = 0; component_idx < n_components; ++component_idx)
        {
          compositional_field_indices.push_back(component_indices_list[component_idx].bulk_index);
        }
      // check if the length of compositional_field_indices equals to n_components
      AssertThrow(compositional_field_indices.size() == n_components,
                  ExcMessage("The number of compositional fields for components do not match the expected number of components."));

      const FluidPressureInputs<dim> *fluid_pressure_input =
        in.template get_additional_input<FluidPressureInputs<dim>>();
      
      
      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          const double pressure_for_material
            = (fluid_pressure_input != nullptr
               && fluid_pressure_input->fluid_pressure.size() == in.n_evaluation_points()
               && std::isfinite(fluid_pressure_input->fluid_pressure[q]))
              ? fluid_pressure_input->fluid_pressure[q]
              : in.pressure[q];

          // if (this->get_parameters().use_operator_splitting)
          //   {
          //     const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
          //     old_porosity[q] = in.composition[q][porosity_idx];
          //   }
 
          if (this->include_melt_transport())
            {
              // // we temporarily don't distinguish between melt fraction (mass) and porosity (volume)
              // const double old_melt_fraction = old_porosity[q];
              std::vector<double> bulk_concentrations(n_components);
              for (unsigned int _i = 0; _i < n_components; ++_i)
                {
                  bulk_concentrations[_i] = in.composition[q][compositional_field_indices[_i]];
                }
              const double temperature_for_equilibrium_calculation = in.temperature[q] - ZERO_CELSIUS_IN_KELVIN;
              melt_fractions[q] = this->solve_eq_melt_fraction(temperature_for_equilibrium_calculation,
                                                               std::max(0.0, pressure_for_material),
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
      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
      std::vector<double> old_porosity(in.n_evaluation_points(), 0.0);
      std::vector<std::vector<double>> old_fields = in.composition; 
      // initialize old_fields with current_linearization_point as a fallback, 
      // in case we don't get the old state input or the size doesn't match
      // however we should not use current_linearization_point as old fields
      // because linearal guessing and picard iteration
      // will update the current_linearization_point and make it different from the old solution

      const OldFieldStateInputs<dim> *old_state_input =
        in.template get_additional_input<OldFieldStateInputs<dim>>();
      if (old_state_input != nullptr
          && old_state_input->old_porosity.size() == in.n_evaluation_points()
          && old_state_input->old_fields.size() == in.n_evaluation_points())
        {
          old_porosity = old_state_input->old_porosity;
          old_fields = old_state_input->old_fields;
        }
      else
        {
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            // we should not use current_linearization_point as old fields
            old_porosity[i] = in.composition[i][porosity_idx];
        }

      const CompositionFieldGradientsInputs<dim> *composition_gradient_input =
        in.template get_additional_input<CompositionFieldGradientsInputs<dim>>();

      const FluidPressureInputs<dim> *fluid_pressure_input =
        in.template get_additional_input<FluidPressureInputs<dim>>();

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();
      PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>();
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();
      ComponentPhaseExchangeOutputs<dim> *phase_exchange_out = out.template get_additional_output<ComponentPhaseExchangeOutputs<dim>>();
      const bool use_pressure_temperature_viscosity =
        (viscosity_activation_energy > 0.0 || viscosity_activation_volume > 0.0);

      // const Quadrature<dim> &quadrature_formula_compositional_fields = this->introspection().quadratures.compositional_fields;
      // const unsigned int n_q_points_compositional_fields = quadrature_formula_compositional_fields.size();
      // FEValues<dim> fe_values_for_grad_C (this->get_mapping(),
      //                                     this->get_fe(),
      //                                     quadrature_formula_compositional_fields,
      //                                     update_values | update_gradients | update_quadrature_points);
      // // The following "grad_C_values" is actually the gradients of all compositional fields
      // // the index of this vector is the index for compositional fields.
      // std::vector<std::vector<Tensor<1,dim>>> grad_C_values; 
      // // to avoid quadrature point mismatch, we define cell average gradients
      // std::vector<Tensor<1,dim>> grad_C(num_compositional_fields);
      // // we need to judge if the cell is valid
      // if (this->include_melt_transport() && in.current_cell.state() == IteratorState::valid)
      //   {     
      //     fe_values_for_grad_C.reinit (in.current_cell);
      //     for (unsigned int c=0; c<num_compositional_fields; ++c)
      //       {
      //         grad_C_values.push_back(std::vector<Tensor<1,dim>>(n_q_points_compositional_fields));
      //         fe_values_for_grad_C[this->introspection().extractors.compositional_fields[c]].get_function_gradients (
      //           this->get_current_linearization_point(), grad_C_values[c]);
      //       }
      //     // to avoid quadrature point mismatch, we calculate and use cell average gradients
      //     for (unsigned int c=0; c<num_compositional_fields; ++c)
      //       {
      //         for (unsigned int q=0; q<n_q_points_compositional_fields; ++q)
      //           grad_C[c] += grad_C_values[c][q];
      //         grad_C[c] /= n_q_points_compositional_fields;
      //       }
      //   }

      std::vector<std::vector<Tensor<1,dim>>> grad_C_values;
      grad_C_values.assign(in.n_evaluation_points(), std::vector<Tensor<1,dim>>(num_compositional_fields));
       if (composition_gradient_input != nullptr
          && composition_gradient_input->composition_gradients.size() == in.n_evaluation_points())
        {
          for (unsigned int c=0; c<num_compositional_fields; ++c)
            {
              for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
                {
                  grad_C_values[q][c] = composition_gradient_input->composition_gradients[q][c];
                }
            }
        }

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          const double pressure_for_material
            = (fluid_pressure_input != nullptr
               && fluid_pressure_input->fluid_pressure.size() == in.n_evaluation_points()
               && std::isfinite(fluid_pressure_input->fluid_pressure[i]))
              ? fluid_pressure_input->fluid_pressure[i]
              : in.pressure[i];

          // calculate density first, we need it for the reaction term
          // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
          double temperature_dependence = 1.0;
          if (this->include_adiabatic_heating ())
            temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                      * thermal_expansivity;
          else
            temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;

          double pressure_dependence = 1.0 + compressibility * (pressure_for_material - this->get_surface_pressure());
          // std::exp(melt_compressibility * (pressure_for_material - this->get_surface_pressure()));

          // // keep density constant for now
          // out.densities[i] = reference_rho_s;

          // now we enable the dependencies of density on T and P
          // TODO: maybe we can also include the dependence of density on composition in the future
          out.densities[i] = reference_rho_s * temperature_dependence * pressure_dependence;

          out.viscosities[i] = eta_0;
          // By default, no melting or freezing --> set all reactions to zero
          for (unsigned int c=0; c<in.composition[i].size(); ++c)
            {
              out.reaction_terms[i][c] = 0.0;

              if (this->get_parameters().use_operator_splitting && reaction_rate_out != nullptr)
                reaction_rate_out->reaction_rates[i][c] = 0.0;
            }

          if (phase_exchange_out != nullptr)
            {
              phase_exchange_out->effective_latent_heat[i] = 0.0;
              phase_exchange_out->partial_phi_partial_T[i] = 0.0;
              phase_exchange_out->melting_rate[i] = 0.0;
            }

          if (this->include_melt_transport())
            {
              const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
              const double porosity = std::min(1.0, std::max(in.composition[i][porosity_idx],0.0));

              // Keep the original porosity-weakening relation for shear viscosity.
              out.viscosities[i] *= std::exp(- alpha_phi * porosity);

              // initialize prescribed field outputs if needed
              if (prescribed_field_out != nullptr)
                {
                  for (unsigned int c=0; c<num_compositional_fields; ++c)
                    {
                      prescribed_field_out->prescribed_field_outputs[i][c] = in.composition[i][c];
                    }
                }

              if (include_melting_and_freezing && (in.requests_property(MaterialProperties::reaction_terms) ||
                                                   in.requests_property(MaterialProperties::reaction_rates) ||
                                                   in.requests_property(MaterialProperties::viscosity)))
                {
                  // Keep both divergence retrieval pathways for future A/B tests:
                  // 1) original trace(in.strain_rate), 2) additional input velocity_divergence.
                  double trace_strain_rate = numbers::signaling_nan<double>();
                  double velocity_divergence = numbers::signaling_nan<double>();
                  if (this->get_timestep_number() > 0)
                    {
                      Assert(std::isfinite(in.strain_rate[i].norm()),
                             ExcMessage("Invalid strain_rate in the MaterialModelInputs. This is likely because it was "
                                        "not filled by the caller."));
                      trace_strain_rate = trace(in.strain_rate[i]);

                      const bool has_velocity_divergence =
                        (composition_gradient_input != nullptr
                         && composition_gradient_input->velocity_divergence.size() == in.n_evaluation_points()
                         && std::isfinite(composition_gradient_input->velocity_divergence[i]));

                      if (has_velocity_divergence)
                        // bad news: interpolate method don't fill the additional input for the velocity divergence
                        velocity_divergence = composition_gradient_input->velocity_divergence[i];
                      else
                        velocity_divergence = trace_strain_rate;
                    }
                  // (void) trace_strain_rate;
                  // (void) velocity_divergence;

                  // we temporarily don't distinguish between melt fraction (mass) and porosity (volume)
                  const double old_melt_fraction = old_porosity[i];
                  std::vector<double> bulk_concentrations(n_components);
                  if (enable_equilibrium_calculation)
                    {
                      for (unsigned int _i = 0; _i < n_components; ++_i)
                        {
                          bulk_concentrations[_i] = in.composition[i][component_indices_list[_i].bulk_index];
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
                                                                                std::max(0.0, pressure_for_material),
                                                                                bulk_concentrations) : 
                                                   old_melt_fraction);
                  double porosity_change = 0.0;
                  if (prescribed_field_out != nullptr)
                    {
                      prescribed_field_out->prescribed_field_outputs[i][porosity_idx] = eq_melt_fraction;

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
                  porosity_change = eq_melt_fraction - old_porosity[i];
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
                          melting_points[comp_idx] = temperature_melting(std::max(0.0, pressure_for_material),
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
                  
                  // before getting into the compositional fields loop
                  // we need to resize the custom variables in MeltOutputs
                  if (melt_out != nullptr && enable_equilibrium_calculation)
                    {
                      // melt_out->concentrations_in_phases 
                      // is a standard vector with a length corresponding to the number of quadrature point
                      // containing standard vectors with a length corresponding to the number of compositional fields
                      // (Yes, only the composition field corresponding to the chemical components will be filled.)
                      // containing pairs of concentration values in solid and fluid phases
                      // first: solid phase, second: fluid phase
                      // we've resized the outer vector according to the number of quadrature points
                      // now we need to resize the inner vectors according to the number of compositional fields
                      melt_out->concentrations_in_phases[i].resize(this->introspection().n_compositional_fields);

                      // melt_out->concentration_gradients_in_phases
                      // is a standard vector with a length corresponding to the number of quadrature point
                      // containing standard vectors with a length corresponding to the number of compositional fields
                      // (Yes, only the composition field corresponding to the chemical components will be filled.)
                      // containing pairs of dealii::Tensor<1,dim> for concentration gradient values in solid and fluid phases
                      // first: solid phase, second: fluid phase
                      // we've resized the outer vector according to the number of quadrature points
                      // now we need to resize the inner vectors according to the number of compositional fields
                      melt_out->concentration_gradients_in_phases[i].resize(this->introspection().n_compositional_fields);
                    }
                  
                  // to take variables out of the loop
                  double temporary_c_l_an_eq = 0.0;
                  Tensor<1,dim> temporary_grad_c_l_an_eq;
                  
                  // now we can calculate the equilibrium concentrations of each component
                  for (unsigned int c = 0; c < in.composition[i].size(); ++c)
                    {
                      if (c == porosity_idx && this->get_timestep_number() > 0)
                        // this term will be treated as melting rate, it should multiply density and divide timestep
                        // we decide to calculate reaction term of porosity field
                        // whether equilibrium calculation is enabled or disabled
                        // for we've already deal with the disabled case above for porosity

                        // // we don't need reaction term for porosity under extended Boussinesq approximation
                        // out.reaction_terms[i][c] = porosity_change * out.densities[i] / this->get_timestep();
                        out.reaction_terms[i][c] = 0.0;

                      else if (true /* this->get_timestep_number() > 0 */ && enable_equilibrium_calculation) 
                      {
                        // reaction term for bulk concentration fields
                        // if equilibrium calculation is disabled
                        // the method "match_component_index_to_composition_index" would not work
                        // and compositional_field_indices would be meaningless

                        // 1st we need c_s_i_eq and c_l_i_eq
                        double c_s_i_eq = 0.0;
                        double c_l_i_eq = 0.0;
                        for (unsigned int component_idx = 0; component_idx < n_components; ++component_idx)
                          {
                            if (c == component_indices_list[component_idx].bulk_index)
                              {
                                c_s_i_eq = c_solid_eq_values[component_idx];
                                c_l_i_eq = c_liquid_eq_values[component_idx];
                                break;
                              }
                          }
                        temporary_c_l_an_eq = c_liquid_eq_values[1]; 

                        // 2nd we need gradient of c_s_i_eq and c_l_i_eq
                        // the new method is to calculate the gradients of c_s_i_eq and c_l_i_eq as prescribed compositional fields in the input of material model evaluation, and then we can directly extract them without numerical derivatives
                        Tensor<1,dim> grad_c_s_i_eq;
                        Tensor<1,dim> grad_c_l_i_eq;                        
                        for (unsigned int component_idx = 0; component_idx < n_components; ++component_idx)
                          {
                            if (c == component_indices_list[component_idx].bulk_index)
                              {
                                if (component_indices_list[component_idx].solid_index != numbers::invalid_unsigned_int)
                                  grad_c_s_i_eq = grad_C_values[i][component_indices_list[component_idx].solid_index];
                                if (component_indices_list[component_idx].liquid_index != numbers::invalid_unsigned_int)
                                  grad_c_l_i_eq = grad_C_values[i][component_indices_list[component_idx].liquid_index];
                                break;
                              }
                          }
                        temporary_grad_c_l_an_eq = grad_C_values[i][component_indices_list[1].liquid_index]; // anorthite component index is 1

                        // we don't fill reaction_terms of bulk concentration fields
                        out.reaction_terms[i][c] = 0.0;
                        // instead we fill our costume output variable in MeltOutputs object.
                        if (melt_out != nullptr)
                          {
                            // // we've already resized the outer and inner vectors above
                            melt_out->concentrations_in_phases[i][c] = std::make_pair(c_s_i_eq, c_l_i_eq);
                            melt_out->concentration_gradients_in_phases[i][c] = std::make_pair(grad_c_s_i_eq, grad_c_l_i_eq);
                          }
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
                          if (reaction_rate_out != nullptr)
                            {
                              // placeholder
                            }
                          out.reaction_terms[i][c] = 0.0;
                        }
                    } // end of looping over compositional fields

                  // To calculate the latent heat term, we need to calculate and output
                  // the effective_latent_heat and partial_phi_partial_T.
                  // We need to do a numerical derivation to get them.

                  // To do this, we calculate an equilibrium state at a perturbed temperature,
                  // which is slightly higher than the current temperature.
                  const double temperature_perturbation = 1e-3;
                  const double perturbed_temperature_for_equilibrium_calculation
                    = temperature_for_equilibrium_calculation + temperature_perturbation;

                  double eq_melt_fraction_perturbed = eq_melt_fraction;
                  std::vector<double> eq_consts_perturbed(n_components);
                  std::vector<double> c_solid_eq_values_perturbed(n_components);
                  std::vector<double> c_liquid_eq_values_perturbed(n_components);

                  if (enable_equilibrium_calculation)
                    {
                      eq_melt_fraction_perturbed
                        = this->solve_eq_melt_fraction(perturbed_temperature_for_equilibrium_calculation,
                                                       std::max(0.0, pressure_for_material),
                                                       bulk_concentrations);

                      for (unsigned int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                        {
                          // The melting point depends on pressure only, so we can reuse melting_points.
                          eq_consts_perturbed[comp_idx]
                            = equilibrium_constant(perturbed_temperature_for_equilibrium_calculation,
                                                   latent_heat_values[comp_idx],
                                                   tuning_parameter_values[comp_idx],
                                                   melting_points[comp_idx]);
                          c_solid_eq_values_perturbed[comp_idx]
                            = calculate_concentration_solid(bulk_concentrations[comp_idx],
                                                            eq_melt_fraction_perturbed,
                                                            eq_consts_perturbed[comp_idx]);
                          c_liquid_eq_values_perturbed[comp_idx]
                            = calculate_concentration_liquid(bulk_concentrations[comp_idx],
                                                             eq_melt_fraction_perturbed,
                                                             eq_consts_perturbed[comp_idx]);
                        }
                    }

                  // // Keep latent-heat related outputs as placeholders for now.
                  // // The perturbed-state variables above are prepared for future use.
                  // (void) eq_melt_fraction_perturbed;
                  // (void) eq_consts_perturbed;
                  // (void) c_solid_eq_values_perturbed;
                  // (void) c_liquid_eq_values_perturbed;

                  std::vector<double> perturbational_heat_effects(n_components, 0.0);
                  double sum_perturbational_heat_effect = 0.0;
                  double perturbational_porosity_change = 0.0;
                  double latent_heat_eff = 0.0;
                  double partial_phi_partial_T_local = 0.0;
                  if (enable_equilibrium_calculation)
                    {
                      // calculate L_eff
                      for (unsigned int component_idx = 0; component_idx < n_components; ++component_idx)
                        {
                          perturbational_heat_effects[component_idx] = latent_heat_values[component_idx] 
                            * ( eq_melt_fraction_perturbed * c_liquid_eq_values_perturbed[component_idx] 
                              - eq_melt_fraction * c_liquid_eq_values[component_idx]);
                          sum_perturbational_heat_effect += perturbational_heat_effects[component_idx];
                        }
                      perturbational_porosity_change = eq_melt_fraction_perturbed - eq_melt_fraction;
                      if (std::abs(perturbational_porosity_change) < 1e-12)
                        {
                          latent_heat_eff = 0.0;
                        }
                      else
                        {
                          latent_heat_eff = sum_perturbational_heat_effect / perturbational_porosity_change;
                          // latent_heat_eff = 450000.0; // temporary constant
                        }  
                      // calculate partial_phi_partial_T
                      partial_phi_partial_T_local = perturbational_porosity_change / temperature_perturbation;
                    }
                  
                  // here we calculate the melting rate
                  double melting_rate = 0.0;
                  // for debug out, copy "partial_phi_partial_time" to a prescribed field output
                  double temp_partial_phi_partial_time = 0.0;
                  if (enable_equilibrium_calculation && this->get_timestep_number() > 0)
                    {
                      // order-1 euler backward method for now.
                      const double partial_phi_partial_time = (eq_melt_fraction - old_melt_fraction) / this->get_timestep();
                      const double phi_advection = in.velocity[i] * grad_C_values[i][porosity_idx];
                      const double phi_compaction = - (1 - eq_melt_fraction) * velocity_divergence; // this is a placeholder, we can also use trace_strain_rate here
                      melting_rate = partial_phi_partial_time + phi_advection + phi_compaction;
                      temp_partial_phi_partial_time = partial_phi_partial_time;
                      // this->get_pcout() 
                      // << "[MM eval]"
                      // << "  partial_phi_partial_time: " 
                      // << partial_phi_partial_time 
                      // << "  phi_advection: " 
                      // << phi_advection 
                      // << "  phi_compaction: " 
                      // << phi_compaction 
                      // // << "  velocity_divergence: "
                      // // << velocity_divergence
                      // // << "  trace_strain_rate: "
                      // // << trace_strain_rate
                      // << "  melting_rate: "
                      // << melting_rate
                      // << "  temp_partial_phi_partial_time: "
                      // << temp_partial_phi_partial_time
                      // << std::endl;
                    }
                  
                  // fill phase exchange outputs if the model includes melt transport and the outputs are requested
                  if (phase_exchange_out != nullptr && enable_equilibrium_calculation)
                    {
                      phase_exchange_out->effective_latent_heat[i] = latent_heat_eff;
                      phase_exchange_out->partial_phi_partial_T[i] = partial_phi_partial_T_local;
                      phase_exchange_out->melting_rate[i] = melting_rate;
                    }

                  if (prescribed_field_out != nullptr)
                    {
                      // temporary output c_bulk_an for debugging
                      // TODO: write a paragraph to find and match debug_field_indices
                      if (fill_debug_fields)
                        {
                          std::vector<unsigned int> debug_field_indices;
                          debug_field_indices.push_back(this->introspection().compositional_index_for_name("debug_1"));
                          debug_field_indices.push_back(this->introspection().compositional_index_for_name("debug_2"));
                          // prescribed_field_out->prescribed_field_outputs[i][porosity_idx] = eq_melt_fraction; // placeholder
                          unsigned int checking_component_index = 1;
                          double T_m_0 = melting_point_0_values[checking_component_index];
                          double m_c_A = melting_curve_coefficient_A_values[checking_component_index];
                          double m_c_B = melting_curve_coefficient_B_values[checking_component_index];
                          const double temp_T_m_morb = temperature_melting(pressure_for_material, T_m_0, m_c_A, m_c_B);
                          double L_teq = latent_heat_values[checking_component_index];
                          double r_teq = tuning_parameter_values[checking_component_index];
                          const double temp_K_morb = equilibrium_constant(temperature_for_equilibrium_calculation, L_teq, r_teq, temp_T_m_morb);
                          prescribed_field_out->prescribed_field_outputs[i][debug_field_indices[0]] = temp_K_morb;
                          checking_component_index = 2;
                          T_m_0 = melting_point_0_values[checking_component_index];
                          m_c_A = melting_curve_coefficient_A_values[checking_component_index];
                          m_c_B = melting_curve_coefficient_B_values[checking_component_index];
                          const double temp_T_m_cmorb = temperature_melting(pressure_for_material, T_m_0, m_c_A, m_c_B);
                          L_teq = latent_heat_values[checking_component_index];
                          r_teq = tuning_parameter_values[checking_component_index];
                          const double temp_K_cmorb = equilibrium_constant(temperature_for_equilibrium_calculation, L_teq, r_teq, temp_T_m_cmorb);
                          prescribed_field_out->prescribed_field_outputs[i][debug_field_indices[1]] = temp_K_cmorb; 
                        }
                      // fill phase concentration fields
                      for (unsigned int component_idx = 0; component_idx < n_components; ++component_idx)
                        {
                          unsigned int solid_field_index = component_indices_list[component_idx].solid_index;
                          unsigned int liquid_field_index = component_indices_list[component_idx].liquid_index;
                          if (solid_field_index != numbers::invalid_unsigned_int)
                            prescribed_field_out->prescribed_field_outputs[i][solid_field_index] = c_solid_eq_values[component_idx];
                          if (liquid_field_index != numbers::invalid_unsigned_int)
                            prescribed_field_out->prescribed_field_outputs[i][liquid_field_index] = c_liquid_eq_values[component_idx];
                        }
                      // fill named melting_rate field
                      if (fill_prescribed_melting_rate_field)
                        {
                          if (this->introspection().compositional_name_exists("melting_rate"))
                            {
                              unsigned int melting_rate_field_index = this->introspection().compositional_index_for_name("melting_rate");
                              prescribed_field_out->prescribed_field_outputs[i][melting_rate_field_index] = melting_rate;
                            }
                        }
                    }
                  
                }  
                
            } // end of "if include_melt_transport()"

          out.entropy_derivative_pressure[i]    = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = 0.0;

          if (use_pressure_temperature_viscosity)
            out.viscosities[i] *= pressure_temperature_viscosity_factor(in.temperature[i], pressure_for_material);
          else
            {
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
        } // end of loop over evaluation points

      // // fill melt outputs if they exist
      // MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();
      // To use melt outputs earlier, we move this part into the loop above

      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              const double pressure_for_material
                = (fluid_pressure_input != nullptr
                   && fluid_pressure_input->fluid_pressure.size() == in.n_evaluation_points()
                   && std::isfinite(fluid_pressure_input->fluid_pressure[i]))
                  ? fluid_pressure_input->fluid_pressure[i]
                  : in.pressure[i];

              double porosity = std::max(in.composition[i][porosity_idx],0.0);
              melt_out->fluid_viscosities[i] = eta_f; // TODO: add compositional dependence
              // here we set the permeability to be proportional to phi^3 instead of phi^3*(1-phi)^2 temporarily.
              // because we want to keep it the same to the solitary wave benchmark.
              melt_out->permeabilities[i] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<0>(1.0-porosity);
              melt_out->fluid_density_gradients[i] = Tensor<1,dim>();

              // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
              double temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                          * thermal_expansivity;
              else
                temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
              melt_out->fluid_densities[i] = reference_rho_f * temperature_dependence
                                             * std::exp(melt_compressibility * (pressure_for_material - this->get_surface_pressure()));

              melt_out->compaction_viscosities[i] = xi_0 * std::exp(- alpha_phi * porosity);

              if (use_pressure_temperature_viscosity)
                melt_out->compaction_viscosities[i] *= pressure_temperature_viscosity_factor(in.temperature[i], pressure_for_material);
              else
                {
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

              // fill melt outputs for concentration fields if equilibrium calculation is enabled
              if (enable_equilibrium_calculation)
                {
                  // placeholder
                }

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
          prm.declare_entry ("Viscosity activation energy", "0.0",
                             Patterns::Double (0.),
                             "Activation energy $E$ for the pressure-temperature dependent "
                             "Arrhenius viscosity factor. If either this value or the activation "
                             "volume is nonzero, the same factor is applied to both shear and "
                             "compaction viscosities. Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Viscosity activation volume", "0.0",
                             Patterns::Double (0.),
                             "Activation volume $V$ for the pressure-temperature dependent "
                             "Arrhenius viscosity factor. The reference state is defined by the "
                             "surface pressure and reference temperature. Units: "
                             "\\si{\\meter\\cubed\\per\\mole}.");
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
          prm.declare_entry ("Enable chemical reaction rate", "true",
                             Patterns::Bool (),
                             "Choosing if the material model gives out reaction_rate vector"
                             "for advection equations. If false, the reaction_rate vector would be zeros."
                             "declare only, no parse.");
          prm.declare_entry ("Fill debug fields", "true",
                             Patterns::Bool (),
                             "Whether to fill debug fields for melt model. "
                             "If true, some debug fields will be filled in the prescribed_field_outputs. "
                             "These fields can be used for debugging purposes.");
          prm.declare_entry ("Fill prescribed melting rate field", "true",
                             Patterns::Bool (),
                             "Whether to fill the prescribed melting_rate field in the prescribed_field_outputs. "
                             "If true, the melting_rate field will be filled in the prescribed_field_outputs. "
                             "This field can be used for debugging purposes or for coupling with other models. "
                             "If false, the melting_rate field will not be filled and will be set to zero.");

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
                             "Temperature in degree Celsius! "
                             "T_m = T_m0 + A * P + B * P^2. "
                            //  "Units: \\si{\\kelvin}" // description
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
                             "A small background porosity added to avoid zero permeability."
                             "declare only, no parse.");    

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
          viscosity_activation_energy       = prm.get_double ("Viscosity activation energy");
          viscosity_activation_volume       = prm.get_double ("Viscosity activation volume");
          thermal_expansivity               = prm.get_double ("Thermal expansion coefficient");
          alpha_phi                         = prm.get_double ("Exponential melt weakening factor");
          compressibility                   = prm.get_double ("Solid compressibility");
          melt_compressibility              = prm.get_double ("Melt compressibility");

          thermal_conductivity              = prm.get_double ("Thermal conductivity");
          include_melting_and_freezing      = prm.get_bool ("Include melting and freezing");
          melting_time_scale                = prm.get_double ("Melting time scale for operator splitting");

          enable_equilibrium_calculation    = prm.get_bool ("Enable equilibrium calculation");
          // enable_chemical_reaction_rate     = prm.get_bool ("Enable chemical reaction rate");
          fill_debug_fields                 = prm.get_bool ("Fill debug fields");
          fill_prescribed_melting_rate_field = prm.get_bool ("Fill prescribed melting rate field");

          equilibrium_solving_method = prm.get ("Equilibrium solving method");

          AssertThrow(this->get_parameters().use_operator_splitting == false,
                      ExcMessage("Error: Material model Melt thermodynamic equilibrium "
                                 "can not support operator splitting."));

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
          
          // background_porosity = prm.get_double("Background porosity");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt thermodynamic equilibrium with Thermal viscosity exponent can not have reference_T=0."));

          if ((viscosity_activation_energy > 0.0 || viscosity_activation_volume > 0.0)
              && reference_T <= 0.0)
            AssertThrow(false,
                        ExcMessage("Error: Material model Melt thermodynamic equilibrium with pressure-temperature "
                                   "dependent viscosity requires a positive reference temperature."));

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
                  const bool match_result = match_component_index_to_composition_index();
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
    MeltThermodynamicEquilibrium<dim>::
    fill_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &input,
                                          const LinearAlgebra::BlockVector        &solution,
                                          const FEValuesBase<dim>                 &fe_values,
                                          const Introspection<dim>                &introspection) const
    {
      if (input.template get_additional_input<OldFieldStateInputs<dim>>() == nullptr)
        input.additional_inputs.push_back(
          std::make_unique<OldFieldStateInputs<dim>>(input.n_evaluation_points(),
                                                     introspection.n_compositional_fields));

      if (input.template get_additional_input<CompositionFieldGradientsInputs<dim>>() == nullptr)
        input.additional_inputs.push_back(
          std::make_unique<CompositionFieldGradientsInputs<dim>>(introspection.n_compositional_fields));

      if (this->include_melt_transport()
          && input.template get_additional_input<FluidPressureInputs<dim>>() == nullptr)
        input.additional_inputs.push_back(
          std::make_unique<FluidPressureInputs<dim>>(input.n_evaluation_points()));

      for (unsigned int i=0; i<input.additional_inputs.size(); ++i)
        {
          OldFieldStateInputs<dim> *old_state_input =
            dynamic_cast<OldFieldStateInputs<dim> *>(input.additional_inputs[i].get());

          if (old_state_input != nullptr)
            old_state_input->fill(this->get_old_solution(), fe_values, introspection);
          else
            input.additional_inputs[i]->fill(solution, fe_values, introspection);
        }
    }


    template <int dim>
    void
    MeltThermodynamicEquilibrium<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<MaterialModel::PrescribedFieldOutputs<dim>>() == nullptr)
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

  } // end of namespace MaterialModel
} // end of namespace aspect                  

# endif // end of the new implementation


# ifndef ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS
// The original implementation, advecting concentration fields in phases.
#  pragma message("ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS is OFF, using the original implementation advecting concentration fields in phases.")

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

    template <int dim>
    double
    MeltThermodynamicEquilibrium<dim>::
    pressure_temperature_viscosity_factor (const double temperature,
                                           const double pressure) const
    {
      const double universal_gas_constant = 8.31446261815324;
      const double clamped_temperature = std::max(temperature, 1.0);
      const double clamped_reference_temperature = std::max(reference_T, 1.0);
      const double local_pressure = std::max(0.0, pressure);
      const double reference_pressure = std::max(0.0, this->get_surface_pressure());

      const double log_viscosity_reference =
        (viscosity_activation_energy + viscosity_activation_volume * reference_pressure)
        / (universal_gas_constant * clamped_reference_temperature);
      const double log_viscosity_local =
        (viscosity_activation_energy + viscosity_activation_volume * local_pressure)
        / (universal_gas_constant * clamped_temperature);

      const double bounded_log_factor =
        std::max(std::min(log_viscosity_local - log_viscosity_reference, 50.0), -50.0);

      return std::exp(bounded_log_factor);
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

      const FluidPressureInputs<dim> *fluid_pressure_input =
        in.template get_additional_input<FluidPressureInputs<dim>>();
      
      for (unsigned int q=0; q<in.n_evaluation_points(); ++q)
        {
          const double pressure_for_material
            = (fluid_pressure_input != nullptr
               && fluid_pressure_input->fluid_pressure.size() == in.n_evaluation_points()
               && std::isfinite(fluid_pressure_input->fluid_pressure[q]))
              ? fluid_pressure_input->fluid_pressure[q]
              : in.pressure[q];

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
                                                               std::max(0.0, pressure_for_material),
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
      const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");
      std::vector<double> old_porosity(in.n_evaluation_points(), 0.0);
      std::vector<std::vector<double>> old_fields = in.composition;

      const OldFieldStateInputs<dim> *old_state_input =
        in.template get_additional_input<OldFieldStateInputs<dim>>();
      if (old_state_input != nullptr
          && old_state_input->old_porosity.size() == in.n_evaluation_points()
          && old_state_input->old_fields.size() == in.n_evaluation_points())
        {
          old_porosity = old_state_input->old_porosity;
          old_fields = old_state_input->old_fields;
        }
      else
        {
          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            old_porosity[i] = in.composition[i][porosity_idx];
        }

      ReactionRateOutputs<dim> *reaction_rate_out = out.template get_additional_output<ReactionRateOutputs<dim>>();
      PrescribedFieldOutputs<dim> *prescribed_field_out = out.template get_additional_output<PrescribedFieldOutputs<dim>>();
      const FluidPressureInputs<dim> *fluid_pressure_input =
        in.template get_additional_input<FluidPressureInputs<dim>>();
      const bool use_pressure_temperature_viscosity =
        (viscosity_activation_energy > 0.0 || viscosity_activation_volume > 0.0);

      for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
        {
          const double pressure_for_material
            = (fluid_pressure_input != nullptr
               && fluid_pressure_input->fluid_pressure.size() == in.n_evaluation_points()
               && std::isfinite(fluid_pressure_input->fluid_pressure[i]))
              ? fluid_pressure_input->fluid_pressure[i]
              : in.pressure[i];

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

              // Keep the original porosity-weakening relation for shear viscosity.
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
                                                                                std::max(0.0, pressure_for_material),
                                                                                bulk_concentrations) : 
                                                   old_melt_fraction);
                  double porosity_change = 0.0;
                  if (prescribed_field_out != nullptr)
                    {
                      prescribed_field_out->prescribed_field_outputs[i][0] = std::min(eq_melt_fraction, 1.0);

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
                  porosity_change = std::min(eq_melt_fraction, 1.0) - old_porosity[i];
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
                          melting_points[comp_idx] = temperature_melting(std::max(0.0, pressure_for_material),
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
                        if (enable_chemical_reaction_rate) 
                          out.reaction_terms[i][c] = (c_eq - old_fields[i][c]);
                        else
                          out.reaction_terms[i][c] = 0.0;
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
                      // TODO: write a paragraph to find and match debug_field_indices
                      if (fill_debug_fields)
                        {
                          prescribed_field_out->prescribed_field_outputs[i][5] = (c_liquid_eq_values[0] - old_fields[i][liquid_indices[0]]);
                          prescribed_field_out->prescribed_field_outputs[i][6] = (c_liquid_eq_values[1] - old_fields[i][liquid_indices[1]]);
                          prescribed_field_out->prescribed_field_outputs[i][7] = (eq_melt_fraction - old_porosity[i]);
                        }
                    }

                }  
            }

          out.entropy_derivative_pressure[i]    = 0.0;
          out.entropy_derivative_temperature[i] = 0.0;
          out.thermal_expansion_coefficients[i] = thermal_expansivity;
          out.specific_heat[i] = reference_specific_heat;
          out.thermal_conductivities[i] = thermal_conductivity;
          out.compressibilities[i] = 0.0;

          if (use_pressure_temperature_viscosity)
            out.viscosities[i] *= pressure_temperature_viscosity_factor(in.temperature[i], pressure_for_material);
          else
            {
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
        }    

      // fill melt outputs if they exist
      MeltOutputs<dim> *melt_out = out.template get_additional_output<MeltOutputs<dim>>();

      if (melt_out != nullptr)
        {
          const unsigned int porosity_idx = this->introspection().compositional_index_for_name("porosity");

          for (unsigned int i=0; i<in.n_evaluation_points(); ++i)
            {
              const double pressure_for_material
                = (fluid_pressure_input != nullptr
                   && fluid_pressure_input->fluid_pressure.size() == in.n_evaluation_points()
                   && std::isfinite(fluid_pressure_input->fluid_pressure[i]))
                  ? fluid_pressure_input->fluid_pressure[i]
                  : in.pressure[i];

              double porosity = std::max(in.composition[i][porosity_idx],0.0);

              // // seems unnecessary to fill melt_out for porosity here
              // // for we've fill it as the prescribed field.
              // melt_out->porosities[i] = 0.0;

              melt_out->fluid_viscosities[i] = eta_f;
              // here we set the permeability to be proportional to phi^3 instead of phi^3*(1-phi)^2 temporarily.
              // because we want to keep it the same to the solitary wave benchmark.
              melt_out->permeabilities[i] = reference_permeability * Utilities::fixed_power<3>(porosity) * Utilities::fixed_power<0>(1.0-porosity);
              melt_out->fluid_density_gradients[i] = Tensor<1,dim>();

              // temperature dependence of density is 1 - alpha * (T - T(adiabatic))
              double temperature_dependence = 1.0;
              if (this->include_adiabatic_heating ())
                temperature_dependence -= (in.temperature[i] - this->get_adiabatic_conditions().temperature(in.position[i]))
                                          * thermal_expansivity;
              else
                temperature_dependence -= (in.temperature[i] - reference_T) * thermal_expansivity;
              melt_out->fluid_densities[i] = reference_rho_f * temperature_dependence
                                             * std::exp(melt_compressibility * (pressure_for_material - this->get_surface_pressure()));

              melt_out->compaction_viscosities[i] = xi_0 * std::exp(- alpha_phi * porosity);

              if (use_pressure_temperature_viscosity)
                melt_out->compaction_viscosities[i] *= pressure_temperature_viscosity_factor(in.temperature[i], pressure_for_material);
              else
                {
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
          prm.declare_entry ("Viscosity activation energy", "0.0",
                             Patterns::Double (0.),
                             "Activation energy $E$ for the pressure-temperature dependent "
                             "Arrhenius viscosity factor. If either this value or the activation "
                             "volume is nonzero, the same factor is applied to both shear and "
                             "compaction viscosities. Units: \\si{\\joule\\per\\mole}.");
          prm.declare_entry ("Viscosity activation volume", "0.0",
                             Patterns::Double (0.),
                             "Activation volume $V$ for the pressure-temperature dependent "
                             "Arrhenius viscosity factor. The reference state is defined by the "
                             "surface pressure and reference temperature. Units: "
                             "\\si{\\meter\\cubed\\per\\mole}.");
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
          prm.declare_entry ("Enable chemical reaction rate", "true",
                             Patterns::Bool (),
                             "Choosing if the material model gives out reaction_rate vector"
                             "for advection equations. If false, the reaction_rate vector would be zeros.");
          prm.declare_entry ("Fill debug fields", "true",
                             Patterns::Bool (),
                             "Whether to fill debug fields for melt model. "
                             "If true, some debug fields will be filled in the prescribed_field_outputs. "
                             "These fields can be used for debugging purposes.");

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
                             "A small background porosity added to avoid zero permeability."
                             "declare only, no parse.");    

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
          viscosity_activation_energy       = prm.get_double ("Viscosity activation energy");
          viscosity_activation_volume       = prm.get_double ("Viscosity activation volume");
          thermal_expansivity               = prm.get_double ("Thermal expansion coefficient");
          alpha_phi                         = prm.get_double ("Exponential melt weakening factor");
          compressibility                   = prm.get_double ("Solid compressibility");
          melt_compressibility              = prm.get_double ("Melt compressibility");

          thermal_conductivity              = prm.get_double ("Thermal conductivity");
          include_melting_and_freezing      = prm.get_bool ("Include melting and freezing");
          melting_time_scale                = prm.get_double ("Melting time scale for operator splitting");

          enable_equilibrium_calculation    = prm.get_bool ("Enable equilibrium calculation");
          enable_chemical_reaction_rate     = prm.get_bool ("Enable chemical reaction rate");
          fill_debug_fields                 = prm.get_bool ("Fill debug fields");

          equilibrium_solving_method = prm.get ("Equilibrium solving method");

          AssertThrow(this->get_parameters().use_operator_splitting == false,
                      ExcMessage("Error: Material model Melt thermodynamic equilibrium "
                                 "can not support operator splitting."));

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
          
          // background_porosity = prm.get_double("Background porosity");

          if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
            AssertThrow(false, ExcMessage("Error: Material model Melt thermodynamic equilibrium with Thermal viscosity exponent can not have reference_T=0."));

          if ((viscosity_activation_energy > 0.0 || viscosity_activation_volume > 0.0)
              && reference_T <= 0.0)
            AssertThrow(false,
                        ExcMessage("Error: Material model Melt thermodynamic equilibrium with pressure-temperature "
                                   "dependent viscosity requires a positive reference temperature."));

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
    MeltThermodynamicEquilibrium<dim>::
    fill_additional_material_model_inputs(MaterialModel::MaterialModelInputs<dim> &input,
                                          const LinearAlgebra::BlockVector        &solution,
                                          const FEValuesBase<dim>                 &fe_values,
                                          const Introspection<dim>                &introspection) const
    {
      if (input.template get_additional_input<OldFieldStateInputs<dim>>() == nullptr)
        input.additional_inputs.push_back(
          std::make_unique<OldFieldStateInputs<dim>>(input.n_evaluation_points(),
                                                     introspection.n_compositional_fields));

      if (this->include_melt_transport()
          && input.template get_additional_input<FluidPressureInputs<dim>>() == nullptr)
        input.additional_inputs.push_back(
          std::make_unique<FluidPressureInputs<dim>>(input.n_evaluation_points()));

      for (unsigned int i=0; i<input.additional_inputs.size(); ++i)
        {
          OldFieldStateInputs<dim> *old_state_input =
            dynamic_cast<OldFieldStateInputs<dim> *>(input.additional_inputs[i].get());

          if (old_state_input != nullptr)
            old_state_input->fill(this->get_old_solution(), fe_values, introspection);
          else
            input.additional_inputs[i]->fill(solution, fe_values, introspection);
        }
    }


    template <int dim>
    void
    MeltThermodynamicEquilibrium<dim>::create_additional_named_outputs (MaterialModel::MaterialModelOutputs<dim> &out) const
    {
      if (out.template get_additional_output<MaterialModel::PrescribedFieldOutputs<dim>>() == nullptr)
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

# endif // end of the original implementation

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
