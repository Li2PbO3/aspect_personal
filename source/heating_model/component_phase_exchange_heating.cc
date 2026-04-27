/*
  Copyright (C) 2015 - 2026 by the authors of the ASPECT code.

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


#include <aspect/heating_model/component_phase_exchange_heating.h>

#if defined(ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS)

#include <aspect/material_model/melt_thermodynamic_equilibrium.h>
#include <aspect/plugins.h>

#include <cmath>

# pragma message("Compiling component_phase_exchange_heating.cc with ASPECT_MELT_ADVECTING_BULK_CONCENTRATIONS ON")

namespace aspect
{
  namespace HeatingModel
  {
    template <int dim>
    void
    ComponentPhaseExchangeHeating<dim>::
    initialize ()
    {
      AssertThrow(this->include_melt_transport(),
                  ExcMessage("Heating model component phase exchange heating requires melt transport to be enabled."));

      const auto &material_model = this->get_material_model();
      AssertThrow(Plugins::plugin_type_matches<const MaterialModel::MeltThermodynamicEquilibrium<dim>>(material_model),
                  ExcMessage("Heating model component phase exchange heating requires material model "
                             "'melt thermodynamic equilibrium' with bulk concentration advection enabled."));
    }


    template <int dim>
    void
    ComponentPhaseExchangeHeating<dim>::
    evaluate (const MaterialModel::MaterialModelInputs<dim> &material_model_inputs,
              const MaterialModel::MaterialModelOutputs<dim> &material_model_outputs,
              HeatingModel::HeatingModelOutputs &heating_model_outputs) const
    {
      // // Parameter wiring is in place; detailed source-form implementation will be added later.
      // (void) use_heat_source_formulation;

      Assert(heating_model_outputs.heating_source_terms.size() == material_model_inputs.position.size(),
             ExcMessage("Heating outputs need to have the same number of entries as the material model inputs."));

      const MaterialModel::ComponentPhaseExchangeOutputs<dim> *phase_exchange_out
        = material_model_outputs.template get_additional_output<MaterialModel::ComponentPhaseExchangeOutputs<dim>>();

      AssertThrow(phase_exchange_out != nullptr,
                  ExcMessage("Need ComponentPhaseExchangeOutputs from the material model for component phase exchange heating."));

      AssertThrow(phase_exchange_out->effective_latent_heat.size() == material_model_inputs.position.size(),
                  ExcMessage("The size of effective_latent_heat must match the number of evaluation points."));
      AssertThrow(phase_exchange_out->partial_phi_partial_T.size() == material_model_inputs.position.size(),
                  ExcMessage("The size of partial_phi_partial_T must match the number of evaluation points."));

      const unsigned int melting_rate_field_index =
        this->introspection().compositional_name_exists("melting_rate")
        ? this->introspection().compositional_index_for_name("melting_rate")
        : numbers::invalid_unsigned_int;

      if (use_heat_source_formulation)
        {
          AssertThrow(melting_rate_field_index != numbers::invalid_unsigned_int,
                      ExcMessage("Compositional field 'melting_rate' must exist when using heat-source "
                                 "formulation in component phase exchange heating."));
          AssertThrow(material_model_inputs.composition.size() == material_model_inputs.position.size(),
                      ExcMessage("The size of composition inputs must match the number of evaluation points."));
        }

      for (unsigned int q = 0; q < heating_model_outputs.heating_source_terms.size(); ++q)
        {
          heating_model_outputs.heating_source_terms[q] = 0.0;
          heating_model_outputs.lhs_latent_heat_terms[q] = 0.0;
          heating_model_outputs.rates_of_temperature_change[q] = 0.0;

          const double effective_latent_heat = phase_exchange_out->effective_latent_heat[q];
          const double partial_phi_partial_T = phase_exchange_out->partial_phi_partial_T[q];

          AssertThrow(std::isfinite(effective_latent_heat),
                      ExcMessage("effective_latent_heat contains a non-finite value."));
          AssertThrow(std::isfinite(partial_phi_partial_T),
                      ExcMessage("partial_phi_partial_T contains a non-finite value."));
          if (use_heat_source_formulation)
            {
              // this->get_pcout() << "[HeatingModel] " 
              //                   << "Using heat source formulation for latent heat. "
              //                   << std::endl;
              AssertIndexRange(melting_rate_field_index,
                               material_model_inputs.composition[q].size());

              // Read melting rate from current_linearization_point composition.
              const double melting_rate =
                material_model_inputs.composition[q][melting_rate_field_index];
              AssertThrow(std::isfinite(melting_rate),
                          ExcMessage("melting_rate contains a non-finite value."));
              const double solid_density = material_model_outputs.densities[q];
              (void) solid_density;
              const MaterialModel::MeltOutputs<dim> *melt_out = material_model_outputs.template get_additional_output<MaterialModel::MeltOutputs<dim>>();
              const double liquid_density = melt_out->fluid_densities[q];

              /**
               * The latent heat term is applied as a source term on the right hand side of the temperature equation, 
               * so the unit of the latent heat term is W/m^3, which is the same as the unit of heating source terms.
               * 
               * melting absorbs heat, so the latent heat term is negative when melting_rate is positive.
               * The effective_latent_heat is positive when melting, and negative when freezing, 
               * so the latent heat term is negative when melting_rate is positive, and positive when freezing, 
               * which is consistent with the physical expectation that melting absorbs heat and freezing releases heat.
               */
              // here the melting_rate is the changing rate of melt fraction or porosity.
              // this means the melting_rate is within the unit of 1/s or 1/yr
              heating_model_outputs.heating_source_terms[q] -= effective_latent_heat * liquid_density * melting_rate;
              // this->get_pcout()
              // << "[HeatingModel] "
              // << "heating_source_terms is " << heating_model_outputs.heating_source_terms[q]
              // << " at position " << material_model_inputs.position[q]
              // << std::endl;
            }
          else
            {
              /**
               * The latent heat term is applied to the left hand side of the temperature equation, 
               * added with the "density times specific heat" term.
               * So the unit of the latent heat term is J/m^3/K, which is the same as the unit of density times specific heat.
               * 
               * The densities here is of the solid phase, for this latent heat is being with the solid term.
               */
              heating_model_outputs.lhs_latent_heat_terms[q]
                = material_model_outputs.densities[q] * effective_latent_heat * partial_phi_partial_T;
            }
        }
    }


    template <int dim>
    void
    ComponentPhaseExchangeHeating<dim>::
    create_additional_material_model_outputs(MaterialModel::MaterialModelOutputs<dim> &outputs) const
    {
      if (outputs.template get_additional_output<MaterialModel::ComponentPhaseExchangeOutputs<dim>>() != nullptr)
        return;

      const unsigned int n_points = outputs.n_evaluation_points();
      outputs.additional_outputs.push_back(
        std::make_unique<MaterialModel::ComponentPhaseExchangeOutputs<dim>>(n_points));
    }


    template <int dim>
    void
    ComponentPhaseExchangeHeating<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Component phase exchange heating");
        {
          prm.declare_entry ("Use heat source formulation", "true",
                             Patterns::Bool(),
                             "Select how latent heat is represented in the temperature equation. "
                             "If true, latent heat is intended to be added as a source term on the RHS. "
                             "If false, latent heat is intended to be added to the apparent heat capacity "
                             "on the LHS. "
                             "Default is true (heat-source form).");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    ComponentPhaseExchangeHeating<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Heating model");
      {
        prm.enter_subsection("Component phase exchange heating");
        {
          use_heat_source_formulation = prm.get_bool("Use heat source formulation");
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
  namespace HeatingModel
  {
    ASPECT_REGISTER_HEATING_MODEL(ComponentPhaseExchangeHeating,
                                  "component phase exchange heating",
                                  "Heating model that adds latent heat in an apparent-heat-capacity form "
                                  "to the left hand side of the temperature equation.")
  }
}

#endif
