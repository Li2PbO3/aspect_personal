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
