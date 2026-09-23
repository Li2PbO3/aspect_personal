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


#ifndef _aspect_postprocess_melting_rate_instantaneous_h
#define _aspect_postprocess_melting_rate_instantaneous_h

#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>

#include <deal.II/base/table_handler.h>


namespace aspect
{
  namespace Postprocess
  {
    /**
     * A postprocessor that computes the melting rate $\Gamma/\rho_l$ and the
     * three component-wise values $\Gamma^i/\rho_l$ by perturbing the
     * instantaneous thermodynamic equilibrium with a small convective
     * perturbation.
     *
     * This is the "external" single-step estimator of the carbon_melting
     * postprocessing pipeline, but evaluated on the in-memory (double
     * precision) solution and with the material model's own equilibrium solver.
     * It therefore avoids the float32 output round trip of the external
     * workflow and does not require solving an additional advection equation.
     *
     * The liquid mass balance reads
     * @f[
     *   \frac{\Gamma^i}{\rho_l}
     *   = \frac{\partial (\phi c_l^i)}{\partial t}
     *     + \nabla\cdot(\phi c_l^i \vec v_l).
     * @f]
     * The first term is approximated by a finite difference of the equilibrium
     * map $\Phi$ along the rates of change of the state variables,
     * @f[
     *   \frac{\partial(\phi c_l^i)}{\partial t}
     *   \approx \frac{\Phi^i(x+\tau\dot x)-\Phi^i(x)}{\tau},
     * @f]
     * where $x = (P, T, \bar c)$ and $\dot x$ is the advective rate of change of
     * the bulk composition (from the model's own bulk mass flux) and the
     * finite-difference rates of temperature and fluid pressure.
     *
     * The result is written to `<output directory>/melting_rate_instantaneous/`
     * as one ASCII file per time step, in double precision. Physical units: 1/s.
     */
    template <int dim>
    class MeltingRateInstantaneous : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Destructor.
         */
        ~MeltingRateInstantaneous () override = default;

        /**
         * Execute this postprocessor.
         */
        std::pair<std::string,std::string>
        execute (TableHandler &statistics) override;

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
         * The virtual time step tau [yr] used to build the perturbation.
         */
        double tau_yr;

        /**
         * Whether the bulk composition is perturbed.
         */
        bool perturb_composition;

        /**
         * Whether the temperature is perturbed.
         */
        bool perturb_temperature;

        /**
         * Whether the fluid pressure is perturbed.
         */
        bool perturb_pressure;

        /**
         * Whether a central (instead of a forward) difference is used.
         */
        bool use_central_difference;

        /**
         * Whether the perturbation direction (rate of change of the bulk
         * composition) is evaluated in the same cell-averaged discrete form
         * that the advection system solves (with the arithmetic cell average
         * of div(V) and the c_l*div(V) - phi*(u_f-V).grad(c_l) source),
         * instead of the pointwise flux divergence.
         */
        bool use_cell_averaged_form;

        /**
         * Whether the melting rate is set to zero where there is no melt.
         */
        bool zero_where_no_melt;

        /**
         * The melt fraction below which the melting rate is set to zero.
         */
        double melt_fraction_threshold;

        /**
         * End members of the tie line of the "uniform phase" calibration case,
         * used for the tie-line-consistency diagnostic written to the output
         * file. For a state on the tie line the three recovered melt fractions
         * f_i = (c_bar_i - S_i)/(L_i - S_i) agree; their spread measures the
         * violation of the tie line.
         */
        std::vector<double> tie_line_solid;

        /**
         * @copydoc tie_line_solid
         */
        std::vector<double> tie_line_liquid;
    };

  }
}

#endif
