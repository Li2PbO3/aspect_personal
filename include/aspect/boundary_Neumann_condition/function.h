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

// Personaly attached head file.

#ifndef _aspect_boundary_Neumann_condition_function_h
#define _aspect_boundary_Neumann_condition_function_h

#include <aspect/boundary_Neumann_condition/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

#include <deal.II/base/parsed_function.h>

#pragma message("Compiling boundary_Neumann_condition/function.h")

namespace aspect
{
    namespace BoundaryNeumannCondition
    {
        using namespace dealii;

        /** 
         * A class that implements Neumann boundary conditions based on a
         * functional description provided in the input file.
         * 
         * @ingroup BoundaryNeumannCondition
         */
        template <int dim>
        class Function : public Interface<dim>, public SimulatorAccess<dim>
        {
            public:
                /**
                 * Constructor.
                 */
                Function ();

                /**
                 * Return the boundary velocity gradient as a function of position. 
                 * The (outward) normal vector to the domain is also provided as a second argument.
                 */
                Tensor<1,dim>
                boundary_velocity_gradient (const types::boundary_id boundary_indicator,
                                            const Point<dim> &position,
                                            const Tensor<1,dim> &normal_vector) const override;

                /**
                 * A function that is called at the beginning of each time step to
                 * indicate what the model time is for which the boundary values will
                 * next be evaluated. For the current class, the function passes to 
                 * the parsed function what the current time is.
                 */
                void update() override;

                /**
                 * Declare the parameters this class takes through input files.
                 * The default implementation of this function does not describe any parameters. 
                 * Consequently, derived classes do not have to overload
                 * this function if they do not take any runtime parameters.
                 */
                static void declare_parameters(ParameterHandler &prm);

                /**
                 * Read the parameters this class declares from the parameter file.
                 * The default implementation of this function does not read any parameters.
                 * Consequently, derived classes do not have to overload
                 * this function if they do not take any runtime parameters.
                 */
                void parse_parameters(ParameterHandler &prm) override;
            
            private:
                /** 
                 * A function object that describes the components of the gradient of 
                 * the velocity at the boundary as a function of position. 
                 */
                Functions::ParsedFunction<dim> boundary_velocity_gradient_function;

                /**
                 * The coordinate representation to evaluate the function.
                 * Possible choices are depth cartesian and spherical coordinates. 
                 * The default is cartesian coordinates.
                 */
                Utilities::Coordinates::CoordinateSystem coordinate_system;

                /**
                 * Whether to specify the gradient in x, y, z components, or r, phi, theta components.
                 * The default is to specify the gradient in x, y, z components.
                 */
                bool use_spherical_unit_vectors;
        };
    }
}
#endif