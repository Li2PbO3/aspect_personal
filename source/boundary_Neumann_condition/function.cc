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

#include <aspect/boundary_Neumann_condition/function.h>
#include <aspect/utilities.h>
#include <aspect/global.h>

#pragma message("Compiling boundary_Neumann_condition/function.cc")

namespace aspect
{
    namespace BoundaryNeumannCondition
    {
        template <int dim>
        Function<dim>::Function ()
            :
            boundary_velocity_gradient_function (dim)
        {}


        template <int dim>
        Tensor<1,dim>
        Function<dim>::boundary_velocity_gradient (const types::boundary_id /*boundary_indicator*/,
                                                   const Point<dim> &position,
                                                   const Tensor<1,dim> &/*normal_vector*/) const
        {
            Tensor<1,dim> gradient;
            Utilities::NaturalCoordinate<dim> point =
                this->get_geometry_model().cartesian_to_other_coordinates(position, coordinate_system);
            for (unsigned int d=0; d<dim; ++d)
                gradient[d] = boundary_velocity_gradient_function.value(Utilities::convert_array_to_point<dim>(point.get_coordinates()),d);
            if (use_spherical_unit_vectors)
                gradient = Utilities::Coordinates::spherical_to_cartesian_vector(gradient, position);

            return gradient;
        }


        template <int dim>
        void
        Function<dim>::update()
        {
            // we get time passed as seconds (always) but may want
            // to reinterpret it in years
            if (this->convert_output_to_years())
                boundary_velocity_gradient_function.set_time(this->get_time() / year_in_seconds);
            else
                boundary_velocity_gradient_function.set_time(this->get_time());
        }


        template <int dim>
        void
        Function<dim>::declare_parameters (ParameterHandler &prm)
        {
            prm.enter_subsection("Boundary Neumann condition model");
            {
                prm.enter_subsection("Function");
                {
                    prm.declare_entry("Coordinate system", "cartesian",
                                      Patterns::Selection("depth|cartesian|spherical"),
                                      "The coordinate system in which the function is evaluated. Possible choices are depth, cartesian, spherical, and ellipsoidal coordinates. The default is cartesian coordinates.");
                    prm.declare_entry("Use spherical unit vectors", "false",
                                      Patterns::Bool(),
                                      "Whether to specify the gradient in x, y, z components, or r, phi, theta components. The default is to specify the gradient in x, y, z components.");
                    Functions::ParsedFunction<dim>::declare_parameters (prm, dim);
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();
        }


        template <int dim>
        void
        Function<dim>::parse_parameters (ParameterHandler &prm)
        {
            prm.enter_subsection("Boundary Neumann condition model");
            {
                prm.enter_subsection("Function");
                {
                    coordinate_system = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
                    use_spherical_unit_vectors = prm.get_bool("Use spherical unit vectors");
                    if (use_spherical_unit_vectors)
                        AssertThrow (this->get_geometry_model().natural_coordinate_system() == Utilities::Coordinates::spherical,
                                     ExcMessage ("Spherical unit vectors should not be used "
                                                 "when geometry model is not spherical."));
                    try
                    {
                        boundary_velocity_gradient_function.parse_parameters (prm);
                    }
                    catch (...)
                    {
                        std::cerr << "ERROR: FunctionParser failed to parse\n"
                                  << "\t'Boundary Neumann condition model.Function'\n"
                                  << "with expression\n"
                                  << "\t'" << prm.get("Function expression") << "'"
                                  << "More information about the cause of the parse error \n"
                                  << "is shown below.\n";
                        throw;
                    }
                }
                prm.leave_subsection();
            }
            prm.leave_subsection();
        }
    } // end namespace BoundaryNeumannCondition
} // end namespace aspect

// explicit instantiations
namespace aspect
{
    namespace BoundaryNeumannCondition
    {
        ASPECT_REGISTER_BOUNDARY_NEUMANN_CONDITION_MODEL(Function,
                                                        "function",
                                                        "Implementation of a model in which the gradient of the velocity at the boundary is given by a user-specified function. "
                                                        "The function can be specified as a function of position and time, and the coordinate system in which the function is evaluated can be selected from cartesian, depth, and spherical coordinates. "
                                                        "The function expression should be given in terms of x, y, z or r, phi, theta depending on the coordinate system selected. "
                                                        "For more information about how to specify the function expression, please refer to the documentation of Functions::ParsedFunction in the ASPECT manual.");
    }
}