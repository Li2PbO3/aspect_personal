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

#include <aspect/boundary_Neumann_condition/ascii_data.h>
#include <aspect/global.h>

#pragma message("Compiling boundary_Neumann_condition/ascii_data.cc")

namespace aspect
{
    namespace BoundaryNeumannCondition
    {
        template <int dim>
        AsciiData<dim>::AsciiData () = default;

        template <int dim>
        void AsciiData<dim>::initialize ()
        {
            for (const auto &bv : this->get_boundary_Neumann_condition())
                if (bv.second.get() == this)
                    boundary_ids.insert(bv.first);

            AssertThrow(*(boundary_ids.begin()) != numbers::invalid_boundary_id,
                        ExcMessage("Did not find the boundary indicator for the Neumann condition ascii data plugin."));

            Utilities::AsciiDataBoundary<dim>::initialize(boundary_ids, 1);
        }

        template <int dim>
        Tensor<1,dim>
        AsciiData<dim>::boundary_velocity_gradient (const types::boundary_id boundary_indicator,
                                                    const Point<dim> &position,
                                                    const Tensor<1,dim> &normal_vector) const
        {
            const double gradient_magnitude = Utilities::AsciiDataBoundary<dim>::get_data_component(boundary_indicator,
                                                                                          position,
                                                                                          0);
            return -gradient_magnitude * normal_vector;
        }

        template <int dim>
        void AsciiData<dim>::update()
        {
            Interface<dim>::update();
            Utilities::AsciiDataBoundary<dim>::update();
        }

        template <int dim>
        void AsciiData<dim>::declare_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Boundary Neumann condition model");
            {
                Utilities::AsciiDataBoundary<dim>::declare_parameters(prm,
                                                                      "$ASPECT_SOURCE_DIR/data/boundary-Neumann-condition/ascii-data/test/",
                                                                    "box_2d_%s.%d.txt");
            }
            prm.leave_subsection();
        }

        template <int dim>
        void AsciiData<dim>::parse_parameters(ParameterHandler &prm)
        {
            prm.enter_subsection("Boundary Neumann condition model");
            {
                Utilities::AsciiDataBoundary<dim>::parse_parameters(prm);
            }
            prm.leave_subsection();
        }
    }
}

// explicit instantiations
namespace aspect
{
    namespace BoundaryNeumannCondition
    {
        ASPECT_REGISTER_BOUNDARY_NEUMANN_CONDITION_MODEL(AsciiData,
                                                        "ascii data",
                                                        "Implementation of a model in which the gradient of the velocity at the boundary is given by values read from ascii data files. "
                                                        "The data files should be formatted as described in the documentation of Utilities::AsciiDataBoundary. "
                                                        "The data files should give the magnitude of the velocity gradient normal to the boundary, "
                                                        "and the direction of the velocity gradient is assumed to be normal to the boundary and pointing inward.");
    }
}