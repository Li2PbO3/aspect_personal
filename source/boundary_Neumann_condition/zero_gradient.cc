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

#include <aspect/boundary_Neumann_condition/zero_gradient.h>

#pragma message("Compiling boundary_Neumann_condition/zero_gradient.cc")

namespace aspect
{
    namespace BoundaryNeumannCondition
    {
        template <int dim>
        Tensor<1,dim>
        ZeroGradient<dim>::
        boundary_velocity_gradient (const types::boundary_id,
                                    const Point<dim> &,
                                    const Tensor<1,dim> &) const
        {
            // return a zero tensor regardless of position
            return Tensor<1,dim>();
        }
    }
}

// explicit instantiations
namespace aspect
{
    namespace BoundaryNeumannCondition
    {
        ASPECT_REGISTER_BOUNDARY_NEUMANN_CONDITION_MODEL(ZeroGradient,
                                                        "zero gradient",
                                                        "Implementation of a model in which the gradient of "
                                                        "the velocity at the boundary is zero. ")
    }
}