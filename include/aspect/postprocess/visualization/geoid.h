/*
  Copyright (C) 2016 - 2021 by the authors of the ASPECT code.

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


#ifndef _aspect_postprocess_visualization_geoid_h
#define _aspect_postprocess_visualization_geoid_h

#include <aspect/simulator_access.h>
#include <aspect/postprocess/visualization.h>
#include <aspect/postprocess/geoid.h>


namespace aspect
{
  namespace Postprocess
  {
    namespace VisualizationPostprocessors
    {
      /**
       * A class derived from DataPostprocessorScalar that computes a variable
       * that represents the geoid topography. This quantity, strictly
       * speaking, only makes sense at the surface of the domain. Thus, the
       * value is set to zero in all the cells inside of the domain.
       *
       * The member functions are all implementations of those declared in the
       * base class. See there for their meaning.
       */
      template <int dim>
      class Geoid
        : public DataPostprocessorScalar<dim>,
          public SimulatorAccess<dim>,
          public Interface<dim>
      {
        public:
          Geoid();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program after parse_parameters is run and after
           * the SimulatorAccess (if applicable) is initialized.
           */
          void
          initialize () override;

          /**
           * @copydoc dealii::DataPostprocessor<dim>::evaluate_vector_field()
           */
          void
          evaluate_vector_field(const DataPostprocessorInputs::Vector<dim> &input_data,
                                std::vector<Vector<double>> &computed_quantities) const override;


          /**
           * Let the postprocessor manager know about the other postprocessors
           * this one depends on. Specifically, the Geoid postprocessor.
           */
          std::list<std::string>
          required_other_postprocessors() const override;

        private:

      };
    }
  }
}

#endif
