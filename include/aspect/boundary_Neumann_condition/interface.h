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

#ifndef _aspect_boundary_Neumann_condition_interface_h
#define _aspect_boundary_Neumann_condition_interface_h

#include <aspect/plugins.h>
#include <aspect/geometry_model/interface.h>

#include <deal.II/base/point.h>
#include <deal.II/base/parameter_handler.h>

#pragma message("Compiling boundary_Neumann_condition/interface.h")

namespace aspect
{
    /**
     * A namespace in which we define everything that has to do with defining
     * Neumann boundary conditions for the Stokes equations.
     * 
     * The namespace "BoundaryNeumannCondition" and its contents are defined for
     * implementing Neumann boundary conditions for velocity.
     * 
     * @ingroup BoundaryNeumannCondition
     */
    namespace BoundaryNeumannCondition
    {
        using namespace dealii;

        /** 
         * A base class for parameterizations of Neumann boundary conditions.
         * 
         * @ingroup BoundaryNeumannCondition
         */
        template <int dim>
        class Interface
        {
            public:
                /** 
                 * Destructor. Made virtual to enforce that derived classes also have
                 * virtual destructors.
                 */
                virtual ~Interface() = default;

                /**
                 * Initialization function. This function is called once at the
                 * beginning of the program after parse_parameters is run and after the
                 * SimulatorAccess (if applicable) is initialized.
                 */
                virtual void initialize();

                /**
                 * A function that is called at the beginning of each time step.
                 * The default implementation of the function does nothing, but
                 * derived classes can overload it to update the boundary condition
                 * values as a function of time.
                 */
                virtual void update();

                /**
                 * Return the value of the boundary velocity gradient at a particular position
                 * on the boundary of the domain.
                 * 
                 * @param boundary_indicator The boundary indicator of the boundary part 
                 * on which we ask for the boundary Neumann condition.
                 * @param position The position of the point at which we ask for the
                 * boundary Neumann condition.
                 * @param normal_vector The (outward) normal vector to the boundary
                 * of the domain.
                 * 
                 * @return The value of the boundary velocity gradient at position @p position.
                 */
                virtual Tensor<1,dim> boundary_velocity_gradient(const types::boundary_id boundary_indicator,
                                                                 const Point<dim> &position,
                                                                 const Tensor<1,dim> &normal_vector) const = 0;
                
                /** 
                 * Declare the parameters this class takes through input files.
                 * The default implementation of this function does not describe any
                 * parameters. Consequently, derived classes do not have to overload
                 * this function if they do not take any runtime parameters. 
                 */
                static void declare_parameters(ParameterHandler &prm);

                /** 
                 * Read the parameters this class declares from the parameter file.
                 * The default implementation of this function does not read any
                 * parameters. Consequently, derived classes do not have to overload
                 * this function if they do not take any runtime parameters. 
                 */
                virtual void parse_parameters(ParameterHandler &prm);

            protected:
                /** 
                 * Pointer to the geometry object in use.
                 */
                const GeometryModel::Interface<dim> *geometry_model;
        }; // end of class Interface




        /** 
         * Register a Neumann boundary conditions model so that it can be
         * selected from the parameter file.
         * 
         * @param name A string that identifies the Neumann boundary conditions model. 
         * This string will be listed in the documentation of the parameter file.
         * @param description A text description of what this model does and that
         * will be listed in the documentation of the parameter file.
         * @param declare_parameters_function A pointer to a function that can be
         * used to declare the parameters that this Neumann boundary conditions
         * model wants to read from input files.
         * @param factory_function A pointer to a function that can be used to
         * create an object that describes this Neumann boundary conditions model.
         */
        template <int dim>
        void register_boundary_Neumann_condition(const std::string &name,
                                                const std::string &description,
                                                void (*declare_parameters_function)(ParameterHandler &),
                                                std::unique_ptr<Interface<dim>> (*factory_function)());
                                        
        /**
         * A function that given the name of a model returns a pointer to an
         * object that describes it. Ownership of the pointer is transferred to
         * the caller.
         * 
         * @param name The name of the Neumann boundary conditions model as given in the parameter file.
         * @return A pointer to an object that describes the Neumann boundary conditions model with the given name.
         * 
         * The model object returned is not yet initialized and has not read its runtime parameters yet.
         * 
         * @ingroup BoundaryNeumannConditions
         */
        template <int dim>
        std::unique_ptr<Interface<dim>> create_boundary_Neumann_condition(const std::string &name);

        /**
         * Return a list of names of all implemented Neumann boundary conditions models,
         * separated by '|' so that it can be used in an object of type Patterns::Selection.
         */
        template <int dim>
        std::string get_names();

        /**
         * Declare the parameters this class takes through input files.
         * The default implementation of this function does not describe any parameters. 
         * Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        template <int dim>
        void declare_parameters(ParameterHandler &prm);

        /** 
         * For the current plugin, write a graph of all of the plugins we know about,
         * in the format that the programs dot and neato understand.
         * This allows for a visualization of how all of the plugins that ASPECT knows about are
         * interconnected, and connect to other part of the ASPECT code.
         * 
         * @param output_stream The stream to which the graph should be written.
         */
        template <int dim>
        void write_plugin_graph(std::ostream &output_stream);



        /**
         * Given a class name, a name, and a description for the parameter file
         * for a Neumann boundary conditions model, register it with the
         * functions that can declare their parameters and create these objects.
         * 
         * @ingroup BoundaryNeumannCondition
         */
#define ASPECT_REGISTER_BOUNDARY_NEUMANN_CONDITION_MODEL(classname, name, description) \
        template class classname<2>; \
        template class classname<3>; \
        namespace ASPECT_REGISTER_BOUNDARY_NEUMANN_CONDITION_MODEL_ ## classname \
        { \
          aspect::internal::Plugins::RegisterHelper<aspect::BoundaryNeumannCondition::Interface<2>, classname<2>> \
          dummy_ ## classname ## _2d (&aspect::BoundaryNeumannCondition::register_boundary_Neumann_condition<2>, name, description); \
          aspect::internal::Plugins::RegisterHelper<aspect::BoundaryNeumannCondition::Interface<3>, classname<3>> \
          dummy_ ## classname ## _3d (&aspect::BoundaryNeumannCondition::register_boundary_Neumann_condition<3>, name, description); \
        }
    } // end of namespace BoundaryNeumannCondition
} // end of namespace aspect

# endif // end of include guard