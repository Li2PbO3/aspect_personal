(sec:extending:plugin-types:geometry-models)=
# Geometry models

The geometry model is responsible for describing the domain in which we want
to solve the equations. A domain is described in
deal.II by a coarse mesh and, if necessary, an object
that characterizes the boundary. Together, these two suffice to reconstruct
any domain by adaptively refining the coarse mesh and placing new nodes
generated by refining cells onto the surface described by the boundary object.
The geometry model is also responsible for marking different parts of the
boundary with different *boundary indicators* for which one can then, in the
input file, select whether these boundaries should be Dirichlet-type (fixed
temperature) or Neumann-type (no heat flux) boundaries for the temperature,
and what kind of velocity conditions should hold there. In
deal.II, a boundary indicator is a number of (integer) type
`types::boundary_id`, but since boundaries are hard to remember and get right
in input files, geometry models also have a function that provide a map from
symbolic names that can be used to describe pieces of the boundary to the
corresponding boundary indicators. For example, the simple `box` geometry
model in 2d provides the map
`{"left" -> 0, "right" -> 1, "bottom" -> 2,"top" -> 3}`,
and we have consistently used these symbolic names in the input files used in
this manual.

To implement a new geometry model, you need to overload the
[aspect::GeometryModel::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1GeometryModel.html)
class and use the
`ASPECT_REGISTER_GEOMETRY_MODEL` macro to register your new class. The
implementation of the new class should be in namespace
`aspect::GeometryModel`.

The functions you need to overload are extensively
discussed in the documentation of this interface class at
[aspect::GeometryModel::Interface](https://aspect.geodynamics.org/doc/doxygen/namespaceaspect_1_1GeometryModel.html).
Note that the `create_coarse_mesh()` function you need to provide does not only create the actual mesh (i.e.,
the locations of the vertices of the coarse mesh and how they connect to
cells) but it must also set the boundary indicators for all parts of the
boundary of the mesh. The deal.II glossary
describes the purpose of boundary indicators as follows:

> In a `Triangulation` object, every part of the boundary is associated with a
> unique number (of type `types::boundary_id`) that is used to identify which
> boundary geometry object is responsible to generate new points when the mesh
> is refined. By convention, this boundary indicator is also often used to
> determine what kinds of boundary conditions are to be applied to a
> particular part of a boundary. The boundary is composed of the faces of the
> cells and, in 3d, the edges of these faces.
>
> By default, all boundary indicators of a mesh are zero, unless you are
> reading from a mesh file that specifically sets them to something different,
> or unless you use one of the mesh generation functions in namespace
> `GridGenerator` that have a 'colorize' option. A typical piece
> of code that sets the boundary indicator on part of the boundary to
> something else would look like this, here setting the boundary indicator to
> 42 for all faces located at $x=-1$:
>
> ```{code-block} c++
> for (auto &cell : triangulation.active_cell_iterators())
>     for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
>       if (cell->face(f)->at_boundary())
>         if (cell->face(f)->center(true)[0] == -1)
>           cell->face(f)->set_boundary_indicator (42);
> ```
>
> This calls functions `TriaAccessor::set_boundary_indicator`. In 3d, it may
> also be appropriate to call `TriaAccessor::set_all_boundary_indicators`
> instead on each of the selected faces. To query the boundary indicator of a
> particular face or edge, use `TriaAccessor::boundary_indicator`.
>
> The code above only sets the boundary indicators of a particular part of the
> boundary, but it does not by itself change the way the Triangulation class
> treats this boundary for the purposes of mesh refinement. For this, you need
> to call `Triangulation::set_boundary` to associate a boundary object with a
> particular boundary indicator. This allows the Triangulation object to use a
> different method of finding new points on faces and edges to be refined; the
> default is to use a `StraightBoundary` object for all faces and edges. The
> results section of step-49 has a worked example that shows all of this in
> action.
>
> The second use of boundary indicators is to describe not only which geometry
> object to use on a particular boundary but to select a part of the boundary
> for particular boundary conditions. *\[...\]*
>
> **Note:** Boundary indicators are inherited from mother faces and edges to
> their children upon mesh refinement. Some more information about boundary
> indicators is also presented in a section of the documentation of the
> Triangulation class.

Two comments are in order here. First, if a coarse triangulation's faces
already accurately represent where you want to pose which boundary condition
(for example to set temperature values or determine which are no-flow and
which are tangential flow boundary conditions), then it is sufficient to set
these boundary indicators only once at the beginning of the program since they
will be inherited upon mesh refinement to the child faces. Here, *at the
beginning of the program* is equivalent to inside the `create_coarse_mesh()`
function of the geometry module shown above that generates the coarse mesh.

Secondly, however, if you can only accurately determine which boundary
indicator should hold where on a refined mesh &ndash; for example because the
coarse mesh is the cube $[0,L]^3$ and you want to have a fixed velocity
boundary describing an extending slab only for those faces for which
$z>L-L_{\text{slab}}$ &ndash; then you need a way to set the boundary
indicator for all boundary faces either to the value representing the slab or
the fluid underneath *after every mesh refinement step*. By doing so, child
faces can obtain boundary indicators different from that of their parents.
deal.II triangulations support this kind of
operations using a so-called *post-refinement signal*. In essence, what this
means is that you can provide a function that will be called by the
triangulation immediately after every mesh refinement step.

The way to do this is by writing a function that sets boundary indicators and
that will be called by the `Triangulation` class. The triangulation does not
provide a pointer to itself to the function being called, nor any other
information, so the trick is to get this information into the function. C++
provides a nice mechanism for this that is best explained using an example:

```{code-block} c++
#include <functional>

    template <int dim>
    void set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation)
    {
      ... set boundary indicators on the triangulation object ...
    }

    template <int dim>
    void
    MyGeometry<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      ... create the coarse mesh ...

      coarse_grid.signals.post_refinement.connect
        (std::bind (&set_boundary_indicators<dim>,
                     std::ref(coarse_grid)));

    }
```

What the call to `std::bind` does is to produce an object that can be
called like a function with no arguments. It does so by taking the address of
a function that does, in fact, take an argument but permanently fix this one
argument to a reference to the coarse grid triangulation. After each
refinement step, the triangulation will then call the object so created which
will in turn call `set_boundary_indicators<dim>` with the reference to the
coarse grid as argument.

This approach can be generalized. In the example above, we have used a global
function that will be called. However, sometimes it is necessary that this
function is in fact a member function of the class that generates the mesh,
for example because it needs to access run-time parameters. This can be
achieved as follows: assuming the `set_boundary_indicators()` function has
been declared as a (non-static, but possibly private) member function of the
`MyGeometry` class, then the following will work:

```{code-block} c++
    #include <functional>

    template <int dim>
    void
    MyGeometry<dim>::
    set_boundary_indicators (parallel::distributed::Triangulation<dim> &triangulation) const
    {
      ... set boundary indicators on the triangulation object ...
    }

    template <int dim>
    void
    MyGeometry<dim>::
    create_coarse_mesh (parallel::distributed::Triangulation<dim> &coarse_grid) const
    {
      ... create the coarse mesh ...

      coarse_grid.signals.post_refinement.connect
        (std::bind (&MyGeometry<dim>::set_boundary_indicators,
                     std::cref(*this),
                     std::ref(coarse_grid)));
    }
```

Here, like any other member function, `set_boundary_indicators` implicitly
takes a pointer or reference to the object it belongs to as first argument.
`std::bind` again creates an object that can be called like a global function
with no arguments, and this object in turn calls `set_boundary_indicators`
with a pointer to the current object and a reference to the triangulation to
work on. Note that because the `create_coarse_mesh` function is declared as
`const`, it is necessary that the `set_boundary_indicators` function is also
declared `const`.

:::{note}
For reasons that have to do with the way the `parallel::distributed::Triangulation` class is
implemented, functions that have been attached to the post-refinement signal of the triangulation
are called more than once, sometimes several times, every time the triangulation is actually refined.
:::
