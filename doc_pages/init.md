title: Initialization

On this page the initialization of the complete Ateles solver is described.
Basicly this process is done with the subroutine
[[atl_initialize_module:atl_initialize]].
During initialization,
the configuration is read from the given Lua configuration file,
all necessary data structures for the computation are build
and the required values are set.

[TOC]

@todo: Describe what is read, how it is read and where it is stored.

## Equation

The equation description gets initialized with the call of
[[atl_equation_init_module:atl_init_equation]].
This subroutine reads the equation system to solve from the configuration file
and calls the equation specific initialization routine
to create the necessary data structures
and to set the required values accordingly.

The currently supported equation systems to load are

-   [Benjamin–Bona–Mahony equation](https://en.wikipedia.org/wiki/Benjamin–Bona–Mahony_equation)
    with [[atl_eqn_bbm_module:atl_load_BBMEM]],
-   [Euler equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics))
    with [[atl_eqn_euler_hlp_module:atl_eqn_euler_init]] or with
    [[atl_eqn_linearEuler_hlp_module:atl_eqn_linearEuler_init]] for linear Euler,
-   [Navier–Stokes equations](https://en.wikipedia.org/wiki/Navier–Stokes_equations)
    with [[atl_eqn_nvrstk_hlp_module:atl_eqn_nvrstk_init]] or with
    [[atl_eqn_filNvrstk_hlp_module:atl_eqn_filtered_nvrstk_init]]
    for filtered Navier-Stokes,
-   [Acoustic wave equation](https://en.wikipedia.org/wiki/Acoustic_wave_equation)
    with [[atl_eqn_acoustic_hlp_module:atl_eqn_acoustic_init]],
-   [Maxwell's equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations)
    with [[atl_eqn_maxwell_hlp_module:atl_eqn_maxwell_init]],
-   [Nernst–Planck equation](https://en.wikipedia.org/wiki/Nernst–Planck_equation)
    with [[atl_eqn_nerplanck_var_module:atl_init_nerplanck_vars]],
-   [Advection](https://en.wikipedia.org/wiki/Advection#The_advection_equation)
    with [[atl_eqn_advection_1d_hlp_module:atl_eqn_advection_1d_init]] and
-   [Heat equation](https://en.wikipedia.org/wiki/Heat_equation)
    with [[atl_eqn_heat_hlp_module:atl_eqn_heat_init]].

## Variable System

The variable system contains all variables used throughout the simulation.
It is constructed during the initialization
based on the data defined by the equation system
and on the user defined variables provided by the configuration file.
This includes also source terms, material properties,
and penalization parameters.

### User Defined Variables

The user can define variables that can serve different purposes.
E.g. they can be used to be the data source for source terms,
material properties or just to track quantities that are not provided
by the used equation system.
User defined variables are specified in the lua configuration file.
They reside in it's own table called `variable`
and follow the notation specified by [[tem_variable_module]].
This table is loaded
and gets added to the variable system during initialization.

### Source Terms

To load the source terms for a given variable system,
we start with creating a list of possible source variables
the variable system can make use of.
Thus these variables depend on the specific equation system,
this is done in the equation specific routines
called by [[atl_equation_init_module:atl_init_equation]].
There the name and the number of components for each variable is added
to the list `initSource%eval_source` of type
[[atl_source_types_module:atl_init_source_type]].
This list is taken to extract the corresponding variables for every source term
from the lua configuration file where they are stored in the table `source`.
Inside this table, the names of user defined variables are assigned
to the names of source terms.
Briefly, the `source` table is a dictionary with the name of the source terms
and the names of the corresponding user defined variables.
With these assignements and the list of possible source terms,
a list of available source terms is created,
which is stored in a list to be added to the variable system.

To get more information about how to set up scource terms
in the configuration file please read the
[source term tutorial](./tutorial/source_term.html).

@todo pv: Correct the initialization description

There all possible variables that are actually notated
in the configuration file and added to `equation%sources`.
There all information from the configuration file
as well as from the expected ovariable list are contained in two arrays.
The first is `equation%sources%varname`,
which contains the names of the source variables.
The second array is called `equation%sources%variable`.
This array contains instances of type [[tem_variable_type]],
which itself contains a label, the number of components,
the type of the variable and a reference to the function
that provides the variable value when called.
The label is the same as the variable name, but prefixed with 'source\_'.
The components are not taken from the configuration file but
from the list of expected variables.
The type of the variable depends on how the source variable
is specified in the configuration file.
And last but not least the reference to the function is added,
which is also dependend on the variable notation in the configuration file.

After extracting all variables from the configuration,
we need to add them to the variable system.
Therefore we have to add the lua variable
which is called `source_<variablename>` to the variable system
and later on add the variable `<variablename>` to the varSys as well.
The latter one needs the former one as input varibale
and it also needs to define specific get\_element and get\_point routines
to actually get the data from the input varibale.

Additionally, the `equation_type` has a list `eval_source`,
which contains function pointers to the routines
that actually evaluate the source term and apply it to the right hand side.
To avoid unnecessary projections,
these evaluation routines are now called in with modal values
and have to do the projection internally, if needed.

### Material

The configuration of materials is done according to source terms. Each equation
system is able to use (or not to use) a different set of material parameters.
Therefore, the material configuration is a part of the equation configuration.
Inside the equation table a material table is to be defined. We handle
materials as variables. For each material parameter the equation system expects
a reference to a user-defined variable:

```lua
variable = {
  {
    name = "permeability",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { ... }
  },
  {
    name = "permittivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { ... }
  },
  {
    name = "conductivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { ... }
  }
}

equation = {
  name = 'maxwell',
  material = {
    permeability = 'permeability',
    permittivity = 'permittivity',
    conductivity = 'conductivity'
  }
}
```

Following this approach, every material parameter is defined by at least one
space-time function. However, there are situations where more than one
space-time function is needed. This is, for example, when there are different
materials for the background and a special region within the domain. To get
different materials for different regions, it is possible to add more than one
space-time function to a material variable:

```lua
variable = {
  {
    name = 'permeability',
    ncomponents = 1,
    vartype = 'st_fun',
    evaltype = 'first',
    st_fun = {
      {
        const = { 1.0 }
      },
      {
        const = { 2.0 },
        shape = {
          kind = 'canoND',
          object = {
            origin = { -0.5, -0.5, -0.25 },
            vec = {
              { 0.5 - scatter_eps, 0.0, 0.0 },
              { 0.0, 0.5 - scatter_eps, 0.0 },
              { 0.0, 0.0, 0.5 - scatter_eps },
            },
            segments = { 100, 100, 100 }
          }
        }
      }
    }
  }
}
```

With more than one space-time function for a variable, it is possible to get
more than one value for one point. Therefore the evaltype setting was
introduced. It defines how the evaluation routines should handle multiple
values. Possible values are `add` (default), `first` and `firstonly_asglobal`.

* For variables with the evaltype `add`, all values for a given point are added.
  Example: There is a background material with the constant value 1 and there
  is a special shape with the constant value 2. For all points outside the shape
  the variable returns 1, for all points inside the shape it returns 3.
* If the variable was defined with the evaltype `first`, the value of a point
  depends on the space-time function's definition order. The example from
  above (first space-time function with global shape returns 1, second
  space-time function for a maller shape returns 2) will result in all points
  having the value 1 as it is the first space-time function being defined.
  When the order is reversed, the points within the shape will have a value of
  2, all points outside the shape keep their value 1.
* TODO: The evaltype `firstonly_asglobal`

After reading all variables providing data for material parameters, Ateles will
loop over all the subtrees of the different material variables and assign that
space-time function to the element it has to evaluate. In the end, every element
knows which space-time function it has to ask for the value of a single material
parameter. This means, that for maxwell, each element has three space-time
functions it needs to ask for the data, one for permeability, one for
permittivity, and one for conductivity.

The user has to take care of defining all material parameters for the whole
domain. With the mechanism of different priorities based on the sorting
described above, the user can define a global material parameter by providing it
as the last space time function with global shape. All other space-time
functions with a non-global shape will act according to the selected evaltype.

After initializing all material space-time functions and assigning them to the
elements, Ateles will check for elements that have no space-time function
assigned (which means, the index to the space time function is 0). If there is
at least one element that has at least one space-time function missing, Ateles
will stop with a meaningful error message.

For implementation details on materials, see [Material](./pages/material.html).
