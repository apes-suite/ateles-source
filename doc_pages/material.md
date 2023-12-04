title: Material

Material is a concept to assign different parameters to different areas of the
domain.

@TODO Extend the description

[TOC]

## Initialization

For the initialization of the material, please refer to
[Initialization](./init.html#material).

## Data structures

Information about material are stored in two main data structures. One is
[[atl_material_type]], the other is [[atl_materialFun_type]].
[[atl_materialFun_type]] is
contained in equations_type and stores information about the lua variables the
material consists of, their representation in the global variable system and
the number of scalars of each variable. [[atl_material_type]] is located in
[[atl_element_container_type]]%[[atl_cube_container_type]] and is stored as an
levelwise array. This type knows about the elements, faces, boundaries affected
by it and is used to evaluate the material during runtime.

## Implementation

The variable table is read in [[atl_initialize]] by calling
[[tem_variable_load]]. The variables are stored in the local array userVars of
type [[tem_variable_type]]. All variables in this array are then given to
[[tem_varSys_append_luaVar]], where they are added to the global variable
system. The space-time functions are added to the global list of space-time
functions stored as [[atl_equations_type]]%[[atl_equations_type:stFunList]].
The variable's method_data is set to point to the corresponding space-time
function. Doing his, the evaluation of a variable based on a space-time
function can just dereference it's method_data to have the correct space-time
function available. Right after that, the element lists for all space-time
functions in [[atl_equations_type]]%[[atl_equations_type:stFunList]] are
created. Based on these lists Ateles afterwards has to assign the latest valid
material to the elements.

### Reading the configuration file

Each equation system expects a set of material and penalization parameters.
Some of these parameters are constant in time and space, some may vary. The
constant material parameters are explicitly defined in the equation table, e.g.
the isentropic coefficient in navier stokes:

```lua
equation = {
  name = 'navier_stokes',
  ...
  isen_coef = 1.4,
  ...
}
```

As long as there is no need to vary those parameter in time or space, they can
stay as an equation global constant. But all other parameters that are not
constant have to be defined using variables. As these variables are configured
in the variable table, we need a connection between a user defined variable in
the variable table and a material parameter. Therefore each equation system has
to announce which materials it can make use of. During the equation
initialization, an instance of [[atl_init_material_type]] is filled. It
contains a property [[atl_init_material_type:poss_materialVars]] of type
[[tem_possible_variable_type]], where all material parameters the equation can
use are stored in.

This list is used in [[tem_variable_loadMapping]] to find the material
parameters in the equation table. If they are found, the material parameter's
name is stored along with the name of the corresponding variable in the
variable table in a dictionary
[[atl_init_material_type]]%[[atl_init_material_type:materialDict]].

### Storing the data in the global variable system

After loading these mappings, the [[atl_init_material_type]] is passed to
[[atl_append_newMaterialVars]] to add a new variable to the global variable
system for each configured material parameter. This variable gets the name of
the material parameter prefixed with mat_. Following this naming scheme, the
user is able to track the material parameter as well as the space-time function
providing the data for this material parameter. The input variable for this new
variable is the variable that was added for the space-time function. The
position of the new material variable is stored along with the number of
components this variable provides in an instance of [[atl_materialFun_type]].
This list is stored inside the [[atl_equations_type]], thus accessible during
runtime for easy access to all material variables.

That said, the sequence to get data for a material parameter looks like:

1. call `varSys%method%val(<position of mat_permeability>)%get_point`
2. inside `mat_permeability%get_point`: call for the one and only input variable
3. inside `input_variable%get_point`: Cast the `method_data` into a space-time
   function and evaluate it for the requested point(s)

Looking at the second step, one can see the reason why we need a special
get_point and get_element routine for materials. Both routines
([[atl_getMaterialForPoint]] and [[atl_getMaterialForElement]]) are defined in
[[atl_materialIni_module]] and only forward requests on the material variable
to the variable representing the space-time function.

### Preparing the data structures for production

When all variables are added to the variable system, the initialization of the
data structures used for production has to be done. Therefore the routine
[[atl_assign_elem_varProp]] iterates over all material variables, gets the
corresponding Lua variable, extracts the space-time functions and assigns it to
each element in the space-time function's subtree. This is done by using a
pointer to the space-time function. As we usually have more than one material
parameter and thus don't have a single space-time function providing data for
it, we end up with a small list of pointers to space-time functions with the
length of the number of material parameters.

### Evaluating material parameters during runtime

**Material parameters are currently evaluated during runtime only for 2d
kernels.** This is done in
[[atl_modg_2d_kernel_module:atl_preprocess_modg_2d_kernel]], where
[[atl_update_materialParams]] is called. This routine calls the three steps of
updating material parameters:

- The first step is performed in [[atl_evalElemMaterial]]. There the material
  for each element is updated by evaluating the space-time function-pointer
  stored for each element and each material parameter. The new material data is
  stored for each element and material parameter in
  [[atl_material_type]]%[[atl_material_type:material_dat]].
- After this evaluation, the material is projected to the element's faces by
  calling [[atl_evalFaceMaterial]]. The face material is also stored in
  [[atl_material_type]]%[[atl_material_type:material_dat]].
- Now [[atl_create_materialBoundaryList]] is called. There the lists of
  boundary elements with material is created and stored in
  [[atl_material_type]]%[[atl_material_type:material_desc]]
