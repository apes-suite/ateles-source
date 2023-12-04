title: Introducing a New Equation System

Here, we describe the way
to introduce a new equation system into the Ateles framework.
For this purpose we are going to describe the modules,
which need to be created or extended to introduce the one dimensional
[Acoustic wave equation](https://en.wikipedia.org/wiki/Acoustic_wave_equation).
Hence this is a simple example,
you will find comments wherever extension is required for more complex cases.

[TOC]

## Equation Specific

First we have to create a module to describe the acoustic equation system
and in addition we also have to create some supporting helper modules.

*   [[atl_eqn_acoustic_module]]
    \- Module to describe the acoustic equation system
    *   [[atl_eqn_acoustic_module:atl_acoustic_type]]:
        Datatype for the acoustic equation system.
    *   [[atl_eqn_acoustic_module:atl_load_acoustic]]:
        Subroutine to get the required values as background velocity
        from the lua configuration file.
*   [[atl_eqn_acoustic_var_module]]
    \- Module to configure the variables of the acoustic equation
    *   [[atl_eqn_acoustic_var_module:atl_init_acoustic_vars]]:
        Subroutine to initialize which variables has to be used
        and how they are organized in memory.
    *   [[atl_eqn_acoustic_var_module:atl_append_acoustic_vars]]:
        Subroutine to append the variable system by the variables
        for acoustic simulations using
        [[tem_varSys_module:tem_varSys_append_stateVar]].
    *   [[atl_eqn_acoustic_var_module:atl_init_acoustic_sourceTerms]]:
        Subroutine to init source terms for flow simulations
        to apply a source term like sponge layer.
*   [[atl_eqn_acoustic_hlp_module]]
    \- Module to define some helper routines
    *   [[atl_eqn_acoustic_hlp_module:atl_eqn_acoustic_init]]:
        Subroutine to set up the necessary infrastructure.
    *   [[atl_eqn_acoustic_hlp_module:atl_eqn_acoustic_load_bc]]:
        Subroutine for reading the boundary conditions.

The location of the equation specific modules is 'source/equation'.
As you may have noticed the name of the modules
related to the specific eqaution system should have the prefix 'atl\_eqn\_'.
The postfix '\_var\_module' indicates a module in which
the variables of the different equation systems should be defined,
the postfix '\_hlp\_module' indicates a module
in which helping routines are defined
and last but not least the postfix '\_derive\_module' indicates a module
in which routines are defined to calculate derived quantities
of the different equation systems.

To use the subroutines and datatypes of the prescribed modules
we have to extend the module [[atl_equation_init_module]]
by the acoustic testcase to call for
[[atl_eqn_acoustic_hlp_module:atl_eqn_acoustic_init]].
Furthermore [[atl_equation_module:atl_Equations_type]] needs to be extended
by the new defined datatype [[atl_eqn_acoustic_module:atl_acoustic_type]].

## Flux Calculation

To calculate the physical and numerical flux in one dimension
we create the modules [[atl_acoustic_physflux_module]]
and [[atl_acoustic_numflux_module]] respectively in 'source/flux'.
For the acoustic equation, the numerical flux can be calculated analytical.
For nonlinear equation,
the numerical flux is calculated by riemann problem solver
like LaxFriedrich, HLL or Roe.

## Kernel for DG

For the routines and datatypes
of the Modal Discontinuous Galerkin (MODG) scheme for the acoustic equation
we create the module [[atl_modg_acoustic_kernel_module]]
in 'source/scheme/modg'.
This is the main module where the 'numerics' happen.
Hence, the calls to the one dimensional numerical
and physical flux calculation are done.
Take care about the correct rotation
according to the corresponding direction of the flux.
For nonlinear equation,
the physical flux is calculated in the nodal space,
though further projection (modal to nodal and vice versa) are required.

## Penalization and Material Properties

For the penealization and material properties we have to extend some modules.

*   [[atl_eqn_acoustic_module]]:
    In [[atl_eqn_acoustic_module:atl_acoustic_type]]
    we have to make sure to include
    [[atl_materialFun_module:atl_materialFun_type]]
    for the penalization functions.
    This is very important for the penalization module
    [[atl_penalization_module]]
    and the overall computation module [[atl_compute_module]],
    since the elements are split into elements
    with const and non_const material properties
    and these elements define the overall compututation loop.
*   [[atl_materialIni_module]]:
    In the subroutine [[atl_materialIni_module:atl_init_materialParams]]
    the case for the new equation need to be added.
    In this routine the elements are sorted by their material properties
    and the correct list for the later on compute module are created
    (without these lists, the loops in the compute module are not iterated)

Without material parameters,
the computatation does not run into the physical flux calculation.
The material and penalization specific modules are in 'source/material'.

## Computation

For the computation some modules needs to be extended.

*   [[atl_compute_module]]
    \- Module which controls the computation.
*   [[atl_global_time_integration_module]]:
    \- Module with routines to steer the time stepping mechanism in Ateles.
*   [[atl_compute_local_module]]
    \- Module with routines to to calculate element-local RHS.
    The preprocess routine
    [[atl_compute_local_module:preprocess_local_rhs_cubes]]
    needs to be extended by the acoustic case.
*   [[atl_calc_time_module]]
    \- Module which provides the definition
    and methods for time step computations.
*   [[atl_physcheck_module]]
    \- Module for checking physical values of a given state.
    Need to be extended by checking for state variables
    of the new variable system.
