title: Example setups

This is a collection of examples that illustrates the usage of the various
capabilities of Ateles.

Ateles is configured via Lua scripts which need to define several variables.
Most of these variables are tables with multiple components. See any of the
specific examples below for a complete configuration.

At least the following need to be defined:

* `sim_control`, see [[tem_simControl_module]]
* `mesh`, see [[treelmesh_module]]
* `equation`, see [[atl_equation_module]]
* `scheme`, see [[atl_scheme_module]]
* `projection`, see [[atl_load_project_module]]
* `initial_condition`, see [[atl_initial_condition_module]]
* `boundary_condition` (if there are boundaries in the mesh),
  see [[atl_bc_header_module]]

Some other variables may be set to enable optional features or override
defaults:

* `tracking`, see [[ply_sampled_tracking_module]]
* `restart`, see [[atl_restart_module]]
* `logging`, see [[tem_logging_module]]
* `check`, defines how to check for unphysical states,
           see [[atl_physCheck_module]]

Treelm also provides various general settings that may be specified in the
configuration, see [[tem_general_module]].

Please note that you can include other Lua scripts with
[require](https://www.lua.org/pil/8.1.html).
And you can access table components with a dot notation like `equation.name`.

The Lua script will be executed by Ateles and in the end the defined variables
will be used as configuration for the simulation.

A configuration that shows the typical parameters for a flow simulation is
provided in [Euler 3D overview](euler/3D/overview).

This collection of examples is organized by equation systems, please see the
respective subdirectories for specific configuration examples and details for
specific equation systems:

* [(compressible) Navier-Stokes](navierstokes)
* [(inviscid) Euler](euler)
* [Linearized Euler](lineareuler)
* [Acoustic wave propagation](acoustic)
* [Maxwell](maxwell)
* [Maxwell with perfectly matched layers](maxwellpml)
* [Maxwell with divergence correction](maxwelldivcorr)
* [Heat conduction](heat)
