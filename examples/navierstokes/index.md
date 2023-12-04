title: Navier-Stokes Equations

This is a collection of setups for the compressible Navier-Stokes equations,
which allow for the simulation of viscous compressible flows.

See the [[atl_eqn_nvrstk_module]] for details.

We distinguish setups for [2D](2D) (`equation.name = 'navier_stokes_2d`)
and [3D](3D) (`equation.name = 'navier_stokes`).
See these subdirectories for complete examples for the respective
equation systems.

Due to the compressibility, strong gradients may form in the fluid, leading
to oscillations in high-order discretizations.
To stabilize simulations it may, therefore, be necessary to employ some
filtering techniques, see [[atl_stabilization_module]].
The [shear tube example](3D/shear_tube) provides a configuration with a
positivity preserving filter.
