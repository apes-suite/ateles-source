title: Euler Equations

This is a collection of setups for the Euler equations of inviscid compressible
flows.

See the [[atl_eqn_euler_module]] for details.

The Euler equations are implemented in [one](1D) (`equation.name = 'euler_1d'`),
[two](2D) (`equation.name = 'euler_2d'`) and [three](3D)
(`equation.name = 'euler'`) space dimensions.
See these subdirectories for complete examples of the respective equation
systems.

In the gasdynamics modeled by this equation system, discontinuities can arise,
leading to oscillations in high-order discretizations.
These oscillations can result in unphysical states that abort the simulation.
To counter this problem, it may be necessary to employ some stabilization as
provided by the [[atl_stabilization_module]].
A simple example for a minimal spectral filter and covolume is provided in
the [shock stabilization parallel setup](1D/shock_stabilization_parallel).

A template configuration that illustrates the various options to use in Euler
simulations is given in the [3D overview example](3D/overview).
