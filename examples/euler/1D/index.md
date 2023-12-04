title: 1D Euler Equations

The 1D Euler equations represent instationary, inviscid, compressible flows in
a single space dimension, typically to model states in a tube.
It is configured by setting the name in the equation table to `euler_1d`.

Here is an example for the equation table of this equation system:

```lua
  equation = {
    name      = 'euler_1d',
    isen_coef = 1.4,
    r         = 287,
    numflux   = 'hll',
    material = {
      characteristic = 0.0,
      relax_velocity = 0.0,
      relax_temperature = 0.0
    }
  }
  equation.cv = equation.r / (equation.isen_coef - 1.0)
```

The fluid is described by the ideal gas constant (`equation.r`), the
isentropic expansion coefficient (`equation.isen_coef`) and the isochoric
specific heat (`equation.cv`).
See [[atl_eqn_euler_module]] for all options.

Note: you have to use the `modg_1d` scheme to compute this equation
system (`scheme.spatial.name = 'modg_1d'`).

The following example setups are available:

* [Modal Estimate](modalEstimate): illustrates the use of a modal estimation
  for the adaptive time step computation based on the CFL condition.

* [Shock Stabilization](shock_stabilization_parallel): demonstrates the use
  of a spectral and covolume filter to stabilize the simulation for a single
  moving shock.

* [Toro](toro1_x): utilizes a Riemann problem setting and simulates a shock
  tube.

* [Piston](piston): shock formation due to the sudden movement of a piston 
  geometry.

* [modereduction](piston_modereduction): shock formation due to the sudden movement
  of a piston geometry. Reduced computation of the physical fluxes inside the piston
  through the mode reduction feature.
