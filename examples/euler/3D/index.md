title: 3D Euler Equations

The 3D Euler equations model inviscid, compressible flows in two spatial
dimensions.
It is configured by setting the name in the equation table to `euler`.

Here is an example for the equation table of this equation system:

```lua
  equation = {
    name      = 'euler',
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
See [[atl_eqn_euler_module]] for all options, and the
[overview example](overview) for a configuration file with the usual options
that can serve as a template for your own configurations.

Note: you have to use the `modg` scheme to compute this equation
system (`scheme.spatial.name = 'modg'`).

Available setups are:

* [Pulse in density with derived quantity](gauss_densitypulse_deriv_quantity):
  shows the tracking of derived quantities like kinetic energy and velocity.

* [Pulse in density with FPT](gauss_densitypulse_fpt): simple setup with FPT
  projection between modal and nodal representation.

* [Pulse in density with FXT](gauss_densitypulse_fxt): simple setup with FXT
  projection between modal and nodal representation. This uses the fast
  multipole method, implemented in FXTPACK.

* [Pulse in density with L2P](gauss_densitypulse_l2p): illustrates the use
  of the L2 projection for transformations between modal and nodal
  representation.

* [Pulse in density with sponge](gauss_densitypulse_spongelayer): shows the use
  of a sponge to dampen out the state.

* [Modal Estimate](modalEstimate): small setup with a pressure pulse that shows
  the use of modal estimation for the adaptive timestep computation.

* [Multilevel](multilevel): modelling of a wall by penalization in a multilevel
  mesh. Employs the IMEX time integration scheme (needed for efficient handling
  of penalization terms).

* [Overview](overview): this example illustrates the overall options for Euler
  simulations and can serve as a template for your configuration files.

* [P polynomials](p_polynomials): a shear layer in a rectangular tube. This
  setup makes use of the P-Polynomial construction for the multidimensional
  polynomials.

* [Q-Criterion](qCriterion_deriv_quantity): shows tracking of q-criterion and
  lambda2 values.

* [Shear layer Q4](shear_layer_Q4): a simple setup with a velocity layer in
  a periodic channel discretized with fourth order in space.

* [Shear layer Q8](shear_layer_Q8): same as above but with an eighth order
  discretization in space.

* [Shock locally refined](shock_stabilization_localrefinement): a shock that
  travels in a tube along the z-axis through a domain with some local
  refinement. This setup illustrates covolume with spectral filtering.

* [Shock parallel](shock_stabilization_parallel): a shock travelling through
  a channel along the z-axis, simulated with the Runge-Kutta Taylor scheme.
  Utilizes covolume and spectral filtering to stabilize the simulation.

* [Shock periodic](shock_stabilization_periodic): a shock moving through a
  periodic domain, progressed in time with the strong stability preserving
  2-stage Runge-Kutta scheme.

* [Toro1 in X](toro1_x): a Riemann problem in a tube along the x-axis.

* [Toro1 in Y](toro1_y): a Riemann problem in a tube along the y-axis.
* [Toro1 in Z](toro1_z): a Riemann problem in a tube along the z-axis.
* [Track vorticity](vorticity)
