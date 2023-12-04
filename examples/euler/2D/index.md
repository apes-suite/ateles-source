title: 2D Euler Equations

The 2D Euler equations model inviscid, compressible flows in two spatial
dimensions.
It is configured by setting the name in the equation table to `euler_2d`.

Here is an example for the equation table of this equation system:

```lua
  equation = {
    name      = 'euler_2d',
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

Note: you have to use the `modg_2d` scheme to compute this equation
system (`scheme.spatial.name = 'modg_2d'`).

The following example setups are available:

* [Gaussian pulse in density](gauss_densitypulse_spongelayer): Simple setup with
  a circular pulse that gets convected into a sponge layer. Illustrates the
  definition of a sponge.

* [Modal Estimate for timestep restriction](modalEstimate): Illustration of the
  use of modal estimation for the adaptive timestep that uses a pressure pulse
  in the initial condition.

* [Shock Stabilization](shock_stabilization_parallel): Demonstrates the use
  of a spectral and covolume filter to stabilize the simulation for a planar
  moving shock.

* [Toro 1 in x](toro1_x)

* [Toro 1 in y](toro1_y)

* [Toro 2](toro2_x)

* [Toro 3](toro3_x)

* [Toro 4](toro4_x)

* [Vorticity](vorticity)

* [Wedge](wedge): Illustrates the supersonic movement of a sharp geometry (wedge).
* [Modereduction](wedge_modereduction): Illustrates the supersonic movement of a 
  sharp geometry (wedge). Elements inside the wedge are reduced in computation using 
  the modereduction feature.
