title: 2D Navier-Stokes Equations

The 2D Navier-Stokes equations represent instationary, viscous, compressible
flows in two spatial dimensions.
It is configured by setting the name in the equation table to
`navier_stokes_2d`.

Here is an example for the equation table of this equation system:
```lua
  equation = {
    name      = 'navier_stokes_2d',
    isen_coef = 1.4,
    r         = 287,
    -- Viscous parameters
    therm_cond = 0.5,
    mu         = 1.e-5,
    ip_param   = 8/3,
    material = {
      characteristic = 0.0,
      relax_velocity = {0.0, 0.0},
      relax_temperature = 0.0
    }
  }
```

Note: you have to use the `modg_2d` scheme to compute this equation
system (`scheme.spatial.name = 'modg_2d'`).

The following example setups are available:

* [constant state](constant_state): represents the simplest possible setup
  with a constant state and periodic boundary conditions. It mainly serves
  as a check to ensure there is nothing fundamentally flawed in the
  implementation.

* [shear hat](shear_hat): provides a small example with extreme viscosity and
  an initial hat velocity profile in y direction (linear y-velocity profile
  left and right of the x-axis). It serves as a check on the treatment of
  viscous effects in the scheme.

* [viscous vortex](viscous_vortex): provides a 2D vortex setup with additional
  source terms that yield a known analytical solution to the system.
