title: 3D Navier-Stokes Equations

The 3D Navier-Stokes equations represent instationary, viscous, compressible
flows in three spatial dimensions.
They are configured by setting the name in the equation table to
`navier_stokes`.

The equation table for this system may look as follows:
```lua
equation = {
  name      = 'navier_stokes',
  isen_coef = 1.4,
  r         = 287,
  -- Viscous parameters
  therm_cond = 0.5,
  mu         = 1.e-5,
  ip_param   = 8*(degree+2)/(2*(degree+3)),
  material = {
    characteristic = 0.0,
    relax_velocity = {0.0, 0.0, 0.0},
    relax_temperature = 0.0
  }
}
```

Note: you need to use the `modg` scheme to solve this equation system
(`scheme.spatial.name = 'modg'`).

The following example setups are available:
* [shear hat](shear_hat): a simple setup with an initial linear y-velocity
  profile left and right of the x-axis (no variation in z) with periodic
  boundary conditions all around.
* [shear tube](shear_tube): this prescribes a cylindrical jet along the
  x-axis in a medium at rest and checks thereby the momentum transfer across
  the resulting shear layer. Again, boundary conditions are all periodic.
* [Taylor-Green Vortex](tgv): the Taylor-Green Vortex with a fully periodic
  cubical domain
