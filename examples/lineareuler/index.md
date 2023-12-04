title: Linearized Euler Equations

This is a collection of setups for the linearized Euler equations.
The definition of linearized Euler equations is provided in
[[atl_eqn_lineareuler_module]].
The examples are further distinguished into setups for [three](3D) and
[two](2D) dimensions.

For the linearized Euler equation we need to provide a background state
around which the linearization is to be done.
This is given in terms density, velocity and pressure.
The equation definition, therefore, looks as follows:

```lua
  equation = {
    name   = 'linearEuler',
    isen_coef = 1.4,
    background = {
      density = 1.225,
      velocityX = 100.0,
      velocityY = 0.0,
      velocityZ = 0.0,
      pressure = 100000.0
    }
  }
```
