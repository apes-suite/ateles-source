title: 2D Linearized Euler Equations

Examples for the 2D linearized Euler equations.
This equation system is configured by setting the equation name
to `'lineareuler_2d'`.
The full equation table looks as follows:

```lua
  equation = {
    name   = 'linearEuler_2d',
    isen_coef = 1.4,
    numflux = 'godunov',
    background = {
      density = 1.225,
      velocityX = 100.0,
      velocityY = 0.0,
      pressure = 100000.0
    }
  }
```

In the linearized Euler equations the fluid needs to be defined with a
background state around which the linearization is done.

Available setups are:

* [Gaussian pulse](gauss_pulse): basic example with a density pulse that gets
  transported by a constant flow state.
* [Transient background](transientBackground): it is possible to vary the
  background state over time. This shows how it is done.
