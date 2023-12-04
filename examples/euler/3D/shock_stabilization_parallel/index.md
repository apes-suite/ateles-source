title: Shock stabilization for 3D Euler equations with RK-Taylor

This setup simulates a shock moving along the z-axis, highlighting covolume and
spectral filters for stabilization.
The time integration uses the Runge-Kutta Taylor scheme.

The configuration is provided in `ateles.lua`:

```lua
{!examples/euler/3D/shock_stabilization_parallel/ateles.lua!}
```

**Features used**

1. Projection: l2p, Oversampling 2

2. Polynomial representation: Q

3. Filtering: spectral, covolume

4. Timestepping: explicitRungeKuttaTaylor, 4 steps

5. Boundary conditions: `slipwall`, `inflow_normal`, `outflow`
