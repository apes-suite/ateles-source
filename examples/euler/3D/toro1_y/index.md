title: Toro1 Riemann problem in Y direction

Setup of a channel in Y direction. The initial condition prescribes a
Riemann problem.

Find the configuration in `ateles.lua`:

```lua
{!examples/euler/3D/toro1_y/ateles.lua!}
```

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: `slipwall`, `inflow_normal`, `outflow`
