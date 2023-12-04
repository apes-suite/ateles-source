title: Toro1 Riemann problem in X direction

Setup of a channel in X direction. The initial condition prescribes a
Riemann problem.

Find the configuration in `ateles.lua`:

```lua
{!examples/euler/3D/toro1_x/ateles.lua!}
```

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: `slipwall`, `inflow_normal`, `outflow`
