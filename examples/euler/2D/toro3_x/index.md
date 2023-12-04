title: Toro3 Setup along X for 2D Euler equations

Like the Toro2 setup, this configuration simulates a one dimensional Riemann
problem and uses the Godunov numerical flux with a first order scheme (FV).
The specialty of this configuration is a very strong jump in pressure in the
initial condition.

The complete configuration is provided in `ateles.lua`:

```lua
{!examples/euler/2D/toro3_x/ateles.lua!}
```

**Features used**

1. Projection: l2p

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 2 steps

5. Boundary conditions: `inflow_normal`, `outflow`
