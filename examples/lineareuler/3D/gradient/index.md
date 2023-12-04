title: Tracking of gradients in 3D linearized Euler equations

This setup of linearized 3D Euler equations illustrates the computation of
gradients in the variable system for tracking.
See the 'variable' table.

We employ here a Taylor Runge-Kutta time integration scheme. This scheme allows
the use of an arbitrary number of stages and yields an according high order for
linear, autonomous systems. It enables the explicit simulation to achieve
timesteps with a Courant factor greater than 1.

The configuration is found in `ateles.lua`:

```lua
{!examples/lineareuler/3D/gradient/ateles.lua!}
```

1. Projection: l2p

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: explicitRungeKuttaTaylor, 8 steps
