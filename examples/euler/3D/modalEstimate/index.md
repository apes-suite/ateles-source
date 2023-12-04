title: Modal estimate in Euler 3D

This setup illustrates the use of modal min/max estimation in computation
of timestep limitation for 3D Euler equations.
It simulates a pressure pulse in a fully periodic domain.

The `ateles.lua` file contains the configuration of the setup:

```lua
{!examples/euler/3D/modalEstimate/ateles.lua!}
```

1. Projection: l2p

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: -

6. Others: `use_modal_estimate`
