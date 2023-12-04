title: Shear-Layer Q4

This setup utilizes Q-polynomial representation and a spatial scheme order of 4.
It simulates a velocity layer in a rectangular channel with slipwall boundaries.
See the `ateles.lua` for the configuration:

```lua
{!examples/euler/3D/shear_layer_Q4/ateles.lua!}
```

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: slipwall
