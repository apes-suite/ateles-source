title: Shear hat 2D

This example uses a initial velocity field with linear velocity increase on
both sides of the x-axis and a very high viscosity.
It illustrates the effect of viscosity.

The configuration in `ateles.lua` describes the setup:

```lua
{!examples/navierstokes/2D/shear_hat/ateles.lua!}
```

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: -

6. Others: 
   - Derived quantity: Momentum, density and energy

