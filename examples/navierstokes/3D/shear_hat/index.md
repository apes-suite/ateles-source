title: Shear hat 3D

This configuration, given in `ateles.lua`, uses an initial velocity field
with linearly increasing velocity in y direction towards the x-axis.
The high viscosity amplifies the momentum transfer by shearing.
The domain is periodic in all directions.

```lua
{!examples/navierstokes/3D/shear_hat/ateles.lua!}
```

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: -

6. Others: material
