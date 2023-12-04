title: Q-Criterion Tracking in Euler 3D

This setup shows the tracking of q-criterion and lambda2 values in a point of
a periodic domain.

```lua
{!examples/euler/3D/qCriterion_deriv_quantity/ateles.lua!}
```

**Features used**

1. Projection: fpt, blocksize 32

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --

6. Others:
   - Derived quantity: `Q_criterion` and `lambda2`

