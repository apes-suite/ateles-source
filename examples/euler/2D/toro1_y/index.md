title: Toro1 Validation Test Case in Y direction

This a validation testcase for Euler 2D. It simulates a shocktube in
y direction, with different, constant states left and right in the shocktube,
resulting in a discontinuity in the middle.
Thus, a Riemann problem is to be simulated, for which we now the solution and
can compare the result.
In this case there is no internal mesh available to represent a line of elements
along the y axis. The setup uses a seeder generated mesh instead.

The configuration is provided in the `ateles.lua` file:

```lua
{!examples/euler/2D/toro1_y/ateles.lua!}
```

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: `slipwall`, `inflow_normal`, `outflow`
