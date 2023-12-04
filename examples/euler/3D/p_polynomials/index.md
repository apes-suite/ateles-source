title: P-polynomial space

This setup illustrates the use of P-polynomial space in the scheme definition.
It simulates a layer of higher velocity fluid embedded in a fluid at rest inside
a rectangular channel, periodic in X.

P-Polynomials decrease the amount of memory required to represent the state,
without decreasing the spatial convergence order of the scheme.
It also decreases the amount of computations, but access patterns are less
optimal and many operations internally will resort to the full Q-Representation.
The advantage of P-Polynomials gets more pronounced the higher the polynomial
degree.

Find the configuration in `ateles.lua`:

```lua
{!examples/euler/3D/p_polynomials/ateles.lua!}
```


1. Projection: fpt

2. Polynomial representation: P

3. Filtering: -

4. Timestepping:  explicitRungeKutta, 4 steps

5. Boundary conditions: slipwall
