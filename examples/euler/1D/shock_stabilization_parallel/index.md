title: Shock Stabilization Test Case

This 1D Euler setup simulates a single shock travelling from left to right
through the domain.
It makes use of a covolume filter with some weak spectral filtering, to overcome
oscillation residues at element boundaries after the shock passed them.
In the time discretization a Taylor Runge-Kutta scheme is employed, that allows
for an arbitrary number of stages.
The mesh is provided by the internally defined `line_bounded`, which provides
a line of elements with a `west` boundary condition at the left end and a
`east` boundary at the right. All other directions are periodic, but in 1D
anyway irrelevant.

The complete configuration is provided in the `ateles.lua` file:

```lua
{!examples/euler/1D/shock_stabilization_parallel/ateles.lua!}
```

**Features used**

1. Projection: fpt, Oversampling 2

3. Filtering: covolume, spectral

4. Timestepping: explicitRungeKuttaTaylor, 4 steps

5. Boundary conditions: `slipwall`, `inflow_normal`, `outflow`
