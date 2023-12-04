title: Toro2 Setup along X for 2D Euler equations

Like the Toro1 setup, this configuration simulates a one dimensional Riemann
problem. Here we make use of the Godunov flux with an exact Riemann solver for
the numerical fluxes.
The polynomial degree is reduced to 0, resulting in a Finite Volume
discretization with constant states in elements.
The domain is constructed by a line of elements along the x-axis defined by
the internal mesh definition `line_bounded`.
This domain has two boundary conditions, one at the western end and on at the
eastern end. All other directions are periodic.

The complete configuration is provided in `ateles.lua`:

```lua
{!examples/euler/2D/toro2_x/ateles.lua!}
```


**Features used**

1. Projection: l2p

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitSSPRungeKutta, 2 steps

5. Boundary conditions: `inflow_normal`, `outflow`
