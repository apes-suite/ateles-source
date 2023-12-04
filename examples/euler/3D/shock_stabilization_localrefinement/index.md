title: Shock Stabilization in Euler 3D with local refinement

This setup illustrates the use of stabilization filters for the 3D Euler
equations in a domain with local refinement.
A shock travelling along the z-axis is simulated.
Such a travelling shock leads to oscillations that remain at element boundaries
and to eliminate those we employ a covolume filter with spectral viscosity.
This smoothes the oscillations and stabilizes the scheme.

The configuration of this setup is provided in `ateles.lua`:

```lua
{!examples/euler/3D/shock_stabilization_localrefinement/ateles.lua!}
```

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: covolume, spectral

4. Timestepping: explicitSSPRungeKutta, 2 steps

5. Boundary conditions: slipwall
