title: Euler 3D setup with pulse in density and L2 projection

This setup employs the pulse in density, convected through a periodic domain
to show the use of the L2 projection method for transformations between modal
and nodal representations.

The configuration is found in `ateles.lua`:

```lua
{!examples/euler/3D/gauss_densitypulse_l2p/ateles.lua!}
```

**Features used**

1. Projection: l2p, Oversampling 2.0

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --
