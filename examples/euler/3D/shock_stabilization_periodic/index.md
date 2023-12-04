title: Shock stabilization in periodic domain

Setup of a travelling planar shock in 3D Euler equations, stabilized by spectral
and covolume filtering.
The configuration is found in `ateles.lua`:

```lua
{!examples/euler/3D/shock_stabilization_periodic/ateles.lua!}
```

**Features used**

1. Projection: l2p, Oversampling 2

2. Polynomial representation: Q

3. Filtering: spectral, covolume

4. Timestepping: explicitSSPRungeKutta, 2 steps

5. Boundary conditions: --
