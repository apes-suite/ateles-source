title: Viscous Vortex (TGV) 2D

This 2D Navier-Stokes setup utilizes an initial vortex configuration with a
source term such, that an analytical solution is known, which enables a check
for correctness of the simulation.

The configuration is found `ateles.lua`:

```lua
{!examples/navierstokes/2D/viscous_vortex/ateles.lua!}
```

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: -

6. Others: 
   - Derived quantity: Momentum
