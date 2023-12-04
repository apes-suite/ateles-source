title: Vorticity Tracking

This setup illustrates the tracking of derived quantities like velocity and
vorticity in a single point of the domain.

The complete configuration in `ateles.lua`:

```lua
{!examples/euler/2D/vorticity/ateles.lua!}
```

**Features used**

1. Projection: l2p

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --

6. Others: 
   - Derived quantity: velocity, vorticity
