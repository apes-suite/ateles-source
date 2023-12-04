title: Gaussian Density Pulse in 2D

This setup for 2D Euler equations implements a circular Gaussian pulse in
density that is advected along with the fluid that has a constant velocity
in the complete domain.
All boundaries in this domain are periodic, however a sponge is applied
close to the downstream boundary that pulls the fluid to a target state.
The sponge is implemented as a source term that appears on the right side
of the equations.

The configuration is found in `ateles.lua`:

```lua
{!examples/euler/2D/gauss_densitypulse_spongelayer/ateles.lua!}
```

**Features used**

1. Projection: fpt, Blocksize 32

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --

6. Others:
   - Source term: Spongelayer
