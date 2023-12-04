title: Euler 3D pulse in density with sponge

Setup of a spherical pulse in density that is convected through a periodic
domain and dampened at a sponge.
The configuration shows how a sponge can be introduced in the simulation.
See the `ateles.lua` file for the complete configuration of the setup:

```lua
{!examples/euler/3D/gauss_densitypulse_spongelayer/ateles.lua!}
```

**Features used**

1. Projection: fpt, Blocksize 32

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --

6. Others:
   - Source term: Spongelayer
