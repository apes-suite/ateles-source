title: Euler 3D pulse in density with FPT

This example, solving the 3D Euler equations, simulates a spherical Gaussian
pulse in density, that is convected throug a periodic domain.
As a projection method the Fast Polynomial Transform (FPT) is used here.

The configuration is provided in `ateles.lua`:

```lua
{!examples/euler/3D/gauss_densitypulse_fpt/ateles.lua!}
```

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --
