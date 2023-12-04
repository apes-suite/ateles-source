title: Pulse in density with derived quantity

This setup for 3D Euler equations features a Gaussian pulse in density
that is transported through the domain.
It illustrates the tracking of kinetic energy as a derived quantity with
the variable system of the solver.

The configuration is provided in `ateles.lua`:

```lua
{!examples/euler/3D/gauss_densitypulse_deriv_quantity/ateles.lua!}
```

**Features used**

1. Projection: fpt, blocksize 32

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --

6. Others:
   - Derived quantity: `velocity`, `kinetic_energy`

