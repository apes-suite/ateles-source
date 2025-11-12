title: Toro4 setup for 2D Euler equations with a positivity preserving filter

In contrast to the other Toro examples, this setup uses a higher-order
discretization and the positivity preserving filter to simulate the
the Riemann problem with a strong discontinuity and the shock running
from right to left.
The positivity preserving filter on its own does not yield a stable
simulation, a spectral viscosity needs to be applied in addition.
An oversampling by a factor of 3 is used and in the projection we use
Chebyshev-Lobatto points in the nodal representation (as required by
the positivity preserving filter).

The complete configuration is provided in `ateles.lua`:

```lua
{!examples/euler/2D/toro4_x/ateles.lua!}
```

**Features used**

1. Projection: l2p, chebyshev-lobatto

2. Polynomial representation: Q

3. Filtering: `cons_positivity_preserv`, `spectral_viscosity`

4. Timestepping: explicitSSPRungeKutta, 2 steps

5. Boundary conditions: `inflow_normal`, `outflow`
