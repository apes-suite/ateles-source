title: Modal estimate for adaptive timesteps 1D

This setup illustrates the use of modal estimation for the adaptive timestep
according to the CFL condition in the 1D Euler equations.
The adaptive timestep requires the computation of maximal velocities in the
domain and involves an expensive modal to nodal transformation.
With the modal estimation shown in this example it is possible to avoid this
transformation.
However, the estimation is not very accurate and the resulting time step may
become very small with this approach.

The complete configuration is provided in `ateles.lua`:

```lua
{!examples/euler/1D/modalEstimate/ateles.lua!}
```

1. Projection: l2p

2. Polynomial representation: -

3. Filtering: -

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: -

6. Others: `use_modal_estimate`
