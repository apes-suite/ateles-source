title: Constant State

This is a most basic setup for the Navier-Stokes 2D equations.
It initializes the domain with a constant state with periodic boundary
conditions.
The setup checks for fundamental flaws in the computation, as the constant
state should be maintained by the simulation.

The configuration file `ateles.lua` describes the complete simulation:

```lua
{!examples/navierstokes/2D/constant_state/ateles.lua!}
```

1. Projection: fpt, Oversampling 2

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: explicitEuler

5. Boundary conditions: -

6. Others: 
   - Derived quantity: Momentum, density and energy

