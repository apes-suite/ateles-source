title: Transient background state in 2D linearized Euler

This setup illustrates the use of time dependent functions as background state
in the linearized Euler equations.
The setup makes use of a single line of elements and no deviations from the
background state.
The densitiy deviation and the complete state (background + deviation) is
tracked in a single element.

```lua
{!examples/lineareuler/2D/transientBackground/ateles.lua!}
```

1. Projection: l2p

2. Polynomial representation: P

3. Filtering: -

4. Timestepping: explicitRungeKutta, 4 steps
