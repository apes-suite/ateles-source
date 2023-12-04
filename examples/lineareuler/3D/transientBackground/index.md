title: Transient background state in 3D linearized Euler

This setup illustrates the use of time dependent functions as background state
in the linearized Euler equations.
The setup makes use of a single line of elements and no deviations from the
background state.
The densitiy deviation and the complete state (background + deviation) is
tracked in a single element.

```lua
{!examples/lineareuler/3D/transientBackground/ateles.lua!}
```

Covered features:

* time dependent Lua function for background state
* FPT
* RK 4
* Q-Space polynomial
* tracked mean of element
* predefined mesh line
* Physcial Checks
