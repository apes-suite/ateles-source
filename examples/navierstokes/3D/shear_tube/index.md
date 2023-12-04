title: Shear tube 3D

In this example a cylindric velocity jet along the x-axis is simulated in a
periodic domain with a fluid at rest.
It illustrates the use of a positivity preserving filtering to ensure the
stability of the scheme.
To ensure the positivity throughout the domain and over the simulation time,
a strong stability preserving time integration and Lobatto integration points
have to be used.
The configuration is found in `ateles.lua`:

```lua
{!examples/navierstokes/3D/shear_tube/ateles.lua!}
```

1. Projection: fpt, Oversampling 3.0, lobattoPoints

2. Polynomial representation: Q

3. Filtering: cons_positivity_preserv

4. Timestepping: explicitSSPRungeKutta, 2 steps

5. Boundary conditions: -

6. Others: 
   - Derived quantity: Momentum
