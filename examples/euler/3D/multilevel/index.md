title: Material with Multilevel in Euler 3D

This setup illustrates the use of material to model a wall in a multilevel
mesh.
For the Brinkman penalization we use to model the wall, we employ a IMEX
time integration scheme.
An acoustic pulse is simulated as it runs against a wall and gets reflected
by it. The domain uses a finer resolved mesh in the part of the domain with
the acoustic pulse, while the penalized part is discretized with a coarser
mesh resolution.
The same maximal polynomial degree is used in both parts.

Find the configuration in `ateles.lua`:

```lua
{!examples/euler/3D/multilevel/ateles.lua!}
```

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: imexRungeKutta, 4 steps

5. Boundary conditions: wall, outflow, periodic

6. Others: material, multilevel
