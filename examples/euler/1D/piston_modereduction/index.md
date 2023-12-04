title: Formation of a shock wave ahead of a moving piston in 1D

This setup illustrates a small example on the setup of a moving geometry.
The fluid is initially in rest. A piston (geometry) is located in the middle 
area of the domain that has a length of 1.0. The piston moves with Mach 0.4 towards 
the right boundaty, due to the sudden movement, a shock is formed ahead of the piston.
Behind the piston rarefaction can be observed. The piston is modelled 
as a porous medium. An exact solution can be found in literature for this test 
case e.g. in Toro. This test case is used to ensure conservation of mass, momentum
and energy. If a shock is not formed ahead of the piston or it has a wronge velocity,
conservation is not maintained. 

In this test case we reduce the computational cost inside the geometry, where the 
solution is not of interest by means of the modereduction feature. Elements covered
by the geometry with neighbouring elements covered by the geometry as well, compute 
only the zero mode for the physical flux. The the zero mode is the integral mean of 
the Legandre polynomials.

The complete configuration is provided in `ateles.lua`:

```lua
{!examples/euler/1D/piston_modereduction/ateles.lua!}
```

1. Projection: l2p

2. Polynomial representation: Q

3. Filtering: spectral_viscosity 

4. Timestepping: imexRungeKutta, 4 steps

5. Boundary conditions: outflow

6. Others: - porous material (geometry), 
           - over-integration for the geometry (Piston),
           - modereduction
           - lua function for geometry definition
            
