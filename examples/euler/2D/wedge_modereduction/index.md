title: Supersonic movement of a wedge

This setup is an academic test case. A sharp wedge moves with
supersonic speed from the right boundary to the left. The deflection
angle of the wedge is 15 degree. The fluid is initially in rest, due 
to the movement of the wedge a shock is formed ahead of the wedge. The
shock angle can be computed analytically. For this setup, the exact shock
angle is 43.345 degree. 

The coordinates are provided to the solver via a list of points. 

The configuration is found in `ateles.lua`:

```lua
{!examples/euler/2D/wedge/ateles.lua!}
```

**Features used**

1. Projection: l2p

2. Polynomial representation: Q

3. Filtering: spectral viscousity and covolume

4. Timestepping: imexRungeKutta, 4 steps

5. Boundary conditions: --

6. Others:
   - Geometry 
   - Geometry defined through Fortran function
   - Modereduction inside the geometry
