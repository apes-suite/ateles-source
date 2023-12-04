title: Source Term

This example shows how to define a source term
for the [Maxwell testcase](maxwell.html).
Thus, you find the source code for this example commented out
in the configuration file of the Maxwell testcase.

In this Lua file, we define a source term with the keyword `source_terms`.
Moreover we need to give the solver a `shape` (point, line, or plane. See
[Tracking Shape](tracking.html#tracking-shape))
where the source terms will be active.
For all elements within this shape, the function `fun` will be calculated
and then changes the right hand side of the underlying equation system.

In the testcase we use a cube with a point source located in the middle.
This point source changes the current with time according to a sinus function.

At the beginning, we specify some variables which define the source
and will be used in the function.

```lua
-- ... charge of the source
Q = 1.0
-- ... radius of sphere source
r = 0.4
-- ... parameters for the analytic solution
freq = ( 2.0*math.pi/math.sqrt(mu*epsi) )
-- ... the temporal period of the waveguide
T = 2.0*math.pi/freq
```

Next we define the function for the point source.
In Maxwell equations, we are talking about space charges and current densities.
In general they can depend on spatial coordinates and time.
Here we define the current in each direction
with a sinus function with the frequency defined as above.

```lua
-- function for the point source
function currentDensitySpaceTime(x, y, z, t)
  d = math.sqrt(x^2.0+y^2.0+z^2.0)
  if d <= r then
    jx=Q*r*freq*math.sin(freq*t) + 1000.0
    jy=Q*r*freq*math.sin(freq*t)
    jz=Q*r*freq*math.sin(freq*t)
    return {jx,jy,jz}
  else
    return {0.0, 0.0, 0.0}
  end
end
```

In the if statement `d<=r`, `d` is the distance
of the actually coordinates (x,y,z) from the source in the centre of the box
and `r` is the radius which we specified above.
This if statement guarantees that the point source
will be a sphere of radius `r`
(Remember, Ateles is working on voxel cubes defined by the octree).
Thus for all elements which are lying in the `shape` the source term is active.
If the refinement level is not high enough,
the elements might be larger then the source should be.
Hence, this can be controlled by a condition.

Finally the source term is defined.

```lua
-- define the source term
source_terms = {
  name = 'current_density',
  ncomponents = 3,
  current_density = {
    kind = 'lua_fun',
    fun = currentDensitySpaceTime,
    shape = {
      kind = 'canoND',
      object = {
        origin = {-0.1, -0.1, -0.1,},
        vec = {{0.3, 0.0, 0.0}, {0.0, 0.3, 0.0}, {0.0, 0.0, 0.3}},
        segments = {100, 100, 100}
      }
    }
  }
}
```

`ncomponents` specifies how many variables the source term deals with.
Actually, this means how many return values the source term function has.
In the maxwell case there are the three currents jx, jy,jz, so `ncomponents` 
has to be 3.
Next, in current_density (the `name` of the source) we define
which kind of `function` should be used.
As we want to control it from within the lua script, we need to make it a lua
function. Further, we give the name of the function `currentDenitySpaceTime`
which is already specified (see above).
The variable `shape` tells the solver in
which part of the computational domain the source term is active.
Similar to [tracking](tracking.html#tracking-shape),
this shape is defined as `kind = 'canoND'`.
In `object` all the properties of the shape are specified.
Hence, we give a starting point (`origin`) and 3 direction vectors (`vec`)
which define a cube inside the computational domain.
The origin should be added by a small offset to make sure it is loacted
inside an element and not on the border of the octree.
For each direction vector,
we specify an amount of subdevisions of these vector.
Make sure to specify enough segments
to take into account all the required elements in the shape.
