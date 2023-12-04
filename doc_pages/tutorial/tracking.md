title: Tracking

During the runtime of a simulation
sometimes the user is interested in some state variables
but does not want to write down the complete domain (i.e. a restart file)
information on disk.
Therefore the tracking module is implemented to give the opportunity
to evaluate only a certain state variable.

[TOC]

## General Structure

To track a specific variable you need to modify the Lua input file and define

-   the variable `label` and the `folder` where the date should be stored,
-   the `shape`,
-   the `interval` for gathering the values and
-   the `output` with at least the `format` defined in it.

These important variables will be explained further
in the following paragraphs.

## Tracking Variable

Possible variables available for tracking are:

-   All conservative variables defined in the global variable system,
    depending on the equation system to solve.
-   Variables derived from the conservative variables. Can be predefined by
    the equation system, but the user can also define derived variables in
    the configuration's variable table
    (see [Variable System](|page|/treelm/page/features/variables/index.html)).
-   Eductions like difference, sum or average of defined variables.
-   Space-time-functions as an analytical result, can be defined in the Lua file.

As an example for the Maxwell equation,
you can track `displacement_field` and `magnetic_field`,
For the Euler equation there is, besides others, `density`, `momentum` and
`energy`, which you can track.
And another example: For the heat equation (nernst-planck equation)
the variables `concentration` and `diffusiveFlux` are, among others, trackable.

## Tracking Interval

Similar to
[configuring the time settings](general.html#general-simulation-settings)
for the restart output,
you specify when to start (`min`)
and when to stop (`max`),
plus the `interval` for writing the tracking output.

## Tracking Format

Ateles supports five formats for writing the tracked data.

-   `ascii`: Transient data - simple gnuplot-style.
    A simple ascii file which creates a small header
    and then writes a column wise representation of all variables
    and points at each defined time step.
    This is not efficient and can only be used
    for simple point-tracking objects.
    Each process writes its own file.
-   `asciiSpatial`: Spatial Data -
    This format stores solution data along
    with coordinates in single data file for every timestep.
    Therefore, there would be multiple files
    depending on the total number of time-steps.
    It is ideal for tracking line and plane.
-   `harvester`:
    Harvester format writes all participating elements to a TreElm mesh
    and creates a restart file for each required time step containing only the
    tracking data.
    This is a scalable tracking object
    which works regardless of the number of processes.
-   `vtk`:
    [Visualization Toolkit](https://en.wikipedia.org/wiki/VTK)
-   `precice`:
    [Precise Code Interaction Coupling Environment](http://www.precice.org)

## Tracking Shape

The tracking is realized using a simple shape.
The shape is defined as a canonical object.
The different shapes are basically derivatives of a simple cube.
A cube with no length in any direction is a point,
a cube with length along only one axis is a line
whereas a cube with length along two axes
is simply treated as plane by the underlying code.
Thus, the following shapes are implemented:

-   zero-dimensional: point
-   one-dimensional: line
-   two-dimensional: plane

For all these shapes we set the parameter `kind`
to the predefined name `canoND`
and configure the different shapes with the parameter `object`.
The following passages contain examples for each shape.
You can find these examples commented out in the configuration file
of the [Maxwell testcase](maxwell.html),
where we track the
[Electric displacement field](https://en.wikipedia.org/wiki/Electric_displacement_field)
(`variable = {'displacement_field'}`).
The result of the tracking will be stored in a folder called `tracking`
(`folder = 'tracking/'`).
So make sure to create this folder before running the simulation with tracking.

### Point Tracking

The most common thing is
to just track the values of a defined point in the computational domain.

```lua
tracking = {
  label = 'point_probe',
  folder = 'tracking/',
  variable = {'displacement_field'},
  shape = {
    kind = 'canoND',
    object = {
      origin = {0., 0., 0.}
    }
  },
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max / 2.0,
  },
  output = {
    format = 'ascii',
    use_get_point = true
  }
}
```

This part in the ateles.lua input file specifies the point tracking,
since only one point in `origin` is definded.
Ateles now generates an ascii file in the defined folder (`tracking`)
where you can find the results.
With the `use_get_point` option set to `true`
we will get the exact point values at the specified point.

### Line Tracking

For line tracking we need to define a starting point (`origin`)
and one direction vector (`vec`) in the parameter `object`.
The line is represented as a list of points in the tracking facility.
The number (`segments`)
and the `distribution` of the points on the line needs to be set.
Make sure to specify enough segments
to take into account all the required elements on the line.
Elements are added only once,
so it might be beneficial to specify much more segments than actually needed
to consider elements on the highest levels as well.

```lua
tracking = {
  label = 'line_probe',
  folder = 'tracking/',
  variable = {'displacement_field'},
  shape = {
    kind = 'canoND',
    object = {
      origin = {-1.0, 0.0, 0.0},
      vec = {{2.0, 0.0, 0.0}},
      segments = {4},
      distribution = 'equal'
    }
  },
  time_control = {
    min = 0.0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max / 2.0
  },
  output = {
    format = 'asciiSpatial',
    ndofs = 1
  }
}
```

After, the above tracking is executed,
the output is stored inside the defined folder (in this case: `tracking`).
The output gives the results for the elements lying on the defined line.
In the .res file only the integral mean value of the element is stored
(`ndofs = 1`).
Actually for line tracking it is convenient to use the harvester format,
since all degrees of freedom are stored
and you can produce visualization files.

### Plane Tracking

For a plane,
we need to define an origin and the two direction vectors.
The plane is, just as the line,
represented as a list of points in the tracking facility.
The number of subdivisions in each direction
and its distribution of the points on the plane needs to be set.
Make sure to specify enough segments
to take into account all the required elements on the plane.
Elements are added only once,
so it might be beneficial to specify much more segments
than actually needed to consider elements on the highest levels as well.

```lua
tracking = {
  label = 'line_probe',
  folder = 'tracking/',
  variable = {'displacement_field'},
  shape = {
    kind = 'canoND',
    object = {
      origin = {0., 0., 0.},
      vec = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}},
      segments = {500, 300},
      distribution = 'equal'
    }
  },
  time_control = {
    min = 0.0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max / 2.0
  },
  output = {
    format = 'harvester'
  }
}
```

Also for plane tracking, the harvester format is the most convenient one.
Similar to the line tracking,
an ascii-file would only give you the mean integral values for the element
which are lying in this plane.
In the harvester format all degrees of freedom are stored
and can be used in the visualization.

### Global Mesh Tracking

The global mesh tracking is defined
by assinging the predefined name `global` to the parameter `shape`.
This is an easy way to track only one variable in the whole domain
without writing a bunch of restart files.

```lua
tracking = {
  label = 'global_probe',
  folder = 'tracking/',
  variable = {'displacement_field'},
  shape = {
    kind = 'global'
  },
  time_control = {
    min = {
      iter = 0
    },
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max / 2.0
  },
  output = {
    format = 'vtk'
  }
}
```

Here we choose the Visualization Toolkit format (`format = 'vtk'`),
so that we can directly visualize the results with
[Paraview](http://www.paraview.org).
