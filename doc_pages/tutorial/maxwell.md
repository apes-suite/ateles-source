title: Maxwell Testcase

This tutorial provides a step by step introduction
to setting up a simulation of
[Maxwell equations](https://en.wikipedia.org/wiki/Maxwell%27s_equations)
in a homogenus media.
For the simulation domain we will use a periodic cube.

We assume that you already have read the
[General Procedure for Running the Testcases](procedure.html)
and the
[General Input File Information](general.html).

[TOC]

## The Input Lua File

In this tutorial, we will show and examine the input file for the Maxwell
testcase in `tutorial/testcase/maxwell/ateles.lua`.
Because this input file is specifically designed for the Maxwell testcase,
it does not make use of all the available options.
A complete input file with all available settings
can be found in the ateles folder `ateles/ateles.lua`.

## General Settings for the Simulation

First we give all the general settings for the simulation.
Since we are going to simulate a periodic oscillator using modal
[DG](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method)
scheme,
we set the name of the simulation to `posci_modg`.
We will simulate \(\sqrt{2}\) seconds
and every 10th iteration we get information about the current timestep
printed to the screen during the simulation.

```lua
-- global simulation options
simulation_name = 'posci_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    interval = {iter = 10},
    min = 0.0,
    max = math.sqrt(2),        -- final simulation time
  }
}
```

## Geometry

Further, the solver needs a defined geometry.
In this testcase we use a pre-defined geometry (`cube`),
which is just a simple cube with periodic boundaries and no obstacles.
We set the `length` of this cube to be 2 meter
and choose a `refinementLevel` of 2
(the space is cut into 4 (\(= 2^2\)) portions in each dimension). We position
the cube's center to the point x=0, y=0, z=0 by moving it half of it's width to
x-, y- and z-.

```lua
-- mesh definitions
cubeLength = 2.0
mesh = {
  predefined = 'cube',
  origin = {
    (-1.0) * cubeLength / 2.0,
    (-1.0) * cubeLength / 2.0,
    (-1.0) * cubeLength / 2.0
  },
  length = cubeLength,
  refinementLevel = 2
}
```

## Specific Parameters for this Testcase

It is possible to define global parameters in Lua,
which is a convenient way if some parameters are used several times
or need to be changed frequently for e.g. a parameter study.
Hence, we define some parameters, which will be used below in the function set.

```lua
-- some global parameters for the T_{nm} mode testcase
amplX = 1.0  -- the integer number of the mode in x direction
amplY = 1.0  -- the integer number of the mode in y direction
```

## Additional Functions for this Testcase

Lua also supports defining functions in the global scope.
For this testcase,
there exists an analytical solution if \(\epsilon = \mu = 1\). Hence,
we define functions which are time and space depended
for the temporal angular frequency
which is then used in the functions for the electric
and magnetic field in each direction.
Moreover, we specify a function for the initial condition of the electric field
in the z direction.

```lua
-- the analytic solution for this testcase is given by the following functions
-- (only correct for epsi = mu = 1):
-- ... definition of temporal angular frequency
w = math.sqrt(amplX^2 + amplY^2)
-- ... E_x = 0.0
function electricX(x, y, z, t)
  return 0.0
end
-- ... E_y = 0.0
function electricY(x, y, z, t)
  return 0.0  -- math.sin(amplX * math.pi * x) * math.sin(amplY * math.pi * z) * math.cos(w * t)
end
-- ... E_z(x, y, z, t) = sin(amplX \pi x) sin(amplY \pi y) cos(w t)
function electricZ(x, y, z, t)
  return math.sin(amplX * math.pi * x) * math.sin(amplY * math.pi * y) * math.cos(w * t)
end
-- ... B_x(x, y, z, t) = -\frac{\pi n}{w} sin(m \pi x) cos(n \pi y) sin(w t)
function magneticX(x, y, z, t)
  return (-1.0) * (math.pi * amplY/w) * math.sin(amplX * math.pi * x) * math.cos(amplY * math.pi * y) * math.sin(w * t)
end
-- ... B_y(x, y, z, t) = \frac{\pi m}{w} cos(m \pi x) sin(n \pi y) sin(w t)
function magneticY(x, y, z, t)
   return (math.pi * amplX / w) * math.cos(amplX * math.pi * x) * math.sin(amplY * math.pi * y) * math.sin(w * t)
end
-- ... B_z = 0.0
function magneticZ(x, y, z, t)
  return 0.0  -- (math.pi * amplX / w) * math.cos(amplX * math.pi * x) * math.sin(amplY * math.pi * z) * math.sin(w * t )
end
```

## Initial Condition

Next, the solver requires an initial condition.
Since we are considering the Maxwell equations,
the electric and magnetic field in each direction has to be defined.

```lua
-- initial condition
-- ...initial condition function for electric field (z component)
function ic_electricZ(x, y, z)
  return electricZ(x, y, z, 0.0)
end
initial_condition = {
  displacement_fieldX = 0.0,           -- electric field , x component
  displacement_fieldY = 0.0,           -- electric field , z component
  displacement_fieldZ = ic_electricZ,  -- electric field , z component
  magnetic_fieldX = 0.0,           -- magnetic induction , x component
  magnetic_fieldY = 0.0,           -- magnetic induction , y component
  magnetic_fieldZ = 0.0,           -- magnetic induction , z component
}
```

## Scheme

Of course, we need to define the parameter
for the scheme using the variable `scheme`.
The scheme is further split into the spatial and temporal discretization, each
described in a correspondigly named variable.

```lua
-- scheme definitions
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',        -- we use the modal discontinuous Galerkin scheme
    m = 4,                -- the maximal polynomial degree for each spatial direction
    modg_space = 'Q',
  },

  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',       -- the name of the timestep control mechanism
      cfl = 0.095,        -- Courant-Friedrichs-Levy number
    }
  }
}
```

For the spatial discretization scheme,
the parameter `name` specifies which scheme is used
whereby `modg` is the predfined name for modal discontinious Galerkin scheme.
Moreover, the parameter `m` defines the maximal polynomial degree
for each spatial direction, while `modg_space` defines
which polynomial space (P or Q) is used for the approximation.

Similary, for the temporal discretization scheme
the parameter `name` defines which scheme should be used.
Here we use the predefined name `explicitRungeKutta`.
To perform the classical four-step
[Runge-Kutta scheme](https://en.wikipedia.org/wiki/Runge–Kutta_methods#The_Runge.E2.80.93Kutta_method),
the parameter `steps` is set to `4`.
Since it is an explicit time stepping scheme,
each timestep needs to be controlled.
Therefore, we set the `control` variables as above.
In Ateles this is done by the Courant-Friedrichs-Lewy- or CFL-condition.
The parameter `cfl` defines the cfl number
which depends on the spatial discretiztion.

## Defining the Equations

Here, we define which system of equations should be solved.

```lua
-- equation definitions
equation = {
  name = 'maxwell',   -- we solve maxwell's equations
  material = {
    permeability = 'global_maxwell_permeability',
    permittivity = 'global_maxwell_permittivity',
    conductivity = 'global_maxwell_conductivity'
  }
}
```

Since this is the maxwell testcase,
we define `name= 'maxwell'` which is a predefined name in the Ateles solver.
Additionally, for the maxwell equation system,
the magnetic permeability, electric permitivity
and conductivity need to be set.
This is done with a material table.
In this table we are using variables which we define in the next section.

## Global Material
We define the global material for Maxwell.
It consists of three different components:

-   permeability
-   permittivity and
-   conductivity.

Each is a scalar, so we need three sclars for this equation system.
As this is the global fallback material,
we define each material to be a neutral term,
which in this case is 0.

```lua
-- Maxwell material parameters
mu = 1.0    -- the magnetic permeability (vacuum has 4.0 * math.pi * (10.0^(-7.0)))
epsi = 1.0  -- the electric permitivity (vacuum has 8.85418781762 * (10.0^(-12.0)))
conductivity = 0.0
variable = {
  -- This is the global material for Maxwell.
  -- It consists of three different components,
  -- permeability, permittivity, and conductivity, each a scalar,
  -- so that we need three scalars for this equation system.
  -- As this is the global fallback material, we define each material to be a neutral term,
  -- which in this case is 0.
  {
    name = "global_maxwell_permeability",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = mu }
  },
  {
    name = "global_maxwell_permittivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = epsi }
  },
  {
    name = "global_maxwell_conductivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = conductivity }
  }
}
```

## Restart Settings

Finally, we set the parameters for writing restart files.

```lua
-- configuration for the restart file
restart = {
  ---- file to restart from
  --read = './restart/posci_modg_lastHeader.lua',

  -- folder to write restart data to
  write = 'restart/',
  -- temporal definition of restart write
  time_control = {
    min = 0.0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max/5
  }
}
```

## Running the simulation

When you run the simulation, make sure that you are calling the solver Ateles
from the folder where the configuration file resides. This is
`tutorial/testcase/maxwell/ateles.lua`. When calling it from another folder, 
you also need to create the restart folder somewhere else, which will cause
some problems later on. First of all you have to take more care when calling
the postprocessing tool in the next step, but it is also easier to recap in
future, when the simulation input and output are stored in the same folder.
