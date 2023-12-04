title: Euler Testcase

This section describes running the testcase using the
[Euler equations](https://en.wikipedia.org/wiki/Euler_equations_(fluid_dynamics)).
Here we would take a small cubic geometry and define a gauss-pulse in density,
i.e we set the initial density such that
it has the peak in the center of the cube
and it decreases towards the sides making a
[gaussian-pulse](http://en.wikipedia.org/wiki/Gaussian_function).
Then we would give some velocity in one spatial direction
(here, we choose x- direction),
and watch this pulse in density move through the cube.

We assume that you already have read the
[General Procedure for Running the Testcases](procedure.html)
and the
[General Input File Information](general.html).

[TOC]

## The Input Lua File

In this tutorial, we will show
and examine the input file for the Euler testcase.
Because this input file is specifically designed for the Euler testcase,
it does not make use of all the available options.
A complete input file with all available settings
can be found in the ateles folder `ateles/ateles.lua`.

## General Simulation Settings

For clarity we first define some variables, which we will use later.
We will simulate a gaussian pulse in density using modal
[DG](https://en.wikipedia.org/wiki/Discontinuous_Galerkin_method) scheme,
hence here we would just name it as `gPulseDens_euler_modg`.
In this example we are going to simulate \(\frac{1}{200}\) seconds
(`max = cubeLength/velocityX/4`).
Every 10th iteration information
about the current timestep is printed to the screen.

```lua
-- the length of the cube
cubeLength = 2.0

-- the refinement level of the octree
level = 1

-- transport velocity of the pulse in x direction.
velocityX = 100

-- global simulation options
simulation_name = 'gPulseDens_euler_modg'   -- the name of the simualtion
sim_control = {
  time_control = {
    interval = { iter = 10 },
    min = 0,
    max = cubeLength / velocityX / 4,       -- final simulation time
  }
}
```

## Geometry

For the geometry we choose a predefined cube
with periodic bounderaies and no obstacles.

```lua
-- mesh definitions
mesh = {
  predefined = 'cube',
  origin = {
    (-1.0) * cubeLength / 2.0,
    (-1.0) * cubeLength / 2.0,
    (-1.0) * cubeLength / 2.0
  },
  length = cubeLength,
  refinementLevel = level
}
```

## Initial Conditions

Here for our case we will specify an initial condition
for the gauss-pulse in density (i.e set the initial density
such that it has a peak in the middle of the domain,
and decreases quickly towards the sides)
and make it travel in the x-direction.
This could also be a check for your boundary conditions
or just to visualize your numerical error
that smears out this gauss-pulse in time.
So, we set the `velocityY` and `velcoityZ` to zero
and give some value to the `velcoityX` and `pressure`.
For `density` we have a set of predefined variables
providing necessary information for the gauss pulse.

```lua
-- initial conditions
initial_condition = {
  density = {
    predefined = 'gausspulse',
    center = { 0.0, 0.0, 0.0 },  -- gauss pulse center
    halfwidth = 0.20,            -- half width of gauss pulse from center
    amplitude = 2.0,             -- height or magnitude of gauss pulse
    background = 1.225           -- reference value. In case of density, it is reference density.
  },
  pressure = 100000,
  velocityX = velocityX,
  velocityY = 0.0,
  velocityZ = 0.0,
}
```

## Scheme

For the variable `spatial` we use the predefined name `modg`
which represents the modal Discontinuous Galerkin Scheme.
Similarly, for the time stepping scheme we have the variable `temporal`.
We use the predefined `explicitRungeKutta`
and are able to perform multi-step
[Runge-Kutta](https://en.wikipedia.org/wiki/Runge–Kutta_methods)
by setting the `steps`.
To control the timesteps we use the
[CFL number](https://en.wikipedia.org/wiki/Courant–Friedrichs–Lewy_condition)
which can be specified in the `control` variable.

```lua
-- scheme definitions
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',            -- we use the modal discontinuous Galerkin scheme
    m = 6,                    -- the maximal polynomial degree for each spatial direction
    dealiasFactor = 1.0,      -- factor to remove aliases: 1.0 means no additional dealiasing
    blocksize = 32,           -- the minimal blocksize for the FPT
    fftMultiThread = true     -- use multithreaded version of FFTW
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',  -- 'explicitEuler'
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',   -- the name of the timestep control mechanism
      cfl  = 0.8      -- Courant_Friedrichs_Lewy number
    }
  }
}
```

## Defining the Equations

The equations that we solve are the Euler equations
and it can be given through the input `name` in the `equation` table.

```lua
-- equation definitions
equation = {
  name   = 'euler',        -- we solve euler equations
  therm_cond = 2.555e-02,  -- thermal conductivity (default = 6.953195-310)
  isen_coef = 1.4,         -- isentropic coefficient
  r      = 296.0,          -- ideal gas constant
  material = {
    characteristic = 'global_characteristic',
    relax_velocity = 'global_relax_velocity',
    relax_temperature = 'global_relax_temperature'
  }
}
-- (cv) heat capacity
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)
```

## Global Material
Here we define the global material for Euler 3D.
It consists of three different components:

-   characteristics with one scalar,
-   relax\_velocity with three scalars and
-   relax\_temperature with one scalar.

Thus we need five scalars for this equation system.
As this is the global fallback material,
we define each material to be a neutral term,
which in this case is 0.

```lua
-- material variables
variable = {
  {
    name = "global_characteristic",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = 0.0 }
  },
  {
    name = "global_relax_velocity",
    ncomponents = 3,
    vartype = "st_fun",
    st_fun = {
      const = { 0.0, 0.0, 0.0 }
    }
  },
  {
    name = "global_relax_temperature",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = { const = 0.0 }
  }
}
```

## Restart Settings

We want dumping restarts every 10th itertion till the end of the simulation.
Because we are not restarting from a file the variable `read` is commented out.

```lua
-- configuration for the restart file
restart = {
  ---- file to restart from
  --read = './restart/gPulseDens_euler_modg_lastHeader.lua',

  -- folder to write restart data to
  write = './restart/',
  -- temporal definition of restart write
  time_control = {
    interval = { iter = 10 },
    min = 0,
    max = sim_control.time_control.max,
  }
}
```
