title: General Input File Information

This section describes the general information on setting up
your testcases with Ateles.

[TOC]

## The Input Lua file

The input file written in [Lua](http://www.lua.org/) provides all simulation
parameter like geometry, initional condition, scheme, etc. for the Ateles
solver. Lua syntax is fairly easy to follow and we go on describing all the
features we use in detail.
Hint: "--" specifies a comment.
It is also good to know that Lua is case sensitive.
Additional, Lua executes the file in linear fashion,
thus everything you want to use (like functions for initial condition)
need to be defined before use.
The general overall structure of the Lua input file for Ateles
can be seen in the file `ateles.lua` present in the repository root
which we assume you already checked out.
This also is the default name for the input file,
however different names can be used in case you want to use
and keep multiple scripts.
You would just need to pass those input file as a command line argument.
If no input file is specified,
the code checks and runs for the default input file if present.

## General Simulation Settings

Every simulation should have a name,
which is specified by a variable called `simulation_name`.
If you would like for example to set up a testcase for flow through a cube,
it could be called 'flow\_cube'.

    simulation_name = 'flow_cube'

You can use just anything as a simulation name,
just make it clear and descriptive,
as it will help you to identify produced results from the simulation.
The simulation name will appear in all output files,
so if you use good and unique names,
you always know where a particular output came from.

The variable `sim_control` is used to define the stopping criteria or the time
when you would like the simulation to finish.
Currently there are three possible options:

1.  Running simulation for certain specific steps in simulation time:

        sim_control = {
          time_control = {
            min = 0,
            max = {sim = 5.0}
          }
        }

    Here `max` denotes the 'final simulation time'
    (amount of time we want to simulate in seconds).
    In this case it will end after 5 seconds.

1.  Running simulation for certain maximum iterations:

        sim_control = {
          time_control = {
            min = 0,
            max = {iter = 1000}
          }
        }

    Here `max` denotes the 'maximum number of iterations'
    after which the simulation stops.
    In this case the simulation will end after 1000 iterations.

1.  Running simulation for a specific amount of computation time:

        sim_control = {
          time_control = {
            min = 0,
            max = {clock = 3600}
          }
        }

    Here `max` denotes the 'maximum computation time' in seconds.
    In this case it is 1 hour (3600s).

Note that if `max` is not defined as a table, e.g. `max = 1.0`,
it refers to the simulation time.

You can also specify the frequency at which information about the current time
is printed to the screen.
For that you have to define a table called `interval` in `time_control`.
Like `max` (and also `min`) this is a time setting
and can be provided in terms of the prescribed definitions.

## Geometry

Complicated geometries are created by Seeder,
a preprocessing tool with which you can get familiar with by visiting its
[documentation](https://apes-suite.github.io/seeder/index.html).
For the moment, we will consider a simple pre-defined geometry
(in our case a cube with periodic boundries and no obstacles).

```lua
mesh = {
  predefined = 'cube',
  origin = {
    -1.0,
    -1.0,
    -1.0
  },
  length = 2,
  refinementLevel = 4
}
```

We can specify the `length` and position (`origin`) of that cube.
The position of the mesh is only relevant
if you specify positions of other objects, too,
such as initial conditions or [trackers](tracking.html).
Next, the `refinementLevel` tells Ateles how fine the space is discretized.
A refinement level of 4 means that the space is cut into 2^4 portions
in each dimension.

## Scheme

The scheme definitions can be given through the variable `scheme`.
It further contains the variables `spatial` and `temporal`.
It is possible to use predefined names such as 'modg'
(for modal discontinuous Galerking scheme) in the `spatial` variable.
The maximal polynomial degree for each spatial direction is given
as input through the variable `m`.
The `temporal` variable can also take predefined names as an input.

## Defining the Equations

Next we define the equation system and the related parameters.
The possible equations that are implemented in Ateles
can be used following the below mentioned syntax.

-   Euler Equations:

        equation = {
          name = 'euler',              -- we solve euler equations
          numflux = 'lax_friedrich',   -- numerical flux to use
                                       -- possible are 'lax_friedrich', 'godunov' or 'hll'
          isen_coef = 1.4,             -- isentropic coefficient
          r = 1.0,                     -- ideal gas constant
          material = {
            characteristic = 'global_euler3D_characteristic',
            relax_velocity = 'global_euler3D_relax_velocity',
            relax_temperature = 'global_euler3D_relax_temperature'
          }
        }
        -- cv heat capacity
        equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-   Navier Stokes Equations:

        equation = {
          name = 'navier_stokes',      -- we solve navier stokes equations
          mu = 1.831e-05,              -- viscosity  (default = 2.121996-314)
          therm_cond = 2.55e-02,       -- thermal conductivity (default = 6.953195-310)
          isen_coef = 1.4,             -- isentropic coefficient
          r = 296.0,                   -- ideal gas constant
          material = {
            characteristic = 'global_euler3D_characteristic',
            relax_velocity = 'global_euler3D_relax_velocity',
            relax_temperature = 'global_euler3D_relax_temperature'
          }
        }
        --(cv) heat capacity and (r) ideal gas constant
        equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-   Maxwells Equations:

        equation = {
          name = 'maxwell',-- we solve maxwell’s equations
          material = {
            permeability = 'global_maxwell_permeability',
            permittivity = 'global_maxwell_permittivity',
            conductivity = 'global_maxwell_conductivity'
          }
        }

-   Filtered Navier Stokes Equations:

        equation = {
          name = 'filtered_navier_stokes',  -- we solve filtered navier stokes equations
          mu = 1.831e-05,                   -- viscosity
          therm_cond = 2.55e-02,            -- thermal conductivity
          isen_coef = 1.4,                  -- isentropic coefficient
          r = 296.0,                        -- ideal gas constant
          nvar = 5,
          -- Turbulence parameter for Smagorinsky model
          turbulence_model = 'smagorinsky',
          -- Germano-Lilli model
          --turbulence_model = 'germano'
          Cs = 0.04,                        -- suggested value between 0.1 and 0.2
          Ci = 0.005,                       -- Value for Yoshizawa model
          prandtl_sgs = 0.6,                -- suggested value
          material = {
            characteristic = 'global_euler3D_characteristic',
            relax_velocity = 'global_euler3D_relax_velocity',
            relax_temperature = 'global_euler3D_relax_temperature'
          }
        }

-   Maxwell Equations with divergence correction:

        equation = {
          name = 'maxwellDivCorrection',    -- we solve maxwell’s equations with divergence correction
          mu = perma,                       -- the magnetic permeability of vacuum
          epsi = permit,                    -- the electric permitivity of vacuum
          chi = 1.0,                        -- parameter for electric divergence cleaning
          gam = 1.0,                        -- parameter for magnetic divergence cleaning
          material = {
            permeability = 'global_maxwell_permeability',
            permittivity = 'global_maxwell_permittivity',
            conductivity = 'global_maxwell_conductivity'
          }
        }

## Equation specific configuration

Some of the configuration is equation specific, as different equation systems
consist of different terms and variables and therefore need different
configuration settings. Thus we are not going into detail about each and every
setting here.

All of these settings can be defined by constant values, functions in the
Lua configuration file, or as predefined functions implemented in the solver.
It is important to know that constants defined in the configuration file are
stored in the 'solver environment' while initializing the simulation, whereas
Lua functions need to be evaluated at every timestep. Taking this into
consideration, constants or predefined values are handled within the solver,
but Lua functions imply the change of the 'evaluation environment', i.e. we
need to run Lua code, which is inherently slow.

### Initial Conditions

Initial conditions can be specified using the `initial_condition` variable.
Through this all the initial condition
for the state variables can be specified.
These might be predefined functions, Lua functions or just constants.

### Material

Material denotes parameters bound to specific reagions of the domain. With
material, you can e.g. specify a sphere inside your flow domain with a given
porosity. Material is also defined by either constants, Lua functions or
predefined functions.

### Source Terms

@todo Describe source terms in short words

## Restart Settings

Restart files are the main method for writing data during the simulation
and can also used for post-processing e.g. visualization.
For writing, the name of the restart folder
needs to be specified in the variable `write`.
The temporal definitions of the restart files
can be given through the variable `time`.
In the `time_control` table one can specify

-   the point in time to start dumping restarts with the variable `min`,
-   the point in time to stop dumping restarts with the variable `max`,
-   and the time between dumping two restart files
    within the timeframe given by `min` and `max` with the variable `interval`.

Again `min`, `max` and `interval` are time settings
and can be defined in the same way as described in
[General Simulation Settings](#general-simulation-settings).

As the name implies, restart files can also be used to restart a simulation
with the state stored in it. To read a restart file, the variable `read` is
used. It defines the file to read and restart the simulation from.

```lua
restart = {
  ---- file to restart from
  --read = './restart/simulation_header_1.892E-03.lua',

  -- folder to write restart data to
  write = './restart/',

  -- temporal definition of restart write
  time_control = {
    min = 0,
    max = {clock = 3600},
    interval = {iter = 10}
  }
}
```

In the above example the generation of the restart files
starts when the simulation begins
and ends at least after 3600 seconds (one hour) of computation time.
In this timeframe every 10th iteration a restart file is being created.
The variable `read` is commented out for this example.

@note
You should make sure that you create the directory
for storing the restart files in advance.
This needs to be done manually.
