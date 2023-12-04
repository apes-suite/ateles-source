-- Configuration File for Ateles

-- This is a configuration file for the Finite Volume / Discontinuous Galerkin Solver ATELES.
-- It provides a testcase for the simulation of Euler equations in a homogenous media.
-- The simulation domain is a periodic cube with edge length 2.0.
-- Therefore this is a very good way to verify your algorithmic implementations,
-- since this testcase does not involve any boundary conditions.
-- The testcase simulates the temporal development of Gaussian pulse in density.
-- Since we are considering a very simple domain,
-- an analytic solution is well known and given as Lua functions in this script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- This testcase can be run in serial (only one execution unit)
-- or in parallel (with multiple mpi ranks).
-- To specify the number of ranks please modify nprocs variable.
-- To calculate a grid convergence behavior please modify the level variable.
-- An increment of one will half the radius of your elements.


--! [General Simulation Settings]
-- the length of the cube
cubeLength = 2.0

-- the refinement level of the octree
level = 2

-- transport velocity of the pulse in x direction.
velocityX = 100

-- global simulation options
simulation_name = 'gPulseDens_euler_modg'   -- the name of the simualtion
sim_control = {
  time_control = {
    interval = {
      iter = 10
    },
    min = 0,
    max = cubeLength / velocityX / 4,       -- final simulation time
  }
}
--! [General Simulation Settings]


--! [Geometry]
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
--! [Geometry]


--! [Initial Conditions]
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
--! [Initial Conditions]


--! [Scheme]
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
--! [Scheme]


--! [Defining the Equations]
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
--! [Defining the Equations]


--! [Global Material]
-- material variables
variable = {
  {
    name = "global_characteristic",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = 0.0
    }
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
    st_fun = {
      const = 0.0
    }
  }
}
--! [Global Material]


--! [Restart Settings]
-- configuration for the restart file
restart = {
  ---- file to restart from
  --read = './restart/gPulseDens_euler_modg_lastHeader.lua',

  -- folder to write restart data to
  write = './restart/',
  -- temporal definition of restart write
  time_control = {
    interval = {
      iter = 10
    },
    min = 0,
    max = sim_control.time_control.max,
  }
}
--! [Restart Settings]


--projection = {
--  kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
--                 -- for fpt the  nodes are automatically 'chebyshev'
--                 -- for lep the  nodes are automatically 'gauss-legendre'
--                 -- lobattoPoints = false  -- if lobatto points should be used, default = false,
--                 -- only working for Chebyshev points --> fpt
--  factor = 1.0,  -- dealising factpr for fpt
--                 -- oversampling factor to remove aliasing
--                 -- effects by padding, default: 1
--                 -- Note, that for FPT we require a power
--                 -- of 2 in the oversampled polynomial
--                 -- order, if no power of two is reached
--                 -- the padding is automatically increased
--                 -- accordingly.

--  -- blocksize = 32,        -- for fpt, default -1
--  -- fftMultiThread = false -- for fpt, logical, default false
--}
