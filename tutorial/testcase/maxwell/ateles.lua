-- Configuration File for Ateles (Periodic Oscillator)

-- This is a configuration file for the
-- Finite Volume / Discontinuous Galerkin Solver ATELES.
-- It provides a testcase for the simulation of
-- Maxwell equations in a homogenous media.
-- The simulation domain is a periodic cube with edge length 2.0.
-- Therefore this is a very good way to verify your algorithmic
-- implementations, since this testcase does not involve
-- any boundary conditions.
-- The testcase simulates the temporal development of standing waves
-- for electric and magnetic fields.
-- Since we are considering a very simple domain,
-- an analytic solution is well known and given
-- as Lua functions in this script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- To verify diffusion and dispersion relations this testcases allows you
-- to modify the spatial harmonics
-- by varying the integer mode number in x and y direction
-- by modifying the lua variables m and n.
-- Please notice,
-- this testcase is correct only for homogenous media with epsi = mu = 1
-- (see equations table).
-- This testcase can be run in serial (only one execution unit)
-- or in parallel (with multiple mpi ranks).
-- To calculate a grid convergence behavior please modify the level variable.
-- An increment of one will half the radius of your elements.


--! [General Settings for the Simulation]
-- global simulation options
simulation_name = 'posci_modg' -- the name of the simualtion
sim_control = {
  time_control = {
    interval = {
      iter = 10
    },
    min = 0.0,
    max = math.sqrt(2),  -- final simulation time
  }
}
--! [General Settings for the Simulation]


--! [Geometry]
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
--! [Geometry]


--! [Specific Parameters for this Testcase]
-- Some global parameters for the T_{nm} mode testcase
amplX = 1.0  -- the integer number of the mode in x direction
amplY = 1.0  -- the integer number of the mode in y direction
--! [Specific Parameters for this Testcase]


--! [Additional Functions for this Testcase]
-- The analytic solution for this testcase is given by the following functions
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
--! [Additional Functions for this Testcase]


--! [Initial Condition]
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
--! [Initial Condition]


--! [Scheme]
-- scheme definitions
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',  -- we use the modal discontinuous Galerkin scheme
    m =  4,         -- the maximal polynomial degree for each spatial direction
    modg_space = 'Q',
  },

  -- the temporal discretization scheme
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    -- how to control the timestep
    control = {
      name = 'cfl',   -- the name of the timestep control mechanism
      cfl  = 0.095,     -- Courant-Friedrichs-Levy number
    },
  },
}
--! [Scheme]


--! [Defining the Equations]
-- equation definitions
equation = {
  name = 'maxwell',  -- we solve maxwell's equations
  material = {
    permeability = 'global_maxwell_permeability',
    permittivity = 'global_maxwell_permittivity',
    conductivity = 'global_maxwell_conductivity'
  }
}
--! [Defining the Equations]


--! [Global Material]
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
    st_fun = {
      const = mu
    }
  },
  {
    name = "global_maxwell_permittivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = epsi
    }
  },
  {
    name = "global_maxwell_conductivity",
    ncomponents = 1,
    vartype = "st_fun",
    st_fun = {
      const = conductivity
    }
  }
}
--! [Global Material]


--! [Restart Settings]
-- configuration for the restart file
restart = {
  ---- file to restart from
  --read = './restart/posci_modg_lastHeader.lua',

  ---- folder to write restart data to
  write = 'restart/',
  -- temporal definition of restart write
  time_control = {
    min = 0.0,
    max = sim_control.time_control.max,
    interval = sim_control.time_control.max / 5
  }
}
--! [restart]


--! [point tracking]
--tracking = {
--  label = 'point_probe',
--  folder = 'tracking/',
--  variable = {'displacement_field'},
--  shape = {
--    kind = 'canoND',
--    object = {
--      origin = {0., 0., 0.}
--    }
--  },
--  time_control = {
--    min = 0,
--    max = sim_control.time_control.max,
--    interval = sim_control.time_control.max / 2.0
--  },
--  output = {
--    format = 'ascii',
--    use_get_point = true
--  }
--}
--! [point tracking]

--! [line tracking]
--tracking = {
--  label = 'line_probe',
--  folder = 'tracking/',
--  variable = {'displacement_field'},
--  shape = {
--    kind = 'canoND',
--    object = {
--      origin = {-1.0, 0.0, 0.0},
--      vec = {{2.0, 0.0, 0.0}},
--      segments = {4},
--      distribution = 'equal'
--    }
--  },
--  time_control = {
--    min = 0.0,
--    max = sim_control.time_control.max,
--    interval = sim_control.time_control.max / 2.0
--  },
--  output = {
--    format = 'asciiSpatial',
--    ndofs = 1
--  }
--}
--! line tracking

--! [plane tracking]
--tracking = {
--  label = 'line_probe',
--  folder = 'tracking/',
--    variable = {'displacement_field'},
--    shape = {
--      kind = 'canoND',
--      object = {
--        origin = {0., 0., 0.},
--        vec = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}},
--        segments = {500, 300},
--        distribution = 'equal',
--      }
--    },
--    time_control = {
--      min = 0.0,
--      max = sim_control.time_control.max,
--      interval = sim_control.time_control.max / 2.0
--    },
--    output = {format = 'harvester'}
--}
--! [plane tracking]

--! [global tracking]
--tracking = {
--  label = 'global_probe',
--  folder = 'tracking/',
--  variable = {'displacement_field'},
--  shape = {
--    kind = 'global'
--  },
--  time_control = {
--    min = {
--      iter = 0
--    },
--    max = sim_control.time_control.max,
--    interval = sim_control.time_control.max / 2.0
--  },
--  output = {
--    format = 'vtk'
--  }
--}
--! [global tracking]


--! [Source term]
---- ... charge of the source
--Q = 1.0
---- ... radius of sphere source
--r = 0.4
---- ... parameters for the analytic solution
--freq = ( 2.0*math.pi/math.sqrt(mu*epsi) )
---- ... the temporal period of the waveguide
--T = 2.0*math.pi/freq
---- function for the point source
--function currentDensitySpaceTime(x, y, z, t)
--  d = math.sqrt(x^2.0+y^2.0+z^2.0)
--  if d <= r then
--    jx=Q*r*freq*math.sin(freq*t) + 1000.0
--    jy=Q*r*freq*math.sin(freq*t)
--    jz=Q*r*freq*math.sin(freq*t)
--    return {jx,jy,jz}
--  else
--    return {0.0, 0.0, 0.0}
--  end
--end
---- define the source term
--source_terms = {
--  current_density = {
--    fun = currentDensitySpaceTime,
--    shape = {
--      kind = 'canoND',
--      object = {
--        origin = {-0.1, -0.1, -0.1,},
--        vec = {{0.3, 0.0, 0.0}, {0.0, 0.3, 0.0}, {0.0, 0.0, 0.3}},
--        segments = {100, 100, 100}
--      }
--    }
--  }
--}
--! [Source term]

----! [posci_setup_exact]
---- The analytic solution for this testcase is given by the following functions
---- (only correct for epsi = mu = 1):
---- ... definition of temporal angular frequency
--w = math.sqrt(amplX^2 + amplY^2)
---- ... E_x = 0.0
--function electricX(x, y, z, t)
--  return 0.0
--end
---- ... E_y = 0.0
--function electricY(x, y, z, t)
--  return 0.0  --math.sin(amplX * math.pi * x) * math.sin(amplY * math.pi * z) * math.cos(w * t)
--end
---- ... E_z(x, y, z, t) = sin(amplX \pi x) sin(amplY \pi y) cos(w t)
--function electricZ(x, y, z, t)
--  return math.sin(amplX * math.pi * x) * math.sin(amplY * math.pi * y) * math.cos(w * t)
--end
---- ... B_x(x, y, z, t) = -\frac{\pi n}{w} sin(m \pi x) cos(n \pi y) sin(w t)
--function magneticX(x, y, z, t)
--  return (-1.0) * (math.pi * amplY / w) * math.sin(amplX * math.pi * x) * math.cos(amplY * math.pi * y) * math.sin(w * t)
--end
---- ... B_y(x, y, z, t) = \frac{\pi m}{w} cos(m \pi x) sin(n \pi y) sin(w t)
--function magneticY(x, y, z, t)
--   return (math.pi * amplX / w) * math.cos(amplX * math.pi * x) * math.sin(amplY * math.pi * y) * math.sin(w * t)
--end
---- ... B_z = 0.0
--function magneticZ(x, y, z, t)
--  return 0.0  --(math.pi * amplX / w) * math.cos(amplX * math.pi * x) * math.sin(amplY * math.pi * z) * math.sin(w * t)
--end
----! [posci_setup_exact]
--
--
------ PROJECTION
------ In the projection table, the method for projection from nodal to modal space for BC,
------ IC and non-;inear problems is given
----projection = {
----  -- kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
----  -- for fpt the  nodes are automatically 'chebyshev'
----  -- for lep the  nodes are automatically 'gauss-legendre'
----  -- lobattoPoints = false  -- if lobatto points should be used, default = false,
----                         -- only working for Chebyshev points --> fpt
----  factor = 1.0,          -- dealising factpr for fpt
----                         -- oversampling factor to remove aliasing
----                         -- effects by padding, default: 1
----                         -- Note, that for FPT we require a power
----                         -- of 2 in the oversampled polynomial
----                         -- order, if no power of two is reached
----                         -- the padding is automatically increased
----                         -- accordingly.
----  -- approx_terms = 18,     -- Number of terms used to approximate the
----  --                        -- matrix multiplication for blocks, that
----  --                        -- are detached from the diagonal.
----  --                        -- The default of 18 is recommended for
----  --                        -- double precision.
----  -- blocksize = 64,        -- for FPT, default 64. The blocksize
----  --                        -- defines how big the minimal block
----  --                        -- should be that is approximated in
----  --                        -- fast algorithm.
----  --                        -- The smaller it is, the more operations
----  --                        -- are merely approximated.
----  --                        -- Recommended for double precision is a
----  --                        -- setting of 64.
----  --                        -- The fast algorithm will only be used
----  --                        -- for m >= blocksize.
----  --                        -- Note, that this has to be larger than
----  --                        -- 2*approx_terms to provide any
----  --                        -- reduction in operation counts.
----  -- striplen = 256,        -- This provides the length for arrays to
----  --                        -- apply the matrix operation to
----  --                        -- simultaneously.
----  --                        -- Default is the vlen from the tem_compileconf_module.
----  -- subblockingWidth = 8   -- The subblockingWidth is used during the
----  --                        -- unrolling of the diagonal multiplication
----  --                        -- during the projection. By setting this
----  --                        -- value to an appropriate value a better
----  --                        -- cache usage can be achieved.
----  --                        -- Default is 8
----  -- adapt_factor_pow2 = true -- for FPT, default false. Should the
----  --                          -- oversampling factor be adjusted to
----  --                          -- obtain a power of 2 in the
----  --                          -- oversampled order?
----  --                          -- multithreading be activated for FFTW?
----  --
----  -- individual projection methods for source terms
----  -- if non is determined the general projection method is used
----  source_terms = {
----    -- the configuration parameter are similar to the
----    -- general projection method
----    -- kind = 'fpt',
----    factor = 1.0
----  },
----  -- individual projection methods for initial condition
----  -- if non is determined the general projection method is used
----  initial_condition = {
----    -- the configuration parameter are similar to the
----    -- general projection method
----    -- kind = 'fpt',
----    factor = 2.0
----  },
----  -- individual projection methods for boundary condition
----  -- if non is determined the general projection method is used
----  boundary_condition = {
----    -- the configuration parameter are similar to the
----    -- general projection method
----    -- kind = 'fpt',
----    factor = 2.0
----  }
----}
