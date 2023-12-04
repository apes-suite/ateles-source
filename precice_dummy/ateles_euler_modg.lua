-- Configuration file for Ateles --


-- This is a configuration file for the Finite Volume / Discontinuous Galerkin Solver ATELES. 
-- It provides a testcase for the simulation of Euler equations in a homogenous media. The simulation domain
-- is a periodic cube with edge length 2.0. Therefore this is a very good way to verify your algorithmic implementations, 
-- since this testcase does not involve any boundary conditions. 
-- The testcase simulates the temporal development of Gaussian pulse in density. Since we 
-- are considering a very simple domain, an analytic solution is well known and given as Lua functions in this script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- This testcase can be run in serial (only one execution unit) or in parallel (with multiple mpi ranks). 
-- To specify the number of ranks please modify nprocs variable. To calculate a grid convergence behavior please modify the 
-- level variable. An increment of one will half the radius of your elements.

timestep_info = 10

-- Transport velocity of the pulse in x direction.
velocityX = 100
cubeLength = 2.0
-- global simulation options
simulation_name = 'gPulseDens_euler_modg' -- the name of the simualtion
sim_control = {
                time_control = {   
                  min = 0, 
                  max = 7.7e-04 -- final simulation time
                }
              }

-- table for preCICE
precice = {
           accessor = 'Ateles',
           configFile ='config_pythonAction.xml',
           meshname = 'AcousticSurface',
           exchange_data =  {'Density', 'Velocities', 'Pressure'}
           --exchange_data = {'Pressure', 'Velocities'}
          }

-- Mesh definitions --
mesh = 'mesh/'
--mesh = { predefined = 'cube',
--         origin = { 
--                    (-1.0)*cubeLength/2.0,
--                    (-1.0)*cubeLength/2.0,
--                    (-1.0)*cubeLength/2.0
--                  },
--         length = cubeLength,
--         refinementLevel = level
--       }
--

---- Restart settings
estart = { 
--            -- file to restart from
--            read = './restart/maxwell/per_osci_maxwell_modg_lastHeader.lua',                        
--            -- folder to write restart data to
            write = './restart/precice/',                                        
            -- temporal definition of restart write
            time_control = {   
              min = 0, 
              max = sim_control.time_control.max, 
              interval = sim_control.time_control.max/30.0
            }  
          }

-- timing settings (i.e. output for performance measurements, this table is otional)
timing = {
          folder = './',                  -- the folder for the timing results
          filename = 'timing.res'         -- the filename of the timing results
         }

-- Equation definitions --
equation = {
    name   = 'euler',
    therm_cond = 2.555e-02,
    isen_coef = 1.4,
    r      = 296.0,
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
    -- the spatial discretization scheme
    spatial =  {
               name = 'modg',            -- we use the modal discontinuous Galerkin scheme 
               m = 2,                   -- the maximal polynomial degree for each spatial direction
               }, 
    -- the temporal discretization scheme
    temporal = {
               name = 'explicitRungeKutta',  --'explicitEuler',
               steps = 4,
               -- how to control the timestep
               control = {
                          name = 'cfl',   -- the name of the timestep control mechanism
                          cfl  = 0.8,     -- Courant–Friedrichs–Lewy number
                         },
               },
}

function dens(x,y,z)
  return x
end

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              -- for fpt the  nodes are automatically 'chebyshev'
              -- for lep the  nodes are automatically 'gauss-legendre'
           -- lobattoPoints = false  -- if lobatto points should be used, default = false
              factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
           -- blocksize = 32,        -- for fpt, default -1
           -- fftMultiThread = false -- for fpt, logical, default false
             }

-- This is a very simple example to define constant boundary condtions.
initial_condition = { density = 3.0, 
                             --  {
                             --    predefined='gausspulse',
                             --    center={1.8, 1.8, 1.8},
                             --    halfwidth=0.20,
                             --    amplitude=2.0,
                             --    background=1.225
                             --   },
                      -- pressure = {
                      --              predefined='gausspulse',
                      --              center={0.0, 0.0, 0.0},
                      --              halfwidth=0.2,
                      --              amplitude=2.0,
                      --              background=100000
                      --            },
                      pressure = 100000, 
                      velocityX = 0.0,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                    }


-- Tracking              
tracking = {
             label = 'track_momentum_precice_10',
             folder = './',
             variable = {'momentum'},
             shape = {kind = 'canoND', object= { origin ={2.5, 0., 0.} } },
             time_control = {
               min = 0,
               max = sim_control.time_control.max,
               interval = sim_control.time_control.max/10.0
             },
             format = 'ascii'
           }

 -- Boundary definitions
 boundary_condition = {
                         {
                           label = 'wall_1_out',
                           kind = 'outflow',
                           pressure = 100001,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_2_out',
                           kind = 'outflow',
                           pressure = 100001,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_3_out',
                           kind = 'outflow',
                           pressure = 100001,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_4_out',
                           kind = 'outflow',
                           pressure = 100001,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_5_out',
                           kind = 'outflow',
                           pressure = 100001,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_6_out',
                           kind = 'outflow',
                           pressure = 100001,
                           density = 'extrapolate'
                         }
                         ,
                         {
                           label = 'wall_1_in',
                           kind = 'precice',
                         --  pressure = '10000',
                         --  density = 1,
                         --  v_x = 0.0,
                         --  v_y = 0.0,
                         --  v_z = 0.0
                         }
                         ,
                         {
                           label = 'wall_2_in',
                           kind = 'precice',
                         --  pressure = '10000',
                         --  density = 1,
                         --  v_x = 0.0,
                         --  v_y = 0.0,
                         --  v_z = 0.0
                         }
                         ,
                         {
                           label = 'wall_3_in',
                           kind = 'precice',
                         --  pressure = '10000',
                         --  density = 1,
                         --  v_x = 0.0,
                         --  v_y = 0.0,
                         --  v_z = 0.0
                         }
                         ,
                         {
                           label = 'wall_4_in',
                           kind = 'precice',
                         --  pressure = '10000',
                         --  density = 1,
                         --  v_x = 0.0,
                         --  v_y = 0.0,
                         --  v_z = 0.0
                         }
                         ,
                         {
                           label = 'wall_5_in',
                           kind = 'precice',
                         --  pressure = '10000',
                         --  density = 1,
                         --  v_x = 0.0,
                         --  v_y = 0.0,
                         --  v_z = 0.0
                         }
                         ,
                         {
                           label = 'wall_6_in',
                           kind = 'precice',
                         --  pressure = '10000',
                         --  density = 1,
                         --  v_x = 0.0,
                         --  v_y = 0.0,
                         --  v_z = 0.0
                         },
                       }
