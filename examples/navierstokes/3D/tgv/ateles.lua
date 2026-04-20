-- Configuration for Taylor-Green Vortex --
-- This setup implements the Taylor-Green Vortex in 3D that is simulated
-- with the compressible Navier-Stokes equations.

sim_name = 'tgv_nvrstk_modg_3d'

-- Variables to be set for the simulation --
--
-- Reference Reynolds number:
Re = 1600

-- Reference Mach number:
Ma = 0.1

-- Reference Prandtl number:
Pr = 0.71

-- Reference length [m]:
L = 0.001524

-- Fluid properties (ISA standard atmosphere 1976, Troposphere)
-- Reference density [kg/m^3]:
rho_0 = 1.225

-- Reference pressure [Pa]:
p_0 = 101325

-- Reference Temperature [K]:
T_0 = 288.15

-- Ratio of specific heats (isentropic expansion coefficient):
gamma = 1.4


-- Derived quantities --
--
-- Speed of sound:
SoS = math.sqrt(gamma*p_0/rho_0)

-- Velocity:
V_0 = Ma*SoS

-- Ideal gas constant:
R = p_0 / (rho_0 * T_0)

-- Characteristic time:
t_c = L/V_0

-- Initial conditions
function ini_press(x,y,z)
  return p_0 + rho_0*V_0^2/16 * (math.cos(2*x/L) + math.cos(2*y/L)) * (math.cos(2*z/L) + 2)
end
function ini_dens(x,y,z)
  return ini_press(x,y,z) / (R*T_0)
end

function ini_vel_x(x,y,z)
  return V_0 * math.sin(x/L)*math.cos(y/L)*math.cos(z/L)
end
function ini_vel_y(x,y,z)
  return -V_0 * math.cos(x/L) * math.sin(y/L) * math.cos(z/L)
end
ini_vel_z = 0.0

-- logging, with higher level, the solver writes out
-- more information regarding the settings.
logging = {level=10}

-- Control the Simulation times --
-- Set the simultion time and when information should write out.
-- The max defines the max. simulation time, the min. where to
-- start the simulation and the interval, after how many iteraition
-- information should be written out to the output file.
sim_control = {
  time_control = {
    max = 20 * t_c,
  },
  abort_criteria = {
    stop_file = 'stop'
  }
}

-- Check for Nans and unphysical values
check = { interval = 1 }

-- Mesh configuration --
mesh = {
  predefined = 'cube',
  origin = { 0.0, 0.0, 0.0 },
  length = 2*math.pi*L,
  refinementLevel = 2
}

-- Equation definitions --
-- We solve the Navier-Stokes 3D equations
-- Besides the fluid parameters, we need to set an internal penalization
-- parameter that is used for the implementation of the viscous fluxes.
equation = {
  name       = 'navier_stokes',
  isen_coef  = gamma,
  r          = R,
  -- internal penalization for viscous fluxes
  ip_param   = 4.0,
  material = {
    characteristic = 0.0,
    relax_velocity = {0.0, 0.0, 0.0},
    relax_temperature = 0.0
  }
}
equation["mu"] = rho_0*V_0*L/Re
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)
equation["therm_cond"] = equation["mu"] * equation["isen_coef"] * equation["cv"] / Pr

-- Scheme definitions --
-- In the spatial discretization we have to use the modg scheme for 3D.
-- In time we use the classical Runge-Kutta 4 stage scheme with adaptive
-- timesteps that are chosen according to the CFL condition, given by
-- the courant factor as 'cfl'.
-- The viscous timestep limit is considered by the factor in 'cfl_visc'.
scheme = {
  spatial =  {
    name = 'modg',
    m = 3
  },
  temporal = {
    name = 'explicitRungeKutta',
    steps = 4,
    control = {
      name     = 'cfl',
      cfl      = 0.8,
      cfl_visc = 0.6
    }
  }
}

-- Projection type --
-- We consider here fast polynomial
-- transformation with an oversampling
-- factor of 2, which means that we
-- use two times more points for the
-- approximation
projection = {
  kind = 'fpt',
  factor = 1.0
}

-- Define the inital conditions --
-- We need to set density, pressure and the velocity in x and y.
initial_condition = {
  density = ini_dens,
  pressure = ini_press,
  velocityX = ini_vel_x,
  velocityY = ini_vel_y,
  velocityZ = ini_vel_z,
  useFpt = true
}

-- Restart --
--
restart = {
  read = 'restart/lastHeader.lua',
  init_on_missing = true,
  write = 'restart/',
  time_control = {
    min = 0,
    max = sim_control.time_control.max,
    interval = t_c,
    align_trigger = true
  }
}

-- Tracking --
-- We track here a the momentum in one element (selected by the defined point).
-- After each iteration, the current value of the tracked variable will be
-- written on a new line in the corresponding file.
tracking = {
  label = 'track_momentum',
  variable = { 'momentum' },
  shape = {
    kind = 'canoND',
    object= {
      origin ={ mesh["length"]/4, mesh["length"]/4, mesh["length"]/4 }
    }
  },
  time_control = {
    -- We write the tracked variable every iteration over the complete
    -- simulation time.
    min = { iter = 0 },
    max = sim_control.time_control.max,
    interval = { iter = 1 }
  },
  output = { format = 'ascii', use_get_point = true }
}

-- Write out a timing file which includes
-- all timings of the solver
timing = { filename = 'timing.res' }
