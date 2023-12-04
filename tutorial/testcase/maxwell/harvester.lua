-- use the configuration of the original simulation run
require 'ateles'

-- set the restart data to harvest
restart.read = './restart/' .. simulation_name .. '_lastHeader.lua'

-- Subsampling (for subresoloved color information):
ply_sampling = {
  nlevels = 3    -- maximal level to use in subsampling
}

-- tracking to generate vtk files:
tracking = {
  {
    label = 'visu',
    variable = {'displacement_field', 'magnetic_field'},
    shape = {
      kind = 'global'
    },
    folder = 'harvest/',
    output = {
      format = 'vtk'
    }
  }
}
