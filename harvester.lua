-- Use the configuration of the original simulation run.
-- (to get the equation table properly defined).
require 'ateles'
-- Logging level, required to be > 3 since harvester takes the written
-- time to generate the output name and create the pvd file
logging = {level=3}

-- Apply the stabilizing filter before processing the data further?
-- Set this flag to true to apply the filter defined in the stabilization
-- table to the loaded data before further processing it.
-- Default is to not use a filter upon loading the data.
use_post_filter = false

-- Set the restart data to harvest.
restart.read = 'restart_sample_lastHeader.lua'


-- Subsampling for tracking, define a ply_sampling table to activate subsampling
-- for all tracking objects (except those with use_get_point):
-- Note, if you do set nlevels to 0, there will be no subsampling, and
-- only the first degree of freedom in the original mesh will be written out.
ply_sampling = {
  nlevels = 3,    -- maximal level to use in subsampling
                  -- defaults to 0, which deactivates subsampling

  -- method  = 'adaptive', -- method to use for subsampling
                        -- 'fixed': will refine all elements by nlevels
                        -- 'adaptive': adaptive refinement of the
                        --  mesh based on solution (default).

  -- Settings for the adaptive method:
  -- tolerance = 0.1, -- 'tolerance': value for the variation of
                      --  the solution representet in an element,
                      --  decision if an element will be refined
                      --  is based on this value.

  -- ignore_highmodes = false, -- Whether to ignore modes in the parent
                               -- element that exceed the target polynomial
                               -- during each refinement in each direction.
                               -- This can be used as a simple low-pass filter
                               -- by excluding the highest modes from the
                               -- projection to refined elements.

  -- reduction_mode = 'factor', -- how to reduce modes during refinement
                                -- the following modes are available:
                                -- 'factor': the default, in each refinement the
                                --           number of dofs will be multiplied
                                --           by dof_reduction.
                                -- 'decrement': in each refinement the last
                                --              dof_decrement modes will be cut
                                --              off.
  -- dof_reduction = 0.5, -- For reduction_mode 'factor': Reduction of degrees
                          --  of freedom to avoid too drastic
                          --  increase in memory consumption,
                          --  0.5 means half the polynomial degree
                          --  in every sampling step, 1.0 means no
                          --  dof_reduction.
  -- adaptiveDofReduction = false,
                         -- 'adaptiveDofReduction': The dofs will
                         --  be reduced adaptively based on the
                         --  number of elements that will be
                         --  refined in the current sampling step,
                         --  minimum is 'dof_reduction'.
  -- dof_decrement = 1, -- Number of dofs to cut off in the 'decrement'
                        -- reduction mode. Default is 1.
  -- absUpperBoundLevel = 0 -- 'absUpperBoundLevel': Every element
                            --  on or above this level won't be
                            --  refined, suitable for
                            --  multi-level-meshes.
  filter_element = {
    -- Use some filtering during each refinement step before splitting.
    -- This filtering enables filtering that adapts to the polynomials.
    -- The following strategies are available:
    -- - 'none': this is the default and does not apply any filtering.
    strategy = 'none',

    -- strategy = 'oddfract',
    -- - 'oddfract': This chooses the strength of the filter based on the
    --               spectral energy of the weighted odd modes in relation
    --               to the total weighted spectral energy.
    --               Each mode is weight with the polynomial degree,
    --               resembling something like the derivative.
    --               A spectral viscosity damping will be applied to the
    --               modes with the filter order given by:
    --               min + (max-min+1) * oddfract^fract_exponent
    --               Accordingly you can set the following three parameters:
    -- min_order      =  2, -- default is 2
    -- max_order      = 10, -- default is 10
    -- fract_exponent =  2, -- default is 2
  }
}

-- Example tracking to generate vtk files:
tracking = {
  { label = 'visu',
    variable = {'displacement_field','magnetic_field'},
    shape = {kind='global'},
    folder = 'harvest_',
    -- output format
    -- write_pvd=false important to create complete pvd file by harvester
    -- and not single pvd files by the solver
    output = {format = 'vtk', write_pvd=false}
  }
}
