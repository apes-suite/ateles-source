-- Use this file as template. Do not modify this file for running some testcases

outputname = 'gauss_pulse'
outputpreview = true 

-- stl_files: for each stl file put a table into
-- this table of tables. Or just a string for
-- the filname, if the other two parameters should
-- get the default: boundary_type = 1, fileformat = 'binary'
-- refinementlevel is the level to which the stl shall be refined
-- it is mandatory


-- boundingbox: two entries: origin and length in this
-- order, if no keys are used
boundingbox = {origin = {0.0, 0.0, 0.0},
               length = 10.0}

-- refinebox: three entries: origin, length and refinementlevel
refinebox = {
             {origin = {0.05, 0.05, 0.05},
              length = {4.9, 4.9, 4.9},
              refinementlevel = 2,
	      deformable = false} }

-- seed: position of seed 
seed = {5.0,
        5.0,
        5.0
       }

-- general refinement levels
-- minrefine: defines minimum refinement level of the entire domain
--
-- maxrefine: defines maximum refinement
-- gets overwritten from madatory STL levels
-- (deprecated - used only for spacer stuff right now)
-- will become default for stl refinement level once
-- it is not mandatory any more
minrefine = 1
maxrefine = 1
