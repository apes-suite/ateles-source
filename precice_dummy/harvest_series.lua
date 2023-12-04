name = 'gPulseDens_euler_modg'

-- define the input
input = {
          read = 'restart/precice/gPulseDens_euler_modg_header_177.150E-06.lua',

         -- define the subsampling parameters
         subsampling = {
                         levels = 1,  -- the number of subsampling steps, i.e. 1 means we increase the number of elements by 8
                         --projection = 'QLegendre', -- the type of projection we use for the subsampling of the data.
                         --dofReductionFactor = 1.5,  -- the reduction factor for the dofs of each subsampling step (per spatial dir)
                         projection = 'QLegendrePoint', -- the type of projection we use for the subsampling of the data.
                       }
         }


-- define the output
output = {  -- some general information
            folder = 'harvest/precice/',     -- Output location 
           {

    
            format = 'VTU',   -- Output format 

            binary = true,
            vrtx = {           -- when this table is defined set use_vrtx = .true.
--                     updated_vertices  = 'prefix_of_new_vertices',            -- when this is given set  new_vrtx = .true.
--                     old_vertices  = 'prefix_of_old_vertices',                -- when this is given set  old_vrtx = .true.
--                     dump_vertices = 'prefix_of_where_to_dump_vertices'       -- when this is given set dump_vrtx = .true.
                   } 
           }    

         }


