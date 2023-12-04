## harvester input file template
## the filename for read and folder name for output are
## defined by harvest_series.py so they are fixed to
## ${filename}$ and ${folder}$
template = 'harvester.template'

## input restart file which replaces the string
## ${filename}$ in the harvester.template.
## Can be single file or multiple files
#in_files = ['restart/spacer_header*.lua']
in_files = ['restart/precice/*.lua'] 

## output folder where replaces the string 
## ${folder}$ in the harvester.template.
out_folder = 'harvest/precice/' 

## harvester exectuable
path_to_exec = '/home/UNISIEGEN/gk780/harvester/build/harvester'

## command to run in parallel
## run_command = 'mpirun -np 4'
run_command = ''
