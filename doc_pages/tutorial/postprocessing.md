title: Postprocessing with atl_harvesting

After creating some restart files, we want to postprocess those using
our postprocessing tool called atl\_harvesting. It is tightly coupled with the
solver to share as much code as possible, which enables us to use some
sophisticated postprocessing features.

First we desribe the general procedure for running it.
Next we show an example to postprocess the
[Maxwell Testcase](maxwell.html) with atl\_harvesting.
In the end we show how to run atl\_harvesting in series
and generate a visualization file where all the timesteps are combined
and can easily be visualized in [Paraview](http://www.paraview.org).

[TOC]

## Run atl\_harvesting

atl\_harvesting is a postprocessing tool specific for the Ateles solver.
It is build along Ateles into the same build directory `ateles/build/`.
So make sure you successfully compiled Ateles.
To run atl\_harvesting from the ateles main directory, type

```bash
build/atl_harvesting harvester.lua
```

where `build/` is the relative path
from the current working directory to the executable `atl_harvesting`
and `harvester.lua` is the name of the input file
(the default file name atl\_harvesting is looking for is `harvester.lua`).
Similar to Ateles, the input file for atl\_harvesting is written in Lua
and defines the data to process, optional subsampling and the output format.
In the main `ateles` directory
there is an example configuration file named `harvester.lua`.

## Example: Maxwell Testcase

Now we desrcibe an example for postprocessing the Maxwell Testcase.
For this purpose we describe
the atl\_harvesting input file `harvester.lua`
located in the directory `tutorial/testcase/maxwell` in the Ateles main
directory.
We stored the simulation result from the previous step in this folder, and we
are going to refer to the configuration files within this folder, so it makes
life pretty easy when we just change to this folder. The only thing we have to
keep in mind is that we now need to refer to atl\_harvesting with the relativ
path `../../../build/atl_harvesting`.

### Load Configuration of Simulation

First we load the configuration of the original simulation
defined in the Ateles input file `ateles.lua` in the same directory.

```lua
-- Use the configuration of the original simulation run.
require 'ateles'
```

This enables us to use the variables defined in there.

### Set Data to Process

We define the restart data
which we want to process in the parameter `restart.read`.
The correct path and name to the restart files should be defined there.

```lua
-- Set the restart data to harvest.
restart.read = 'restart/' .. simulation_name .. '_lastHeader.lua'
```

Here we are going to use the last timestemp stored in the `restart` folder. The
variable simulation\_name is defined in `ateles.lua`.

### Subsampling

Subsampling is a procedure to get more visualizations points
and thus a better picture of the domain.
The variable `nlevels` defined in the table `ply_sampling`
specifies the maximum number of subsampling steps
(maximum number of elements for visualization should be increased).

```lua
-- Subsampling (for subresoloved color information):
ply_sampling = {
  nlevels = 3,  -- maximal level to use in subsampling
}
```

In our example we increase the number of levels by 3 at most, which results in
\(2^3\) more elements per direction.

### Output

With the parameter `label` 
one can define the beginning of the output files' name.
The table `variable` defines which parameters of the simulation to process.
For our Maxwell testcase we choose `displacement_field` and `magnetic_filed`.
Which variables are available depends on the solved equation system. See
the description of the equation systems for more information about available
variables.
The parameter `folder` refers to the folder where the output files should be
stored. It is important that this folder needs to be created in advance
(similar to the restart folder).
If the folder does not exist, the program will fail with an error.
The variable `format` defines the format for the output files.
It can be
'[vtk](http://www.vtk.org)',
'[ascii](https://en.wikipedia.org/wiki/ASCII)',
'asciiSpatial', 'harvester' or '[precice](http://www.precice.org)'.
We want to visualize our data and set it to 'vtk'.

```lua
-- Example tracking to generate vtk files:
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
```

After running atl\_harvesting with the `harvester.lua` input file
we get a .vtu file in the `harvest` folder.
This file can be visualized with Paraview.

## Run in Series

In this section we describe
how to run atl\_harvesting in series for multiple restart files
and create one visualization file for all timesteps
using the harvest\_series.py script.
The harvest\_series.py script is part of
[TreElm](../../../../doxy/treelm/index.html)
and located in `<ateles>/treelm/peons`.

###  Set up Template

First harvest\_series.py needs a template file to construct
the configuration for atl\_harvesting.
Thus, it should resemble the `harvester.lua` configuration
described above, but may use some replacement markers that
will be filled in by harvest\_series.py

Typically it may look like this:

```lua
{!series_example/harvester.template!}
```

Above is the `harvester.template` located in `series_example` in the Ateles
main folder.
The strings with `$!`at the beginning and `!$` at the ending
are placeholders for variables defined by the harvest\_series.py script.

### Set up Variables

The main arguments to the harvest\_series.py are the restart files to process.
However, several other options can be provided as command line arguments.
These options can also be put into a configuration script, so they do not
have to be stated explicitly every time.
Per default the script will look for this configuration in `series.config`,
but this can be changed via the `-c` option.
A description of all options is available in the help provided by the
`-h` option. Here are some of the major settings:

-   `template`:
    name and path of the template
-   `files`:
    name and path of the input restart files (replaces `'$!filename!$'`)
-   `out`:
    folder for dumping the output (replaces `'$!folder!$`)
-   `harvester`:
    path to the atl\_harvesting executable
-   `run`:
    command to run in parallel

An example can be found in `series_example/series.config` which is:

```
{!series_example/series.config!}
```

Lines starting with a `#` are ignored.

### Run harvest\_series.py

You can run the `harvest_series.py` from main the `ateles` directory with

```bash
python treelm/peons/harvest_series.py
```

where `treelm/peons/` is the relative path to the `harvest_series.py` script
from the main `ateles` directory.

Thus, if you want to postprocess the default `ateles.lua` testcase in the
ateles dirctory you can run:

```bash
python treelm/peons/harvest_series.py -c series_example/series.config
```

in the ateles directory.
The output will be stored for that example in the directory `t_series`, which
is the default for the `out` option.
In this directory you will then find all `.vtu` files for each timestep
and one `.pvd` file (visu.pvd), where all timesteps are collected.
You can open the `.pvd` file in Paraview
to visualize the simulation as a timeseries.
In the .pvd file there is just a linking to the .vtu files,
where the information about the timestep is stored.
So make sure the .pvd and the collected .vtu files are in the same folder.

### Example

Assume we want to postprocess the restart files from the
[Maxwell tutorial](maxwell.html). The restart files are stored in
`<ateles>/tutorial/testcase/maxwell/restart`. When we want to run harvest_series
from the folder `<ateles>/tutorial/testcase/maxwell`, we need to change the
paths in series.config to point to the correct files:

```
files: restart/posci_modg_header_*.lua
lua: ../../../build/debug/treelm/aotus/lua
harvester: ../../../build/debug/atl_harvesting
template: ../../../series_example/harvester.template
```

The harvester.template doesn't need any modifications. To run it, some more
relativ paths are used:

```bash
python ../../../treelm/peons/harvest_series.py -c ../../../series_example/series.config
```

The result is stored in `<ateles>/tutorial/testcase/maxwell/t_series`.
