project: Ateles
summary: Discontinuous Galerkin solver with explicit time stepping.
project_website: https://apes.osdn.io/pages/ateles
src_dir: source/
src_dir: build/ford
src_dir: polynomials/source
exclude_dir: build/ford/treelm
external: aoturl = https://geb.inf.tu-dresden.de/doxy/aotus
external: temurl = https://geb.inf.tu-dresden.de/doxy/treelm
output_dir: docu
page_dir: doc_pages
copy_subdir: media
graph: true
graph_maxdepth: 4
graph_maxnodes: 32
display: public
display: protected
display: private
sort: permission
source: false
author: University of Siegen
title: Mainpage
print_creation_date: True
md_extensions: markdown.extensions.toc

(Adaptive Tree based Efficient and Lithe Equation Solver)

Ateles implements a modal/nodal Discontinuous Galerkin scheme on top of
[TreElm](|temurl|/index.html) data structures.
The dedicated mesh generator [Seeder](https://apes.osdn.io/pages/seeder.html)
provides the possibility to create meshes for the solver.
Various equation systems are supported, each with their dedicated optimized
kernel.
Ateles uses an explicit time integration.

The solver is configured with the help of Lua scripts.
Your simulation setup has to be described in such a script and the name of
the script file has to be passed to the solver as its sole command line
argument. If no argument is provided, the solver will try to open a file
named `ateles.lua` in the current working directory.
Thus, running Ateles takes the following form:

```sh
ateles myconfig.lua
```

To get more information about Ateles and how to use it, please take a look into
the [Documentation](|page|).
See especially the [Examples](|page|/examples) there for specific setups for the
various supported equation systems.

## Downloading Ateles

We use Mercurial (hg), an open source distributed version control system,
for the APES projects.
To obtain the
[sources of Ateles](https://osdn.net/projects/apes/scm/hg/ateles/), you need to
have Mercurial available.
Mercurial may be downloaded from the
[official website](https://www.mercurial-scm.org)
or already available through some packaging system for your machine.

Cloning the Ateles repository will also get you the
[TreElm library](https://osdn.net/projects/apes/scm/hg/treelm/)
as a subrepository in the directory `ateles/treelm`.
Note that you might want to get other supporting tools from the
[APES framework](https://osdn.net/projects/apes/)
or the complete framework itself to make use of other solvers for your
simulations.
