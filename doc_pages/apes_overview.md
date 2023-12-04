title: APES Overview

The APES Suite contains several tools that allow for a flexible computation of
mesh based simulations.
The central component is the [TreElm](treelm/index.html) library that implements
the mesh handling.

In general the workflow in APES is to first generate a mesh (with Seeder), then
run the simulation (for example with Ateles), and finally post-process the
results (harvesting, for Ateles with `atl_harvester`.
Mesh generation can be done by the solver itself for simple meshes.
In this case there is no need to run Seeder beforehand.
Similarily, the post-processing can also already be done online by the solver.
However, including the post-processing in the solver run itself is only
advisable for small setups or smaller subsets from the overall data, as each
process will generate a separate visualization file, which is not suitable for
large scale runs with many processes.

All tools are configured with the help of the Aotus library, which provides
access to Lua scripts.
Lua is a complete scripting language, which allows you to use variables and
functions in the configuration.
This enables us, for example, to explicitly state relations between quantitites.

The components of APES relevant to Ateles:

- [Aotus](https://bitbucket.org/apesteam/aotus): Library to interact with Lua
  scripts.
- [Treelm](https://bitbucket.org/apesteam/treelm): Library to manage the
  octree mesh it also includes the libharvesting for post-processing. Treelm
  incorporates Aotus as a subrepository.
- [Polynomials](https://bitbucket.org/apesteam/polynomials): A library providing
  functions for Legendre polynomials.
- [Seeder](https://bitbucket.org/apesteam/seeder): Mesh generation in the Treelm
  octree format.
- [Ateles](index.html), the Discontinuous Galerkin Solver implemented within
  APES.
