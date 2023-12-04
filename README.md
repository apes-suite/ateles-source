Ateles
======

(Adaptive Tree based Efficient and Lithe Equation Solver)

Ateles implements a high-order modal/nodal discontinuous Galerkin solver within
the APES suite.
It is working on a linearized octree and uses efficient data structures
for large scale parallel computations.
Cubical elements allow for efficient numerical schemes, while geometries
can be represented with a penalization method within the elements.

Documentation
-------------

The documentation is generated with
[FORD](https://github.com/Fortran-FOSS-Programmers/ford), and
[available online](https://geb.inf.tu-dresden.de/doxy/ateles/index.html).
You can also generate the documentation with `waf gendoxy`, configuration
of FORD is detailed in `mainpage.md`.

Configuration
-------------

See the examples directory for configuration files.
A configuration file that illustrates the various options is available
in `examples/euler/3D/overview/ateles.lua`.

How to compile
--------------

You'll need the [parent repository](https://github.com/apes-suite/ateles) for compilation.
This repository only provides the sources, see there for brief instructions
on compilation.

How to run
----------

Ateles, like the other tools in the APES suite, uses Lua scripts for the
configuration. A sample configuration is available in `ateles.lua`. Run
ateles by providing it as sole argument the configuration file to use.
If you do not provide any argument, the solver will attempt to use
the file `ateles.lua` in the current directory and fails if it is not
present.

License
-------

Ateles is freely available under the terms of the
[ISC License](https://opensource.org/licenses/ISC).
Please see the *LICENSE* file for its details.

Main author of the original application code is Jens Zudrop.
Please consider citing his Thesis when using Ateles in your work:

Zudrop, Jens. Efficient Numerical Methods for Fluid- and Electrodynamics
on Massively Parallel Systems. Diss. RWTH Aachen University, Shaker 2019.
ISBN 978-3-8440-4407-2.
