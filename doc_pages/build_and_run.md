title: Build and Run Ateles

This page provides a brief introduction on building
and using the Ateles solver.

[TOC]

## Compilation

Ateles is a MPI-parallel application and you will need to have
the MPI library installed to compile it.
For compilation the waf build tool is used.
The Fortran compiler is specified using the environment variable `FC`.
Thus you should set this environment variable
to the MPI wrapper on your machine, typically `export FC=mpif90`.
With this variable set, you can configure and compile the code by typing

```bash
bin/waf configure build
```

Now, an executable called ateles got created in the folder `ateles/build/`.
The two steps, configuration and building might also be done separately.
Configuration is done with the folowing command

```bash
bin/waf configure
```

whereas the command

```bash
bin/waf distclean
```

is used to clean the current configuration as well as most build output.

If you want to clean only the build output but keep the configuration in place,
use

```bash
bin/waf clean
```

To keep the configuration tidy, it is recommended to clean it before
creating a new one. Therefore it is adviced to call `configure` together with
`distclean` or `cleanall`.

For using q-polynomials you also need to load the fftw module.
For that the `--fftw_path=$FFTWDIR/` option
in the configuration process has to be set

```bash
bin/waf distclean configure --fftw_path=$FFTWDIR/
```

To compile the code, type

```bash
bin/waf build
```

You can also do a faster compilation without doing unit tests
by building only the exectubale with the command

```bash
bin/waf build --target=ateles
```

In the compilation described above you will obtain
an optimized version of the executable.
If you would like to get an executable
that is suitable for debugging, you can run

```bash
bin/waf debug
```

This will produce the executable with debugging symbols in `build/debug`.

Usually the project needs to be configured only once on a system,
but you might recompile it several times.

## Generating this documentation

This documentation is generated with FORD. To create it you need to
have that tool [installed](https://github.com/Fortran-FOSS-Programmers/ford#installation).

With FORD installed, this documentation can be generated with:

```bash
bin/waf docu
```

The generated documentation is then found in

```bash
bin/waf build/ford/atl/docu
```

## Running

In the main `ateles` directory there is an example configuration file
that documents most parameters you might set for the simulation,
and allows you to do a small simulation run right after compilation
to check if everything worked fine.
The configuration is provided in the form of a Lua script,
and the input is therefore named `ateles.lua`.
The application takes just a single parameter on the command line,
which is the name of the configuration file to load.
If no command line parameter is given, this will default to `ateles.lua`.
Said example configuration can therefore be run by just executing

```bash
./build/ateles
```

in the main project directory.
