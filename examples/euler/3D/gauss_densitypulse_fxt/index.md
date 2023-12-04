title: Euler 3D setup of Gaussian pulse in density with FXT

This configuration shows the use of the FXT projection method. This
method uses the fast multipole method implemented in FXTPACK for the
projection between modal and nodal representations.
A simple setup with a pulse in density that is convected through the
periodic domain is used for this example.

The configuration is found in `ateles.lua`:

```lua
{!examples/euler/3D/gauss_densitypulse_fxt/ateles.lua!}
```

**Features used**

1. Projection: fxt, Oversampling 2.0

2. Polynomial representation: Q

3. Filtering: --

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: --
