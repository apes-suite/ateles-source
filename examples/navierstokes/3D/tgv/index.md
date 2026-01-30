title: Taylor-Green Vortex (TGV) 3D

This 3D Navier-Stokes setup utilizes an initial vortex configuration,
widely used to assess the behavior of numerical schemes with respect
to turbulence.

See for example [Solutions of the Taylor-Green Vortex Problem Using
High-Resolution Explicit Finite Difference Methods](https://ntrs.nasa.gov/api/citations/20130011044/downloads/20130011044.pdf)

The problem definition is a cube with full periodicity in
all directions.
It domain length is defined with: 
$$ 2 \pi L $$

And \(L\) indicates the reference length of the problem.
The isothermal initial condition is defined as:

\begin{equation}u = V_0 \sin\left(\frac{x}{L}\right) \cos\left(\frac{y}{L}\right) \cos\left(\frac{z}{L}\right)\end{equation}
\begin{equation}v = -V_0 \cos\left(\frac{x}{L}\right) \sin\left(\frac{y}{L}\right) \cos\left(\frac{z}{L}\right)\end{equation}
\begin{equation}w = 0\end{equation}
\begin{equation}p = p_0 = \frac{\rho_0 V_0^2}{16} \left( \cos\left(\frac{2x}{L}\right) + \cos\left(\frac{2y}{L}\right)\right)\left(\cos\left(\frac{2z}{L} + 2\right)\end{equation}

From the ideal gas assumption the initial condition in density follows accordingly:
\begin{equation}\rho = \frac{p}{R T_0}\end{equation}

Where \( (u,v,w,p,\rho) \) are the primitive state variables of the compressible
Navier-Stokes equations velocity in x, y, z, pressure and density.
\( V_0 \) is the maximum velocity in each component,
\( p_0 \) is the reference pressure,
\( \rho_0 = \frac{p_0}{R T_0} \) is the reference density,
\( T_0 \) is the reference temperature,
\( R \) is  the gas constant,
and \( L \) is the reference length.

The reference Reynolds number is given by:

$$ Re = \frac{\rho_0 V_0 L}{\mu} $$

where \(\mu\) is the dynamic viscosity of the fluid.

And the reference Mach number is given by:

$$ Ma = \frac{V_0}{C_0} $$

where \(C_0\) is the reference speed of sound given by:

$$ C_0 = \sqrt{\gamma R T_0} $$

with \( \gamma = \frac{c_p}{c_v}\) as the relation of specific
heat capacities of the fluid.

Finally the reference Prandtl number is defined by

$$ Pr = \frac{\mu c_p}{\kappa} $$

where \( \kappa \) is the heat conductivity of the fluid.

The characteristic time is given by reference velocity and the
reference length: \( t_c = \frac{L}{V_0} \).

The configuration is found in `ateles.lua`:

```lua
{!examples/navierstokes/3D/tgv/ateles.lua!}
```

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: -

4. Timestepping: explicitRungeKutta, 4 steps

5. Boundary conditions: -

6. Others: 
   - Derived quantity: Momentum
