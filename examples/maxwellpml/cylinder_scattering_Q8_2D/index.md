title: Perfectly matched Layer test case

A perfectly matched layer (PML) is an artificial absorbing layer for wave equations, 
commonly used to truncate computational regions in numerical methods to simulate problems 
with open boundaries, especially in the FDTD and FE methods. The key property of a PML that
distinguishes it from an ordinary absorbing material is that it is designed so that waves 
incident upon the PML from a non-PML medium do not reflect at the interface—this property
allows the PML to strongly absorb outgoing waves from the interior of a computational region
without reflecting them back into the interior.

Test case description: This setup investigates the behaviour of the PML for the maxwell equations by 
means of a scattering cylinder. We consider here a homogenous media. 

This test case is tailored to check the Q-polynomial representation and a scheme order of 8

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: -- 

4. Timestepping: explicitRungeKutta, 4 steps 

5. Boundary conditions: --  

6. Others: 
   - Derived quantity: Displacement field
   - Source term: Current Density
