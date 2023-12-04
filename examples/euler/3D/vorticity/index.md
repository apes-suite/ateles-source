title: Vorticity Test Case

This testcase for Euler 3D investigates the vorticity. We consider vorticity.
Vorticityi is the tendency for elements of the fluid to "spin.“.
Vorticity can be related to the amount of “circulation” or "rotation" 
(or more strictly, the local angular rate of rotation) in a fluid.
The simulation domain is a cube with periodic boundaries.

**Features used**

1. Projection: fpt, Blocksize 32

2. Polynomial representation: Q

3. Filtering: -- 

4. Timestepping: explicitRungeKutta, 4 steps 

5. Boundary conditions: -- 

6. Others: 
   - Derived quantity: Momentum, velocity, vorticity

