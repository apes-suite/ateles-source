title: Periodic_Oscillator

Test case description: This setup investigates the periodic oscillation, solving the maxwell equations. We consider here a homogenous media and a periodic simulation domain, 
hence no boundaries are taken into account. The testcase simulates the temporal development of standing waves for electric and magnetic fields. 
it is a well-known that the approximation of the particle motion introduces numerical errors, hence the charge conservation equation is not satisfied on the dicrete level. 
To avoid the increase of errors in Gauss' law, a divergence cleaning step for a correction potential is often added. This correction is used for this testcase

This test case is tailored to check the P-polynomial representation and a scheme order of 8

**Features used**

1. Projection: l2p, Oversampling of 2

2. Polynomial representation: P

3. Filtering: -- 

4. Timestepping: explicitRungeKutta, 4 steps 

5. Boundary conditions: --  

6. Others: 
   - Derived quantity: Displacement field
