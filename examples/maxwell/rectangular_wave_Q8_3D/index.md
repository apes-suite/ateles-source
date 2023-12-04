title: Rectangle Wave

Test case description: This setup investigates the rectangle wave, solving the maxwell equations. We consider here a homogenous media. The testcase simulates the temporal development of standing waves for electric and magnetic fields. This testcase has an analythic solution, which is given in a lua fuction in the script.

This test case is tailored to check the Q-polynomial representation and a scheme order of 8

**Features used**

1. Projection: fpt

2. Polynomial representation: Q

3. Filtering: -- 

4. Timestepping: explicitRungeKutta, 4 steps 

5. Boundary conditions: Perfectly electric conductor 

6. Others: 
   - Derived quantity: Displacement field
