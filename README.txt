1.What is the effect of the Ks and Kd parameters on the jello? 

Ks will pull particles on two sides of the spring together, trying to keep their distance at a certain leval.
Kd will damp the particles from moving, trying to stop all velocities within the systems.

2. What are the benefits and the drawbacks of the collision system used here? What are some different ways in which it could be improved?

The benefit is that the particle should never go deep inside an obstacle, as Collision test will apply a force to pull it out from an obstacle and Contact will bring it out of the obstacle directly.

3. From lecture, What is an example of a stiff constrained system?

Stiff constrained system will bound all elements together and only allow for minor relative displacements.

4. From lecture, What is the difference between and explicit and implicit integration scheme?

In explicit scheme, we know every value on the left side of an equation and calculate along time while in implicit integration scheme, there exist unkown values on the left hand side, thus we need to estimate it using numerical methods.

5. Does the jello behave realistically? 

Mostly yes. It stays relatively stable and stiff and can react to external forces and collision forces.

About the project
I finished the first extra credit problem by implementing a sphere in the world. The jello can detect the sphere can collide with it.

I tried the second extra credit problem, my thought is :
1. Generate the mouse position once left-clicked with ALT key down;
2. Unproject the mouse into the World coordinate.
3. Detect collision of the ray with world particles with a range of 0.1. (Cout the number of particle and I realize the numbers are not correct, may be the ray is not well unprojected.)
4. Choose the particle closest to the screen and change it's position using mouseMotion function.

