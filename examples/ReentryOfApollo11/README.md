# Reentry_Apollo_11

The reentry problem for Apollo 11 type space vehicles, following Stoer/Bulirsch 1973; see
J. Stoer, R. Bulirsch, Introduction to Numerical Analysis, 3rd ed., Springer-Verlag, New York, §7.3.7 (pp. 565--572).

Based on the German version written by Prof. Folkmar Bornemann, for the lectures on Numerical Mathematics III at Technical University of Munich in the winter term 2002/2003 using MATLAB and Mex-Interfaces to BVPSOL. Translation to English and Julia
by Vishal Sontakke, April 2016.

Read the [Optimal_Control.pdf](./Optimal_Control.pdf) file for details.

#### Exercises:

1. Change the passage in the Maple script which restricts the control 'u' unnecessarily to the interval [-pi/2, pi/2]. What is the condition of the minimum principle now?
2. What will be the maneuver time and the value of J if the vehicle starts the maneuver with a flight path angle of -6.5 degrees?
3. Which flight path angle at time t=0 optimizes J? How large is the gain over that studied in the lecture (-8.1 degrees)?

If you think hard, no more than a handful of changes to the Maple script and the Julia files are to be made for solving (1)-(3).

#### Solutions:

1. The limitation of the range of u is avoided by
	```    
	U [solve]: = arctan (+ lambda [2] * c [3], + lambda [1] * x [1] * c [2]):
	```
	or
  	```
  	U [solve]: = arctan (-Lambda [2] * c [3], - lambda [1] * x [1] * c [2]):
  	```
    The Pontryagin's Minimum principle is always satisfied for the second option as compared to the first, where it is never satisfied. So we choose the second option. This also eliminates the condition `lambda[1] < 0`.

	How to choose now useful initial guesses of the adjoint variables?

2. The checking of the solution of the auxiliary problem for control values u ​​in the range [-pi/2, pi/2] should be removed. The auxiliary problem should be solved with backward shooting method.

3. One has to choose natural boundary conditions for gamma at t=0, that is, the associated adjoint variable has to be zero at t = 0. This requires the following changes:

	in reentry_bc.jl
    ```c
	r[2] = xa[5];
	```
    in reentry.jl

    before the solution of the optimal control problem, force natural boundary condition by
    ```c
    x[5,1] = 0;
	```

    Also, a few more multiple shooting nodes will not hurt.

    That's it.
