# Nonlinear-Finite-Element-Code-for-2D-Problems
A program for solving nonlinear boundary value problems in solid mechanics with nonlinear material behaviour and large deformation. To start, modify and run 'Run_code.mlx'.

## Instruction
The code provides both static and dynamic (alpha method) solver for 2-D boundary value problem with the consideration of large deformation.



Formulation: Updated Lagrangian & Total Lagrangian

Element: Quadrilateral Element

Numerical Quadrature: 2-Points Gauss–Legendre Quadrature

Material: Neo-Hookean Material

Quadratic convergence rate has been achieved with consistent tangent matrix.

## Function List
Run_code.mlx

Finite_Derivative.m 

Finite_Deformation.m

Finite_Stress_Tangent.m

Finite_Internal_Force_Tangent.m

Evaluate_Integral.m

Solver.m (statics)

Solver_for_alpha.m (dynamics, alpha method)

## Flowchart
![image](https://user-images.githubusercontent.com/112973740/215729485-df2b4e52-0fe8-4efc-bafa-56a9231ea4ab.png)

## Reference
[1] Ted Belytschko, Wing Kam Liu, Brian Moran, Khalil I. Elkhodary (2014) Nonlinear Finite Elements for Continua and Structures

[2] Thomas J.R. Hughes (2000) The Finite Element Method
