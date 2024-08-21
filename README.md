# 2D-Upwind-Finite-Difference-Fortran
This repository contains a Fortran implementation of the 2D upwind finite difference method. The method is used for solving partial differential equations (PDEs) by discretizing the equations on a 2D grid and updating the solution based on the upwind scheme. This approach is particularly useful for problems involving advection, where the direction of the flow is critical in determining the numerical scheme.

Key Features:
Modular Design: The code is organized into subroutines within a Fortran module, making it easy to extend and modify.
Boundary Condition Handling: The code includes explicit handling of boundary conditions on all sides of the grid, including special cases for corners.
Double Precision: The calculations are performed using double precision (real(8)), ensuring numerical accuracy.
How to Use:
1-Compile: Use a Fortran compiler, such as gfortran, to compile the code. For example
gfortran -o upwind upwind.f90
2-Run: Execute the compiled program, passing any necessary parameters.
3-Modify: The code is designed to be easily modified for different grid sizes, time steps, and physical scenarios.
Applications:
This code can be applied to a wide range of problems in computational fluid dynamics, weather modeling, and any field where the transport of quantities like mass, energy, or momentum needs to be simulated.
