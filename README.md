# The Heat Equation Solver

## Description
The Heat Equation Solver is a Python project that explores the solution of the parabolic partial differential equation known as the Heat Equation. This equation describes the conductance of temperature through a material over time. In this project, we focus on an aluminium bar and use two numerical methods, Explicit and Implicit (Crank-Nicolson), to solve the Heat Equation over a time period of 600 seconds.

## Introduction
The Heat Equation is given by:

∂T(x,t)/∂t = k * ∂^2T(x,t)/∂x^2

Where T(x,t) is the temperature at position x and time t, and k is the thermal conductivity of the material.

In this project, we model the temperature distribution along an aluminium bar over time using a two-dimensional spatio-temporal grid. The bar has a length of 100 cm and a thermal conductivity of 0.835 cm^2/s. The initial temperature distribution is T(x, 0) = 500°C, and the boundary conditions are T(0, t) = 0°C and T(L, t) = 0°C.

## Analytical Solution of the Heat Equation
To validate the numerical methods, we also present the analytical solution of the Heat Equation. We tune the hyper-parameter N to achieve high accuracy and computational efficiency. The analytical solution is given by:

T(x, t) = Σ (2 * sin(n * π * x / L) * exp(-n^2 * π^2 * k * t / L^2)) / (n * π)

where the summation is over odd integers n = 1, 3, 5, ..., N, and L is the length of the bar.

## Explicit Solution of the Heat Equation
### Background
The Explicit method is a finite-difference method used to approximate the solution to parabolic partial differential equations like the Heat Equation. It is easy to implement and involves forward differences in time and central differences in space.

### Implementation
In this project, we use the Forward Time Centered Space (FTCS) method as the Explicit method. The computational molecule for FTCS is applied along the grid and through time to estimate the future state of the temperature based on known states.

### Experiments
We conduct experiments to analyze the stability and accuracy of the Explicit method. While the method is conditionally stable for the Heat Equation, we show that certain step size configurations produce stable solutions. However, small step sizes are required to reduce numerical errors, resulting in higher computational costs.

## Implicit Solution of the Heat Equation
### Background
The Implicit method we use is the Crank-Nicolson (C-N) method, which is second-order accurate for both space and time. It incorporates midpoints to improve accuracy and is unconditionally stable for any spatio-temporal configuration.

### Implementation
In this project, we implement the Crank-Nicolson method using a matrix-vector product to solve a system of algebraic equations. The method is computationally intensive due to the solution of the system, but it offers higher accuracy and stability compared to the Explicit method.

### Experiments
We conduct experiments to demonstrate the accuracy and stability of the Implicit method. The Crank-Nicolson method shows a significant improvement in accuracy compared to the Explicit method, especially when using larger step sizes.

## Conclusion
The Heat Equation Solver project provides a comprehensive analysis of the Heat Equation's numerical solutions using both the Explicit and Implicit methods. It demonstrates the trade-offs between accuracy, stability, and computational efficiency for each method and offers valuable insights into solving parabolic partial differential equations.