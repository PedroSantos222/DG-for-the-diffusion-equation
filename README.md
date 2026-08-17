This project implements a **Discontinuous Galerkin (DG) method** in Mathematica and C++ for solving the two-dimensional heat equation

$$
u_t - \Delta u = 0,
\qquad (x,y)\in[0,1]^2,
\qquad t\in[0,1].
$$

The C++ implementation runs significantly faster than the Mathematica implementation. Mathematica was mainly used to derive the DG formulation and compute the required coefficients and matrix entries, which involve several symbolic integrals.

We first consider the **stationary problem**, where the DG discretization leads to a linear system that can be solved to obtain the numerical solution. We then introduce the time dependence. After spatial discretization, the PDE is transformed into a **semi-discrete system of ordinary differential equations (ODEs)**, which is solved using the forward Euler method to obtain the numerical solution of the PDE.

