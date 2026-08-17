This project implements a **Discontinuous Galerkin (DG) method** in Mathematica and C++ for solving the two-dimensional heat equation

$$
u_t - \Delta u = 0,
\qquad (x,y)\in[0,1]^2,
\qquad t\in[0,1].
$$

The C++ implementation runs significantly faster than the Mathematica implementation, as expected. Mathematica was mainly used to derive and compute the coefficients of the DG method, since many symbolic integrals are involved.

We start by considering the equation without the time component, and solve a linear system to obtain the numerical solution. Then we consider the time dependent case, and obtain an ODE, for which we use the Euler method to obtain the numerical solution of the PDE.
