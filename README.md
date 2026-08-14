# Discontinuous Galerkin Method for the 2D Heat Equation

## Introduction

This project implements a **Discontinuous Galerkin (DG) method** in Mathematica for solving the two-dimensional heat equation

$$
u_t - \Delta u = 0,
\qquad (x,y)\in[0,1]^2,\quad t\in[0,1].
$$

Equivalently,

$$
\frac{\partial u}{\partial t}
-
\frac{\partial^2 u}{\partial x^2}
-
\frac{\partial^2 u}{\partial y^2}
=0.
$$

The initial condition is chosen as

$$
u(x,y,0)=\sin(\pi x)\sin(\pi y).
$$

The spatial domain is divided into a uniform $3\times3$ mesh, with a bilinear ($Q_1$) polynomial approximation on each element. Since the DG approximation is allowed to be discontinuous across element boundaries, numerical fluxes and an interior penalty term are introduced to couple neighboring elements.

The spatial discretization leads to a system of ordinary differential equations

$$
M\mathbf{U}'(t)+A\mathbf{U}(t)=0,
$$

where $M$ is the global mass matrix, $A$ is the DG stiffness/interface matrix, and $\mathbf{U}(t)$ contains the coefficients of the DG approximation.

The initial condition is projected onto the DG finite-dimensional space by solving

$$
M\mathbf{U}(0)=\mathbf{b},
$$

where the entries of $\mathbf{b}$ are given by

$$
b_k = \int_{\Omega} u_0(x,y)\phi_k(x,y)\,dx\,dy.
$$

The resulting ODE system can be written as

$$
\mathbf{U}'(t)=B\mathbf{U}(t),
\qquad
B=-M^{-1}A.
$$

Forward Euler is then used to advance the solution in time:

$$
\mathbf{U}^{n+1}
=
\mathbf{U}^{n}
+
\Delta t\,B\mathbf{U}^{n}.
$$

Finally, the DG coefficients are used to reconstruct the approximate solution $u_h(x,y,t)$ on each element.

For the chosen initial condition and homogeneous Dirichlet boundary conditions, the exact solution is

$$
u(x,y,t)
=
e^{-2\pi^2t}\sin(\pi x)\sin(\pi y),
$$

which provides a reference solution for evaluating the numerical approximation.
