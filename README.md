# Discontinuous Galerkin Method for the 2D Heat Equation

## Introduction

This project implements a **Discontinuous Galerkin (DG) method** in Mathematica and C++ for solving the two-dimensional heat equation

$$
u_t - \Delta u = 0,
\qquad (x,y)\in[0,1]^2,\quad t\in[0,1].
$$
using a DG method. We note that the C++ code runs much faster (as expected). Mathematica was mainly used to get the computations for the coefficients of the DG method, since many integrals were involved.
