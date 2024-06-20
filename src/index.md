This repository contains the teaching material of the course
[MATH-500 Error control in scientific modelling](https://staging-edu.epfl.ch/coursebook/en/error-control-in-scientific-modelling-MATH-500) at EPFL.

## Outline

1. [Introduction](https://epfl-matmat.github.io/error-control-modelling/01_Introduction.html)
1. [The Dirichlet Laplacian in an Unbounded Domain](https://epfl-matmat.github.io/error-control-modelling/02_Laplace_error_sources.html)
1. [Matrix Eigenproblems](https://epfl-matmat.github.io/error-control-modelling/03_Matrix_eigenproblems.html)
1. [Bounds on Eigenvalues](https://epfl-matmat.github.io/error-control-modelling/04_Matrix_error_bounds.html)
1. [Errors due to floating-point arithmetic](https://epfl-matmat.github.io/error-control-modelling/05_Floating_point_arithmetic.html)
1. [Diagonalisation algorithms](https://epfl-matmat.github.io/error-control-modelling/06_Diagonalisation_algorithms.html)
1. [Matrix Pertubation Theory](https://epfl-matmat.github.io/error-control-modelling/07_Pertubation_theory.html)
1. [Hilbert Spaces](https://epfl-matmat.github.io/error-control-modelling/08_Hilbert_spaces.html)
1. [Operators](https://epfl-matmat.github.io/error-control-modelling/09_Operators.html)
1. [Self-adjoint operators](https://epfl-matmat.github.io/error-control-modelling/10_Spectra_self_adjoint.html)
1. [Periodic Problems](https://epfl-matmat.github.io/error-control-modelling/11_Periodic_problems.html)
1. [Nomenclature](https://epfl-matmat.github.io/error-control-modelling/nomenclature.html)

If you spot an error feel free to make a pull request to the
[github repository](https://github.com/epfl-matmat/error-control-modelling)
generating this website.

## Summary
Errors are ubiquitous in computational science as neither models nor numerical
techniques are perfect. With respect to eigenvalue problems motivated from
materials science and atomistic modelling we discuss,
implement and apply numerical techniques for estimating simulation error.

## Content
* Important eigenvalue problems in materials science
* Motivation for studying errors in eigenvalue problems
* Types of simulation error
* Residual-error relationships for eigenvalue problems
* Perturbation theory and parametrised eigenvalue problems
* Subtleties of infinite-dimensional eigenvalue problems
* Discretisation and discretisation error
* Plane-wave basis sets
* Errors due to uncertain parameters (if time permits)
* Non-linear eigenvalue problems (if time permits)

Algorithm demonstrations and implementations will be based on
the [Julia programming language](https://julialang.org/)
and interactive [Pluto](https://plutojl.org/) notebooks.

## Prerequisites
* Analysis
* Linear algebra
* Exposure to numerical linear algebra
* Some experience with numerical methods for solving differential equations
  (such as finite-element methods, finite-difference approaches, plane-wave methods)
* Exposure to implementing numerical algorithms (e.g. using Python or Julia)

This course delivers a mathematical viewpoint on materials modelling and it is
explicitly intended for an interdisciplinary student audience. To keep it
accessible, the key mathematical and physical concepts will both be revised as
we go along. However, the learning curve will be steep and an interest to learn
about the respective other discipline is required. The problem sheets and the
projects require a substantial amount of work and feature both theoretical
(proof-oriented) and applied (programming-based and simulation-based) components.
While there is some freedom for students to select their respective
focus, students are encouraged to team up across the disciplines
for the course work.

## Literature
There is no single textbook corresponding to the content of the course.
Parts of the lectures have substantial overlap with the following resources,
where further information can be found.

- Youssef Saad. *Numerical Methods for Large Eigenvalue Problems*, SIAM (2011).
  A PDF is available free of charge on
  [Youssef Saad's website](https://www-users.cse.umn.edu/~saad/eig_book_2ndEd.pdf).
- Nicholas J. Higham. *Accuracy and Stability of Numerical Algorithms*, SIAM (2002).
- Peter Arbenz. *Lecture notes on solving large scale eigenvalue problems*, ETHZ.
  A PDF is available from
  [Peter Arbenz' website](https://people.inf.ethz.ch/arbenz/ewp/Lnotes/lsevp.pdf).
- Mathieu Lewin. *Théorie spectrale et mécanique quantique*, Springer (2022).
