# ODEInterface

[![Build Status](https://travis-ci.org/luchr/ODEInterface.jl.svg?branch=master)](https://travis-ci.org/luchr/ODEInterface.jl)

This julia module provides an interface to solvers for 
ordinary differential equations (ODEs) written in Fortran
for solving initial value problems of the form

    x' = rhs(t,x),      x(t₀) = x₀

or (for solvers supporting a "mass matrix" M)

    M⋅x' = rhs(t,x),    x(t₀) = x₀.

## What does "Interface" mean?

This julia module does *not* contain code for solving initial value
problems, but this module does contain code for interacting with
compiled Fortran-solvers. That's the reason, why this module is not called
ODESuite.

## What solvers are currently supported?

Currently the following Fortran-solvers, written by
Prof. E. Hairer and Prof. G. Wanner, are supported:

* dopri5: explicit Runge-Kutta method of order 5(4) due to Dormand & Prince
* dop853: explicit Runge-Kutta method of order 8(5,3) due to Dormand & Prince
* odex: GBS extrapolation-algorithm based on the explicit midpoint rule
* radau5: implicit Runge-Kutta method (Radau IIA) of order 5
* radau: implicit Runge-Kutta method (Radau IIA) of variable order between 5 and 13
* seulex: extrapolation-algorithm based on the linear implicit Euler method

see [Software page of Prof. Hairer](http://www.unige.ch/~hairer/software.html).

Description: [Calling the Solvers](./doc/CallSolvers.md) 

The following features of this solvers are supported by this ODEInterface:

* providing an [output function](./doc/OutputFunction.md) (e.g. 
for dense output or for event location) to the solvers
* providing mass- and jacobi-matrices for the solvers (with support for
banded matrices)
* all the solvers' parameters for fine-tuning them, 
see [Options for Solvers](./doc/SolverOptions.md)
* support for problems with "special structure", 
see [special structure](./doc/SpecialStructure.md)

## What are the requirements for this module

In order to use this module, you have to *compile* the supported
Fortran solvers and provide a shared library for each solver. Just call
`ODEInterface.help_solversupport` for further informations (help topics)
on how to compile the solvers and how to create shared libraries.

## Further help

see `ODEInterface.help_overview` for an overview of some help topics. 

## Contacting the author of this module

The author of this julia module is 

     Dr. Christian Ludwig
     email: ludwig@ma.tum.de
       (Faculty of Mathematics, Technische Universität München)

