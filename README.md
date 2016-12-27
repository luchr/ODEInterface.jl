# ODEInterface

[![Build Status](https://travis-ci.org/luchr/ODEInterface.jl.svg?branch=master)](https://travis-ci.org/luchr/ODEInterface.jl)

This julia module provides an interface to solvers for 
ordinary differential equations (ODEs) written in Fortran
for solving initial value problems (IVP) of the form

    x' = rhs(t,x),      x(t₀) = x₀

or (for solvers supporting a "mass matrix" M)

    M⋅x' = rhs(t,x),    x(t₀) = x₀.

Additionally a boundary value solver (called `bvpsol`) is
supported for boundary value problems (BVP) of the form

    x' = rhs(t,x),      r = bc( xa, xb ) = 0

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
* rodas: Rosenbrock method of order 4(3) (with possibly singular mass matrix)

see [Software page of Prof. Hairer](http://www.unige.ch/~hairer/software.html).

Also supported:

* bvpsol: a boundary value problem solver for highly nonlinear two point
  boundary value problems using either a local linear solver or a global
  sparse linear solver. **Please note: The license for `bvpsol` only 
  covers non commercial use, see [License](./LICENSE.md).**

written by P. Deuflhard, G. Bader, L. Weimann, see
[CodeLib at ZIB](http://elib.zib.de/pub/elib/codelib/en/bvpode.html).

Description: [Calling the Solvers](./doc/CallSolvers.md) 

The following features of the IVP-solvers are supported by this ODEInterface:

* providing an [output function](./doc/OutputFunction.md) (e.g. 
for dense output or for event location) to the solvers
* providing mass- and jacobi-matrices for the solvers (with support for
banded matrices)
* all the solvers' parameters for fine-tuning them, 
see [Options for Solvers](./doc/SolverOptions.md)
* support for problems with "special structure", 
see [special structure](./doc/SpecialStructure.md)

## What are the requirements for this module

This module needs the *compiled* Fortran solvers as shared libraries 
(i.e. `.so`, `.dylib` or `.dll` files, respectively). 
The build script of this module tries to do this compilation. 
It was tested with:

* Linux (64bit) and `gfortran`
* MacOS and `gfortran`
* Windows 7 (64bit) and `gfortran` of mingw-w64 (`x86_64`)

If you want to compile the solvers yourself (perhaps with different
options and/or a different compiler), then just call
`ODEInterface.help_solversupport` for further informations (help topics)
on how to compile the solvers and how to create shared libraries.

## Further help

see `ODEInterface.help_overview` for an overview of some help topics. 

## Examples

* basic examples can be found in [examples/BasicExamples](./examples/BasicExamples/)
* more advanced examples can be found in [examples/AdvancedExamples](./examples/AdvancedExamples/)
* the Reentry Problem of the Apollo11 mission can be found in [examples/ReentryOfApollo11](./examples/ReentryOfApollo11/)

## Contacting the author of this module

The author of this julia module is 

     Dr. Christian Ludwig
     email: ludwig@ma.tum.de
       (Faculty of Mathematics, Technische Universität München)

