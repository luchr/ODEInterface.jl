[ This file was auto-generated from the module's documentation included in the doc-strings. Use julia's help system to get these informations in a nicer output format. ]

## Calling the Solvers

There are two ways of calling the solvers.

1. A calling convention close to the original Fortran-call, trying to provide/expose all the features the Fortran-codes have.
2. A simplified version, closer to odecalls like in MATLAB.

### The full-featured calling-method

All ODE-solvers have the same calling convention:

```
(t,x,retcode,stats) = 
    odesolver(rhs, t0::Real, T::Real,
              x0::Vector, opt::AbstractOptionsODErhs)

function rhs(t::Float64,x::Vector{Float64}) -> Vector{Float64}
    if OPT_RHS_CALLMODE == RHS_CALL_RETURNS_ARRAY

function rhs(t::Float64,x::Vector{Float64},dx::Vector{Float64}) 
    if OPT_RHS_CALLMODE == RHS_CALL_INSITU
```

The input arguments are:

1. a julia function `rhs` for evaluating the right-hand side of the ODE, see below. It's OK to return something, that `convert` can transform to a `Vector{Float64}`.
2. the initial time `t0`. `(t0,x0)` is the initial value of the  initial value problem.
3. the final time `T`.
4. the initial state `x0`. `(t0,x0)` is the initial value of the  initial value problem.
5. further parameters/options in `opt` for the solver and for the interface.  There is a separate section for the explanation of the options, see help_options.

The output arguments are:

1. `t` the *last* time for which the solution has been computed  (if the whole computation was successfull, then `t==T`)
2. `x` the numerical solation at time `t`
3. `retcode` the return code of the solver (interpretation is solver dependent)

There are two possible ways to provide the Julia right-hand side:

```
function rhs(t::Float64,x::Vector{Float64}) -> Vector{Float64}
```

This is used, if `OPT_RHS_CALLMODE == RHS_CALL_RETURNS_ARRAY`. So you can use anonymous functions like `(t,x) -> x` as right-hand sides. But this form has a price: every time the right-hand side is called, a temporary Array is created (the result). The form

```
function rhs(t::Float64,x::Vector{Float64},dx::Vector{Float64}) 
             -> nothing
```

used if `OPT_RHS_CALLMODE == RHS_CALL_INSITU` does not have this problem. But the right-hand side must be a function filling in the values of `x'` in `dx`.

### The simplified version

There is a simplified calling convention (using the methods above) to provide a method like odecalls in MATLAB,  see `odecall`.



```
function odecall(solver, rhs, t::Vector, x0::Vector,
                opt::AbstractOptionsODE)
    -> (tVec,xVec,retcode,stats)
```

Calls `solver` with the given right-hand side `rhs`. There are two cases:

1. `2==length(t)`
2. `2<length(t)`

If `2==length(t)`, then in the output `tVec` consists of the time points the (adaptive) solver has automatically chosen. And the `xVec` has the states at this times. So: `tVec` is a `Vector{Float64}(m)` and `xVec` is a `Array{Float64}(m,length(x0))`.

If `2<length(t)`, then the values in `t` must be strictly ascending or strictly descending. Then a special output function is used to get the numerical solution at the given `t`-values. In this case `tVec` is a `Vector{Float64}(length(t))` and `xVec` is a `Array{Float64}(length(t),length(x0))`.

If in `opt` a output function is given, then this output function is also called/used.



