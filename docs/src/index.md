```@meta
CurrentModule = ODEInterface
```

# ODEInterface.jl

`ODEInterface.jl` provides direct Julia wrappers around the supported Fortran
initial-value ODE solvers. `ODEInterfaceDiffEq.jl` adapts these wrappers to the
SciML `solve` interface; use this package directly when configuring the native
solver options and callback protocol.

```@docs
ODEInterface
```

## Initial-value solvers

Each solver accepts `rhs`, an initial and final time, an initial state, and an
`OptionsODE` object. The right-hand side below uses the default return-array
calling convention.

The names below are `public`, not exported, so bring in the ones you need
explicitly (or reach them as `ODEInterface.dopri5`).

```jldoctest
julia> using ODEInterface: OptionsODE, dopri5

julia> rhs(t, u) = u;

julia> t, u, retcode, stats = dopri5(rhs, 0.0, 1.0, [1.0], OptionsODE("example"));

julia> retcode
1
```

```@docs
OptionsODE
dopri5
dop853
odex
seulex
radau
radau5
rodas
ddeabm
ddebdf
```

## Callback protocol

These values configure the direct-wrapper callback protocol. They are primarily
useful to packages that adapt ODEInterface solvers to another solving API.

```@docs
RHS_CALL_INSITU
OUTPUTFCN_CALL_REASON
OUTPUTFCN_CALL_STEP
OUTPUTFCN_WODENSE
OUTPUTFCN_DENSE
OUTPUTFCN_RET_CONTINUE
OUTPUTFCN_RET_CONTINUE_XCHANGED
```
