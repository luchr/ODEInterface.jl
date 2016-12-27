[ This file was auto-generated from the module's documentation included in the doc-strings. Use julia's help system to get these informations in a nicer output format. ]

## `OPT_OUTPUTMODE`

This option determines if the `OPT_OUTPUTFCN` is called, and if dense output (the `eval_sol_fcn`) is prepared/supported.

  * `OUTPUTFCN_NEVER`: don't call the output function
  * `OUTPUTFCN_WODENSE`: call the output function, but `eval_sol_fcn`

is not used

  * `OUTPUTFCN_DENSE`: call the output function and prepare `eval_sol_fcn`

## `OPT_OUTPUTFCN`

```
 function outputfcn(reason::OUTPUTFCN_CALL_REASON,
  told::Float64,t::Float64, x::Vector{Float64},eval_sol_fcn::Function,
  extra_data::Dict)  -> OUTPUTFCN_RETURN_VALUE
```

A (julia) function that is called 

1. at beginning of the solution process with

`reason == OUTPUTFCN_CALL_INIT`, `told=t0`, `t`=`T`, `x=x0`, `eval_sol_fcn` a dummy function throwing an error if called, `extra_data` a `Dict` persistent until the last call of the output  function. The return value is *ignored*.

1. after every successfull integration step with

`reason == OUTPUTFCN_CALL_STEP`, `[told,t]` the time interval of the last step, `x` the numerical solution at time `t`, `eval_sol_fcn` a function to evaluate the solution in `t1 ∊ [told,t]`, if requested by `OPT_OUTPUTMODE`, otherwise a call to this function will result in an error.

1. at the end (after the last step) with

`reason == OUTPUTFCN_CALL_DONE`. The return value is *ignored*.

With `eval_sol_fcn`

```
    function eval_sol_fcn(t1::Float64) -> Vector{Float64}
```

the numerical solution can be evaluted in the time interval `[told,t]` (if `OPT_OUTPUTMODE == OUTPUTFCN_DENSE`).

If supported by the solver, the numerical solution may be changed in the `outputfcn` (if `reason == OUTPUTFCN_CALL_STEP`) and the solver continues the process with the changed solution. The return value `OUTPUTFCN_RET_CONTINUE_XCHANGED` indicates  this. `OUTPUTFCN_RET_CONTINUE` tells the solver to continue (without changes in `x`) and `OUTPUTFCN_RET_STOP` tells the solver to stop the solver.



