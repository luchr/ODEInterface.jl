# Common Functions for calling the output function and for dense output

"""macro for importing enums, functions, etc. for output function call."""
macro import_outputfcn()
  :(
    using ODEInterface: OUTPUTFCN_CALL_REASON,
          OUTPUTFCN_CALL_INIT, OUTPUTFCN_CALL_STEP, OUTPUTFCN_CALL_DONE,
          OUTPUTFCN_MODE,
          OUTPUTFCN_NEVER, OUTPUTFCN_WODENSE, OUTPUTFCN_DENSE,
          OUTPUTFCN_RETURN_VALUE,
          OUTPUTFCN_RET_STOP,
          OUTPUTFCN_RET_CONTINUE, OUTPUTFCN_RET_CONTINUE_XCHANGED
  )
end

@enum(OUTPUTFCN_CALL_REASON,
      OUTPUTFCN_CALL_INIT, OUTPUTFCN_CALL_STEP, OUTPUTFCN_CALL_DONE)

@doc """
  Possible reasons for calling the `OPT_OUTPUTFCN`
  """
OUTPUTFCN_CALL_REASON

@doc """`OPT_OUTPUTFCN` is called for 1st time."""
OUTPUTFCN_CALL_INIT

@doc """
  `OPT_OUTPUTFCN` is called after a successfull
  integration step
  """
OUTPUTFCN_CALL_STEP

@doc """`OPT_OUTPUTFCN` is called for the last time."""
OUTPUTFCN_CALL_DONE

@enum(OUTPUTFCN_MODE,
      OUTPUTFCN_NEVER, OUTPUTFCN_WODENSE, OUTPUTFCN_DENSE)

@doc """Possible mode for calling the `OPT_OUTPUTFCN`.

see `OPT_OUTPUTFCN` and `OPT_OUTPUTMODE`"""
OUTPUTFCN_MODE

@doc """`OPT_OUTPUTFCN` is never called."""
OUTPUTFCN_NEVER

@doc """
  `OPT_OUTPUTFCN` is called after every successfull step,
  but no support for dense output.
  """
OUTPUTFCN_WODENSE

@doc """
  `OPT_OUTPUTFCN` is called after every successfull step
  and dense output is supported.
  """
OUTPUTFCN_DENSE

@enum(OUTPUTFCN_RETURN_VALUE,
      OUTPUTFCN_RET_STOP,
      OUTPUTFCN_RET_CONTINUE, OUTPUTFCN_RET_CONTINUE_XCHANGED)

@doc """Possible return values of an `OPT_OUTPUTFCN`.

see `OPT_OUTPUTFCN` and `OPT_OUTPUTMODE`"""
OUTPUTFCN_RETURN_VALUE

@doc """tell ODE solver to stop."""
OUTPUTFCN_RET_STOP

@doc """tell ODE solver to continue."""
OUTPUTFCN_RET_CONTINUE

@doc """
  tell ODE solver to continue and inform the solver, that the
  `OUTPUTFCN` has altered the numerical solution.
  """
OUTPUTFCN_RET_CONTINUE_XCHANGED

"""
  Output function that does nothing and returns `OUTPUTFCN_RET_CONTINUE`.
  """
function output_fcn_donothing(args...)
  return OUTPUTFCN_RET_CONTINUE
end

"""
  (Dummy-)eval_sol_function, throwing an error, telling the caller,
  that this is the INIT-call of the output function.
  """
function eval_sol_fcn_init(t)
  throw(FunctionCallNotSupported(string(
    "The julia output function was called with ",
    "reason == OUTPUTFCN_CALL_INIT. During such an init-call there is no",
    "numerical solution available. Hence: no evaluation is possible."
  )))
end

"""
  (Dummy-)eval_sol_function, throwing an error, telling the caller,
  that no evaluation is possible.
  """
function eval_sol_fcn_noeval(t)
  throw(FunctionCallNotSupported(string(
    "The evaluation of the numerical solution is not possible. ",
    "OPT_OUTPUTMODE == OUTPUTFCN_DENSE is needed and the solver ",
    "needs to support dense output. That's not the case here."
  )))
end

"""
  (Dummy-)eval_sol_function, throwing an error, telling the caller,
  that this is the DONE-call of the output function.
  """
function eval_sol_fcn_done(t)
  throw(FunctionCallNotSupported(string(
    "The julia output function was called with ",
    "reason == OUTPUTFCN_CALL_DONE. During such a done-call there is no",
    "numerical solution available. Hence: no evaluation is possible."
  )))
end

"""
        function call_julia_output_fcn{CI<:ODEinternalCallInfos}(cbi::CI,
          reason:: OUTPUTFCN_CALL_REASON, told::Float64,t::Float64,
          x::Vector{Float64},eval_sol_fcn::Function)

  calls the julia output function with the given arguments.

  This is more than a simple call, because this function takes care of
  logging, error-checking, etc.
  """
function call_julia_output_fcn(cbi::CI,
  reason:: OUTPUTFCN_CALL_REASON, told::Float64,t::Float64,
  x::Vector{Float64},eval_sol_fcn::Function) where CI<:ODEinternalCallInfos

  lprefix = "call_julia_output_fcn: "

  (lio,l)=(cbi.logio,cbi.loglevel)
  l_out = l & LOG_OUTPUTFCN>0

  if l_out
    println(lio,lprefix,"calling ",cbi.output_fcn," with reason=",reason,
            " told=",told," t=",t," x=",x, " eval_sol_fcn=",eval_sol_fcn,
            " extra_data=",cbi.output_data)
  end
  result = cbi.output_fcn(reason,told,t,x,eval_sol_fcn,cbi.output_data)
  try
    @assert isa(result,OUTPUTFCN_RETURN_VALUE)
  catch e
    throw(OutputErrorODE("Result was not an OUTPUTFCN_RETURN_VALUE",
                         cbi.output_fcn,e))
  end
  l_out && println(lio,lprefix,cbi.output_fcn," returned ",result)
  return result
end

# vim:syn=julia:cc=79:fdm=indent:
