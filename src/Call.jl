# User-friendly method for calling solvers

"""macro for importing odecall."""
macro import_odecall()
  :(
    using ODEInterface:  odecall
  )
end

"""
     function odecall(solver::Function, rhs::Function, t::Vector, x0::Vector, 
                      opt::AbstractOptionsODE)
          -> (tVec,xVec,retcode,stats)
  
  Calls `solver` with the given right-hand side `rhs`.
  There are two cases:

  1. `2==length(t)`
  1. `2<length(t)`

  If `2==length(t)`, then in the output `tVec` consists of the time points the
  (adaptive) solver has automatically chosen. And the `xVec` has the
  states at this times. So: `tVec` is a `Vector{Float64}(m)` and
  `xVec` is a `Array{Float64}(m,length(x0))`.
  
  If `2<length(t)`, then the values in `t` must be strictly ascending
  or strictly descending. Then a special output function is used to get
  the numerical solution at the given `t`-values. In this case
  `tVec` is a `Vector{Float64}(length(t))` and
  `xVec` is a `Array{Float64}(length(t),length(x0))`.
  
  If in `opt` a output function is given, then this output function is
  also called/used.
  """
function odecall(solver::Function, rhs::Function, t::Vector, x0::Vector, 
                 opt::AbstractOptionsODE)
  tLen = length(t)
  tLen < 2 && throw(ArgumentErrorODE(
              "t must have 2 or more elements; found $tLen",:t))
  try
    t = getVectorCheckLength(t,Real,tLen)
  catch e
    throw(ArgumentErrorODE("Cannot convert to to Vector{Real}",:t,e))
  end

  t0 = t[1];  T = t[end]
  t0==T && throw(ArgumentErrorODE("t[1]==t[end]",:t))
  d = length(x0)
  locopt = OptionsODE(opt)
  orig_outputfcn = getOption(opt, OPT_OUTPUTFCN, nothing)
  orig_outputmode = getOption(opt, OPT_OUTPUTMODE, OUTPUTFCN_NEVER)
  
  if 2==tLen
    tVec = Vector{Float64}()    # we do not know how many steps/points
    xVec = Vector{Float64}()    # the solver will use => tVec and xVec grow

    #             orig    │       local
    #          ───────────┼────────────────
    #             NEVER   │       WODENSE
    #            WODENSE  │       WODENSE
    #             DENSE   │       DENSE
    orig_outputmode != OUTPUTFCN_DENSE && 
      setOption!(locopt,OPT_OUTPUTMODE,OUTPUTFCN_WODENSE)

    # save only grid points of solver in tVec, xVec
    function outputfcn_save(reason::OUTPUTFCN_CALL_REASON,
          told::Float64, tnew::Float64, x::Vector{Float64},
          eval_sol_fcn::Function, extra_data::Dict)

      if reason == OUTPUTFCN_CALL_STEP
        push!(tVec,tnew); append!(xVec,x)
      end
      if orig_outputmode == OUTPUTFCN_NEVER || orig_outputfcn == nothing
        return OUTPUTFCN_RET_CONTINUE
      else
        return orig_outputfcn(reason,told,tnew,x,eval_sol_fcn,extra_data)
      end
    end
    setOption!(locopt,OPT_OUTPUTFCN, outputfcn_save)
  else
    s = sign(T-t0)   # direction: forward or backward in time
    if any( x -> sign(x) != s, diff(t) )
      throw(ArgumentErrorODE(string("Because sign(T-t0)=sign($T-$t0)=$s ",
        "the vector t has to be ",(T≥t0?"ascending":"descending"),
        ". But this is not the case."),:t))
    end
    # We know the size of tVec and xVec (if the solver completes)
    # If not (e.g. stopped by output function) => init with NaN
    tVec = ones(tLen)*NaN; xVec = ones(tLen*d)*NaN;
    setOption!(locopt,OPT_OUTPUTMODE,OUTPUTFCN_DENSE)
    tPos = 1  # next position to check in output fcn

    function outputfcn_givent(reason::OUTPUTFCN_CALL_REASON,
          told::Float64, tnew::Float64, x::Vector{Float64},
          eval_sol_fcn::Function, extra_data::Dict)
      
      if (reason == OUTPUTFCN_CALL_STEP) 
        while (tPos ≤ tLen) &&
              ( (told ≤ t[tPos] ≤ tnew) || (told ≥ t[tPos] ≥ tnew) )
          tVec[tPos] = t[tPos]
          xVec[1+(tPos-1)*d : tPos*d] = eval_sol_fcn(t[tPos])
          tPos += 1
        end
      end
      if orig_outputmode == OUTPUTFCN_NEVER || orig_outputfcn == nothing
        return OUTPUTFCN_RET_CONTINUE
      else
        return orig_outputfcn(reason,told,t,x,eval_sol_fcn,extra_data)
      end
    end
    setOption!(locopt,OPT_OUTPUTFCN, outputfcn_givent)
  end
  
  
  (tout,xout,retcode,stats) = solver( rhs, t0, T, x0, locopt)
  
  return (tVec,reshape(xVec,d,length(xVec)÷d).',retcode,stats)

end

# vim:syn=julia:cc=79:fdm=indent:
