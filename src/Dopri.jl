# Functions common for all Dopri-Solvers

"""
  Type encapsulating all required data for Dopri-Solver-Callbacks.
  
  We have the typical calling stack:

       dopri5/dop853
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           ccall( DOPRI5_/DOP853_ ... )
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW1RHSCallback_c                    │  ⎬ cb. rhs
              │    rhs                                    │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_dopriSoloutCallback_c               │  ⎪
              │    call_julia_output_fcn( ... STEP ...)   │  ⎪ cb. solout
              │        output_fcn ( ... STEP ...)         │  ⎬ with eval
              │            eval_sol_fcn                   │  ⎪
              │                ccall(CONTD5_/CONTD8_ ... )│  ⎪
              └───────────────────────────────────────────┘  ⎭
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)
  """
type DopriInternalCallInfos{FInt} <: ODEinternalCallInfos
  callid       :: Array{UInt64}         # the call-id for this info
  logio        :: IO                    # where to log
  loglevel     :: UInt64                # log level
  # RHS:
  rhs          :: Function              # right-hand-side 
  rhs_mode     :: RHS_CALL_MODE         # how to call rhs
  rhs_lprefix  :: AbstractString        # saved log-prefix for rhs
  # SOLOUT & output function
  output_mode  :: OUTPUTFCN_MODE        # what mode for output function
  output_fcn   :: Function              # the output function to call
  output_data  :: Dict                  # extra_data for output_fcn
  eval_sol_fcn :: Function              # eval_sol_fcn 
  tOld         :: Float64               # tOld and
  tNew         :: Float64               # tNew and
  xNew         :: Vector{Float64}       # xNew of current solout interval
  cont_i       :: Vector{FInt}          # argument to contd5/contd8
  cont_s       :: Vector{Float64}       # argument to contd5/contd8
  cont_con     :: Ptr{Float64}          # saved pointers for contd5/contd8
  cont_icomp   :: Ptr{FInt}             # saved pointers for contd5/contd8
  cont_nd      :: Ptr{FInt}             # saved pointers for contd5/contd8
end

"""
       type DopriArguments{FInt} <: AbstractArgumentsODESolver{FInt}
  
  Stores Arguments for Dopri solver.
  
  FInt is the Integer type used for the fortran compilation:
  FInt ∈ (Int32,Int4)
  """
type DopriArguments{FInt} <: AbstractArgumentsODESolver{FInt}
  N       :: Vector{FInt}      # Dimension
  FCN     :: Ptr{Void}         # rhs callback
  t       :: Vector{Float64}   # start time (and current)
  tEnd    :: Vector{Float64}   # end time
  x       :: Vector{Float64}   # initial value (and current state)
  RTOL    :: Vector{Float64}   # relative tolerance
  ATOL    :: Vector{Float64}   # absolute tolerance
  ITOL    :: Vector{FInt}      # switch for RTOL, ATOL
  SOLOUT  :: Ptr{Void}         # solout callback
  IOUT    :: Vector{FInt}      # switch for SOLOUT
  WORK    :: Vector{Float64}   # double working array
  LWORK   :: Vector{FInt}      # length of WORK
  IWORK   :: Vector{FInt}      # integer working array
  LIWORK  :: Vector{FInt}      # length of IWORK
  RPAR    :: Vector{Float64}   # add. double-array
  IPAR    :: Vector{FInt}      # add. integer-array
  IDID    :: Vector{FInt}      # Status code
    ## Allow uninitialized construction
  function DopriArguments()
    return new()
  end
end

"""
       function create_dopri_eval_sol_fcn_closure{FInt}(cid::UInt64, d::FInt,
                    method_contd::Ptr{Void})
  
  generates a eval_sol_fcn for dopri5 and dop853.
  
  Why is a closure needed? We need a function `eval_sol_fcn`
  that calls `CONTD5_` or `CONTD8_` (with `ccall`).
  But `CONTD?_` needs the informations for the current state. This
  informations were saved by `unsafe_dopriSoloutCallback` in the
  `DopriInternalCallInfos`. `eval_sol_fcn` needs to get this informations.
  For finding this "callback informations" the "call id" is needed.
  Here comes `create_dopri_eval_sol_fcn_closure` into play: this function
  takes the "call id" (and some other informations) and generates
  a `eval_sol_fcn` with this data.
  
  Why doesn't `unsafe_dopriSoloutCallback` generate a closure (then
  the current state needs not to be saved in `DopriInternalCallInfos`)?
  Because then every call to `unsafe_dopriSoloutCallback` would have
  generated a closure function. That's a lot of overhead: 1 closure function
  for every solout call. With the strategy above, we have 1 closure function
  per ODE-solver-call, i.e. 1 closure function per ODE.

  For the typical calling sequence, see `DopriInternalCallInfos`.
  """
function create_dopri_eval_sol_fcn_closure{FInt}(cid::UInt64, d::FInt,
             method_contd::Ptr{Void})
  
  function eval_sol_fcn_closure(s::Float64)
    lprefix = string(int2logstr(cid),"eval_sol_fcn_closure: ")
    cbi = get(GlobalCallInfoDict,cid,nothing)
    cbi==nothing && throw(InternalErrorODE(
        string("Cannot find call-id ",int2logstr(cid[1]),
               " in GlobalCallInfoDict")))

    (lio,l)=(cbi.logio,cbi.loglevel)
    l_eval = l & LOG_EVALSOL>0

    l_eval && println(lio,lprefix,"called with s=",s)
    cbi.cont_s[1] = s
    result = Vector{Float64}(d)
    if s == cbi.tNew
      result[:] = cbi.xNew
      l_eval && println(lio,lprefix,"not calling contd because s==tNew")
    else
      for k = 1:d
        cbi.cont_i[1] = k
        result[k] = ccall(method_contd,Float64,
          (Ptr{FInt}, Ptr{Float64}, Ptr{Float64}, Ptr{FInt}, Ptr{FInt},),
          cbi.cont_i,cbi.cont_s, cbi.cont_con, cbi.cont_icomp, cbi.cont_nd)
      end
    end
    
    l_eval && println(lio,lprefix,"contd returned ",result)
    return result
  end
  return eval_sol_fcn_closure
end

"""
       function unsafe_dopriSoloutCallback{FInt}(nr_::Ptr{FInt}, 
         told_::Ptr{Float64}, t_::Ptr{Float64}, x_::Ptr{Float64}, 
         n_::Ptr{FInt}, con_::Ptr{Float64},
         icomp_::Ptr{FInt}, nd_::Ptr{FInt}, rpar_::Ptr{Float64}, 
         ipar_::Ptr{FInt}, irtrn_::Ptr{FInt})
  
  This is the solout given as callback to Fortran-dopri.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.

  This function saves the state informations of the solver in
  `DopriInternalCallInfos`, where they can be found by
  the `eval_sol_fcn`, see `create_dopri_eval_sol_fcn_closure`.
  
  Then the user-supplied `output_fcn` is called (which in turn can use
  `eval_sol_fcn`, to evalutate the solution at intermediate points).
  
  The return value of the `output_fcn` is propagated to 
  `DOPRI5_` or `DOP853_`.
  
  For the typical calling sequence, see `DopriInternalCallInfos`.
  """
function unsafe_dopriSoloutCallback{FInt}(nr_::Ptr{FInt}, 
  told_::Ptr{Float64}, t_::Ptr{Float64}, x_::Ptr{Float64}, 
  n_::Ptr{FInt}, con_::Ptr{Float64},
  icomp_::Ptr{FInt}, nd_::Ptr{FInt}, rpar_::Ptr{Float64}, 
  ipar_::Ptr{FInt}, irtrn_::Ptr{FInt})

  nr = unsafe_load(nr_); told = unsafe_load(told_); t = unsafe_load(t_)
  n = unsafe_load(n_)
  x = pointer_to_array(x_,(n,),false)
  ipar = pointer_to_array(ipar_,(2,),false)
  irtrn = pointer_to_array(irtrn_,(1,),false)
  cid = unpackUInt64FromVector(ipar)
  lprefix = string(int2logstr(cid),"unsafe_dopriSoloutCallback: ")
  cbi = get(GlobalCallInfoDict,cid,nothing)
  cbi==nothing && throw(InternalErrorODE(
      string("Cannot find call-id ",int2logstr(cid[1]),
             " in GlobalCallInfoDict")))
  
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_sol = l & LOG_SOLOUT>0

  l_sol && println(lio,lprefix,"called with nr=",nr," told=",told,
                               " t=",t," x=",x)
  
  cbi.tOld = told; cbi.tNew = t; cbi.xNew = x;
  cbi.cont_con = con_; cbi.cont_icomp = icomp_; cbi.cont_nd = nd_;
  cbi.output_data["nr"] = nr

  ret = call_julia_output_fcn(cid,OUTPUTFCN_CALL_STEP,told,t,x,
                              cbi.eval_sol_fcn)
  if      ret == OUTPUTFCN_RET_STOP
    irtrn[1] = -1
  elseif  ret == OUTPUTFCN_RET_CONTINUE
    irtrn[1] = 0
  elseif  ret == OUTPUTFCN_RET_CONTINUE_XCHANGED
    irtrn[1] = 2
  else
    throw(InternalErrorODE(string("Unkown ret=",ret," of output function")))
  end

  return nothing
end

"""
  `cfunction` pointer for unsafe_dopriSoloutCallback with 64bit integers.
  """
const unsafe_dopriSoloutCallback_c = cfunction(
  unsafe_dopriSoloutCallback, Void, (Ptr{Int64}, 
    Ptr{Float64}, Ptr{Float64},Ptr{Float64}, 
    Ptr{Int64}, Ptr{Float64},
    Ptr{Int64}, Ptr{Int64}, Ptr{Float64}, 
    Ptr{Int64}, Ptr{Int64}));

"""
  `cfunction` pointer for unsafe_dopriSoloutCallback with 32bit integers.
  """
const unsafe_dopriSoloutCallbacki32_c = cfunction(
  unsafe_dopriSoloutCallback, Void, (Ptr{Int32}, 
    Ptr{Float64}, Ptr{Float64},Ptr{Float64}, 
    Ptr{Int32}, Ptr{Float64},
    Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, 
    Ptr{Int32}, Ptr{Int32}));

"""
       function dopri_extract_commonOpt{FInt}(t0::Real, T::Real, x0::Vector, 
                    opt::AbstractOptionsODE, args::DopriArguments{FInt})
             -> (d,nrdense,rhs_mode,output_mode,output_fcn)
  
  calls solver_extract_commonOpt and additionally sets args.ITOL, args.IOUT 
  """
function dopri_extract_commonOpt{FInt}(t0::Real, T::Real, x0::Vector, 
             opt::AbstractOptionsODE, args::DopriArguments{FInt})
  
  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)

  args.ITOL = [ scalarFlag?0:1 ]
  args.IOUT = [ FInt( output_mode == OUTPUTFCN_NEVER? 0:
                     (output_mode == OUTPUTFCN_DENSE?2:1) )]
  
  return (d,nrdense,rhs_mode,output_mode,output_fcn)
end



# vim:syn=julia:cc=79:fdm=indent:
