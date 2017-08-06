# Functions common for all Dopri-Solvers

"""
  Type encapsulating all required data for Dopri-Solver-Callbacks.
  
  We have the typical calling stack:

       dopri5/dop853
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           ccall( DOPRI5_/DOP853_ ... )
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW1RHSCallback                      │  ⎬ cb. rhs
              │    rhs                                    │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_dopriSoloutCallback                 │  ⎪
              │    call_julia_output_fcn( ... STEP ...)   │  ⎪ cb. solout
              │        output_fcn ( ... STEP ...)         │  ⎬ with eval
              │            eval_sol_fcn                   │  ⎪
              │                ccall(CONTD5_/CONTD8_ ... )│  ⎪
              └───────────────────────────────────────────┘  ⎭
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)
  """
mutable struct DopriInternalCallInfos{FInt<:FortranInt, 
        RHS_F<:Function, OUT_F<:Function} <: ODEinternalCallInfos
  logio        :: IO                    # where to log
  loglevel     :: UInt64                # log level
  # RHS:
  rhs          :: RHS_F                 # right-hand-side 
  rhs_mode     :: RHS_CALL_MODE         # how to call rhs
  rhs_lprefix  :: AbstractString        # saved log-prefix for rhs
  # SOLOUT & output function
  output_mode  :: OUTPUTFCN_MODE        # what mode for output function
  output_fcn   :: OUT_F                 # the output function to call
  output_data  :: Dict                  # extra_data for output_fcn
  out_lprefix  :: AbstractString        # saved log-prefix for solout
  eval_sol_fcn :: Function              # eval_sol_fcn (for output_fcn)
  eval_lprefix :: AbstractString        # saved log-prefix for eval_sol
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
       mutable struct DopriArguments{FInt<:FortranInt} <: 
                AbstractArgumentsODESolver{FInt}
  
  Stores Arguments for Dopri solver.
  
  FInt is the Integer type used for the fortran compilation.
  """
mutable struct DopriArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
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
  IPAR    :: Ref{DopriInternalCallInfos} # misuse IPAR
  IDID    :: Vector{FInt}      # Status code
    ## Allow uninitialized construction
  function DopriArguments{FInt}(dummy::FInt) where FInt
    return new{FInt}()
  end
end

"""
       function create_dopri_eval_sol_fcn_closure( cbi::CI, d::FInt, 
               method_contd::Ptr{Void}) where {FInt<:FortranInt, 
                                               CI<:DopriInternalCallInfos}
  
  generates a eval_sol_fcn for dopri5 and dop853.
  
  Why is a closure needed? We need a function `eval_sol_fcn`
  that calls `CONTD5_` or `CONTD8_` (with `ccall`).
  But `CONTD?_` needs the informations for the current state. This
  informations were saved by `unsafe_dopriSoloutCallback` in the
  `DopriInternalCallInfos`. `eval_sol_fcn` needs to get this informations.
  Here comes `create_dopri_eval_sol_fcn_closure` into play: this function
  takes the call-informations and generates a `eval_sol_fcn` with this data.
  
  Why doesn't `unsafe_dopriSoloutCallback` generate a closure (then
  the current state needs not to be saved in `DopriInternalCallInfos`)?
  Because then every call to `unsafe_dopriSoloutCallback` would have
  generated a closure function. That's a lot of overhead: 1 closure function
  for every solout call. With the strategy above, we have 1 closure function
  per ODE-solver-call, i.e. 1 closure function per ODE.

  For the typical calling sequence, see `DopriInternalCallInfos`.
  """
function create_dopri_eval_sol_fcn_closure( cbi::CI, d::FInt, 
        method_contd::Ptr{Void}) where {FInt<:FortranInt, 
                                        CI<:DopriInternalCallInfos}
  
  function eval_sol_fcn_closure(s::Float64)
    (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.eval_lprefix)
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
       function unsafe_dopriSoloutCallback(
               nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
               x_::Ptr{Float64}, n_::Ptr{FInt}, con_::Ptr{Float64},
               icomp_::Ptr{FInt}, nd_::Ptr{FInt}, rpar_::Ptr{Float64}, 
               cbi::CI, irtrn_::Ptr{FInt}) where {FInt<:FortranInt, 
                                                  CI<:DopriInternalCallInfos}
  
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
function unsafe_dopriSoloutCallback(
        nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
        x_::Ptr{Float64}, n_::Ptr{FInt}, con_::Ptr{Float64},
        icomp_::Ptr{FInt}, nd_::Ptr{FInt}, rpar_::Ptr{Float64}, 
        cbi::CI, irtrn_::Ptr{FInt}) where {FInt<:FortranInt, 
                                           CI<:DopriInternalCallInfos}

  nr = unsafe_load(nr_); told = unsafe_load(told_); t = unsafe_load(t_)
  n = unsafe_load(n_)
  x = unsafe_wrap(Array, x_,(n,),false)
  irtrn = unsafe_wrap(Array, irtrn_,(1,),false)
  
  (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.out_lprefix)
  l_sol = l & LOG_SOLOUT>0

  l_sol && println(lio,lprefix,"called with nr=",nr," told=",told,
                               " t=",t," x=",x)
  
  cbi.tOld = told; cbi.tNew = t; cbi.xNew = x;
  cbi.cont_con = con_; cbi.cont_icomp = icomp_; cbi.cont_nd = nd_;
  cbi.output_data["nr"] = nr

  ret = call_julia_output_fcn(cbi,OUTPUTFCN_CALL_STEP,told,t,x,
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
       function unsafe_dopriSoloutCallback_c(cbi::CI, 
               fint_flag::FInt) where {FInt,CI}
          -> C-callable function pointer
  """
function unsafe_dopriSoloutCallback_c(cbi::CI, 
        fint_flag::FInt) where {FInt,CI}
  return cfunction(unsafe_dopriSoloutCallback, Void, (Ptr{FInt}, 
    Ptr{Float64}, Ptr{Float64},Ptr{Float64}, 
    Ptr{FInt}, Ptr{Float64},
    Ptr{FInt}, Ptr{FInt}, Ptr{Float64}, 
    Ref{CI}, Ptr{FInt}))
end

"""
       function dopri_extract_commonOpt(
               t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
               args::DopriArguments{FInt}) where FInt<:FortranInt
             -> (d,nrdense,rhs_mode,output_mode,output_fcn)
  
  calls solver_extract_commonOpt and additionally sets args.ITOL, args.IOUT 
  """
function dopri_extract_commonOpt(
        t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
        args::DopriArguments{FInt}) where FInt<:FortranInt
  
  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)

  args.ITOL = [ scalarFlag ? 0 : 1 ]
  args.IOUT = [ FInt( output_mode == OUTPUTFCN_NEVER ? 0 :
                     (output_mode == OUTPUTFCN_DENSE ? 2 : 1) )]
  
  return (d,nrdense,rhs_mode,output_mode,output_fcn)
end



# vim:syn=julia:cc=79:fdm=indent:
