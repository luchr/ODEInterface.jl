# Functions for ODE-Solver: Odex

"""Name for Loading odex solver (64bit integers)."""
const DL_ODEX                 = "odex"

"""Name for Loading odex solver (32bit integers)."""
const DL_ODEX_I32           = "odex_i32"

"""macro for import Odex solver."""
macro import_odex()
  :(
    using ODEInterface: odex, odex_i32
  )
end

"""macro for import Odex dynamic lib names."""
macro import_DLodex()
  :(
    using ODEInterface: DL_ODEX, DL_ODEX_I32
  )
end

"""macro for importing Odex help."""
macro import_odex_help()
  :(
    using ODEInterface: help_odex_compile, help_odex_license
  )
end

"""
  Type encapsulating all required data for Odex-Solver-Callbacks.
  
  We have the typical calling stack:

       odex
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           ccall( ODEX_ ... )
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW1RHSCallback                      │  ⎬ cb. rhs
              │    rhs                                    │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_odexSoloutCallback                  │  ⎪
              │    call_julia_output_fcn( ... STEP ...)   │  ⎪ cb. solout
              │        output_fcn ( ... STEP ...)         │  ⎬ with eval
              │            eval_sol_fcn                   │  ⎪
              │                ccall(CONTEX_ ... )        │  ⎪
              └───────────────────────────────────────────┘  ⎭
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)
  """
type OdexInternalCallInfos{FInt<:FortranInt,
       RHS_F<:Function, OUT_F<:Function } <: ODEinternalCallInfos
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
  eval_sol_fcn :: Function              # eval_sol_fcn 
  eval_lprefix :: AbstractString        # saved log-prefix for eval_sol
  tOld         :: Float64               # tOld and
  tNew         :: Float64               # tNew and
  xNew         :: Vector{Float64}       # xNew of current solout interval
  cont_i       :: Vector{FInt}          # argument to contex
  cont_s       :: Vector{Float64}       # argument to contex
  cont_con     :: Ptr{Float64}          # saved pointers for contex
  cont_ncon    :: Ptr{FInt}             # saved pointers for contex
  cont_icomp   :: Ptr{FInt}             # saved pointers for contex
  cont_nd      :: Ptr{FInt}             # saved pointers for contex
end

"""
       type OdexArguments{FInt} <: AbstractArgumentsODESolver{FInt}
  
  Stores Arguments for Odex solver.
  
  FInt is the Integer type used for the fortran compilation.
  """
type OdexArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
  N       :: Vector{FInt}      # Dimension
  FCN     :: Ptr{Void}         # rhs callback
  t       :: Vector{Float64}   # start time (and current)
  tEnd    :: Vector{Float64}   # end time
  x       :: Vector{Float64}   # initial value (and current state)
  H       :: Vector{Float64}   # initial step size
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
  IPAR    :: Ref{OdexInternalCallInfos} # misuse IPAR
  IDID    :: Vector{FInt}      # Status code

    ## Allow uninitialized construction
    (::Type{OdexArguments{T}}){T}() = new{T}()
end

"""
        function unsafe_odexSoloutCallback{FInt<:FortranInt,
                CI<:OdexInternalCallInfos}(
                nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
                x_::Ptr{Float64}, n_::Ptr{FInt}, con_::Ptr{Float64}, 
                ncon_::Ptr{FInt}, icomp_::Ptr{FInt}, nd_::Ptr{FInt}, 
                rpar_::Ptr{Float64}, cbi::CI, irtrn_::Ptr{FInt})
  
  This is the solout given as callback to Fortran-odex.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.

  This function saves the state informations of the solver in
  `OdexInternalCallInfos`, where they can be found by
  the `eval_sol_fcn`, see `create_odex_eval_sol_fcn_closure`.
  
  Then the user-supplied `output_fcn` is called (which in turn can use
  `eval_sol_fcn`, to evalutate the solution at intermediate points).
  
  The return value of the `output_fcn` is propagated to `ODEX_`.
  
  For the typical calling sequence, see `OdexInternalCallInfos`.
  """
function unsafe_odexSoloutCallback{FInt<:FortranInt,
        CI<:OdexInternalCallInfos}(
        nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
        x_::Ptr{Float64}, n_::Ptr{FInt}, con_::Ptr{Float64}, 
        ncon_::Ptr{FInt}, icomp_::Ptr{FInt}, nd_::Ptr{FInt}, 
        rpar_::Ptr{Float64}, cbi::CI, irtrn_::Ptr{FInt})
  
  nr = unsafe_load(nr_); told = unsafe_load(told_); t = unsafe_load(t_)
  n = unsafe_load(n_)
  x = unsafe_wrap(Array, x_,(n,),false)
  irtrn = unsafe_wrap(Array, irtrn_,(1,),false)

  (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.out_lprefix)
  l_sol = l & LOG_SOLOUT>0

  l_sol && println(lio,lprefix,"called with nr=",nr," told=",told,
                               " t=",t," x=",x)
  
  cbi.tOld = told; cbi.tNew = t; cbi.xNew = x;
  cbi.cont_con = con_; cbi.cont_ncon = ncon_; cbi.cont_icomp = icomp_;
  cbi.cont_nd = nd_;
  cbi.output_data["nr"] = nr
  
  ret = call_julia_output_fcn(cbi,OUTPUTFCN_CALL_STEP,told,t,x,
                              cbi.eval_sol_fcn)
  if      ret == OUTPUTFCN_RET_STOP
    irtrn[1] = -1
  elseif  ret == OUTPUTFCN_RET_CONTINUE
    irtrn[1] = 0
  elseif  ret == OUTPUTFCN_RET_CONTINUE_XCHANGED
    throw(FeatureNotSupported(string("Sorry; odex does not support ",
          "to change the solution inside the output function.")))
  else
    throw(InternalErrorODE(string("Unkown ret=",ret," of output function")))
  end
  
  return nothing
end

"""
        function unsafe_odexSoloutCallback_c{FInt,CI}(cbi::CI, 
                fint_flag::FInt)
  """
function unsafe_odexSoloutCallback_c{FInt,CI}(cbi::CI, fint_flag::FInt)
  return cfunction(unsafe_odexSoloutCallback, Void, (Ptr{FInt}, 
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
    Ptr{FInt}, Ptr{Float64}, Ptr{FInt},
    Ptr{FInt}, Ptr{FInt}, Ptr{Float64}, 
    Ref{CI}, Ptr{FInt}))
end

"""
        function create_odex_eval_sol_fcn_closure{FInt<:FortranInt,
                CI<:OdexInternalCallInfos}(
                cbi::CI, d::FInt, method_contex::Ptr{Void})
  
  generates a eval_sol_fcn for odex.
  
  Why is a closure needed? We need a function `eval_sol_fcn`
  that calls `CONTEX_` (with `ccall`).
  But `CONTEX_` needs the informations for the current state. This
  informations were saved by `unsafe_odexSoloutCallback` in the
  `OdexInternalCallInfos`. `eval_sol_fcn` needs to get this informations.
  Here comes `create_odex_eval_sol_fcn_closure` into play: this function
  takes the call informations and generates a `eval_sol_fcn` with this data.
  
  Why doesn't `unsafe_odexSoloutCallback` generate a closure (then
  the current state needs not to be saved in `OdexInternalCallInfos`)?
  Because then every call to `unsafe_odexSoloutCallback` would have
  generated a closure function. That's a lot of overhead: 1 closure function
  for every solout call. With the strategy above, we have 1 closure function
  per ODE-solver-call, i.e. 1 closure function per ODE.

  For the typical calling sequence, see `OdexInternalCallInfos`.
  """
function create_odex_eval_sol_fcn_closure{FInt<:FortranInt,
        CI<:OdexInternalCallInfos}(
        cbi::CI, d::FInt, method_contex::Ptr{Void})
  
  function eval_sol_fcn_closure(s::Float64)
    (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.eval_lprefix)
    l_eval = l & LOG_EVALSOL>0

    l_eval && println(lio,lprefix,"called with s=",s)
    cbi.cont_s[1] = s
    result = Vector{Float64}(d)
    if s == cbi.tNew
      result[:] = cbi.xNew
      l_eval && println(lio,lprefix,"not calling contex because s==tNew")
    else
      for k = 1:d
        cbi.cont_i[1] = k
        result[k] = ccall(method_contex,Float64,
          (Ptr{FInt}, Ptr{Float64}, Ptr{Float64}, Ptr{FInt}, 
           Ptr{FInt}, Ptr{FInt},),
          cbi.cont_i,cbi.cont_s, cbi.cont_con, cbi.cont_ncon,
          cbi.cont_icomp, cbi.cont_nd)
      end
    end
    
    l_eval && println(lio,lprefix,"contex returned ",result)
    return result
  end
  return eval_sol_fcn_closure
end

"""
       function odex(rhs::Function, t0::Real, T::Real, 
                     x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)

  `retcode` can have the following values:

        1: computation successful
        2: computation. successful, but interrupted by output function
       -1: error
  
  main call for using Fortran-odex solver. In `opt` the following
  options are used:
  
      ╔═════════════════╤══════════════════════════════════════════╤═════════╗
      ║  Option OPT_…   │ Description                              │ Default ║
      ╠═════════════════╪══════════════════════════════════════════╪═════════╣
      ║ RTOL         &  │ relative and absolute error tolerances   │    1e-3 ║
      ║ ATOL            │ both scalars or both vectors with the    │    1e-6 ║
      ║                 │ length of length(x0)                     │         ║
      ║                 │ error(xₖ) ≤ OPT_RTOLₖ⋅|xₖ|+OPT_ATOLₖ     │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ OUTPUTFCN       │ output function                          │ nothing ║
      ║                 │ see help_outputfcn                       │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ OUTPUTMODE      │ OUTPUTFCN_NEVER:                         │   NEVER ║
      ║                 │   dont't call OPT_OUTPUTFCN              │         ║
      ║                 │ OUTPUTFCN_WODENSE                        │         ║
      ║                 │   call OPT_OUTPUTFCN, but without        │         ║
      ║                 │   possibility for dense output           │         ║
      ║                 │ OUTPUTFCN_DENSE                          │         ║
      ║                 │   call OPT_OUTPUTFCN with support for    │         ║
      ║                 │   dense output                           │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTEPS        │ maximal number of allowed steps          │   10000 ║
      ║                 │ OPT_MAXSTEPS > 0                         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ EPS             │ the rounding unit                        │ 2.3e-16 ║
      ║                 │ 1e-35 < OPT_EPS < 1.0                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSS           │ maximal step size                        │  T - t0 ║
      ║                 │ OPT_MAXSS ≠ 0                            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ INITIALSS       │ initial step size guess                  │    1e-4 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXEXCOLUMN     │ the maximum number of columns in         │       9 ║
      ║                 │ the extrapolation table                  │         ║
      ║                 │ OPT_MAXEXCOLUMN ≥ 3                      │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ STEPSIZESEQUENCE│ switch for the step size sequence        │       4 ║
      ║                 │ 1: 2, 4,  6,  8, 10, 12, 14, 16, …       │ if      ║
      ║                 │ 2: 2, 4,  8, 12, 16, 20, 24, 28, …       │ OUTPUT- ║
      ║                 │ 3: 2, 4,  6,  8, 12, 16, 24, 32, …       │ MODE == ║
      ║                 │ 4: 2, 6, 10, 14, 18, 22, 26, 30, …       │ DENSE;  ║
      ║                 │ 5: 4, 8, 12, 16, 20, 24, 28, 32, …       │ other-  ║
      ║                 │ 1 ≤ OPT_STEPSIZESEQUENCE ≤ 5             │ wise  1 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTABCHECKS   │ how many times is the stability check    │       1 ║
      ║                 │ activated at most in one line of the     │         ║
      ║                 │ extrapolation table                      │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTABCHECKLINE│ stability check is only activated in     │       1 ║
      ║                 │ the lines 1 to MAXMAXSTABCHECKLINE of    │         ║
      ║                 │ the extrapolation table                  │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ DENSEOUTPUTWOEE │ boolean flag: suppress error estimator   │   false ║
      ║                 │ in dense output                          │         ║
      ║                 │ true is only possible, if                │         ║
      ║                 │      OUTPUTMODE == DENSE                 │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ INTERPOLDEGREE  │ determines the degree of interpolation   │       4 ║
      ║                 │ formula:                                 │         ║
      ║                 │ μ = 2*κ - INTERPOLDEGREE + 1             │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SSREDUCTION     │ step size is reduced by factor if the    │     0.5 ║
      ║                 │ stability check is negative              │         ║
      ║                 │ OPT_EPS < OPT_SSREDUCTION < 1            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SSSELECTPAR1 &  │ parameters for step size selection       │    0.02 ║
      ║ SSSELECTPAR2    │ the new step size for the k-th diagonal  │    4.00 ║
      ║                 │ entry is chosen subject to               │         ║
      ║                 │ FMIN/SSSELECTPAR2 ≤ hnewₖ/hold ≤ 1/FMIN  │         ║
      ║                 │ with FMIN = SSSELECTPAR1^(1/(2*k-1))     │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ ORDERDECFRAC &  │ parameters for the order selection       │     0.8 ║
      ║ ORDERINCFRAC    │ decrease order if                        │     0.9 ║
      ║                 │         W(k-1) ≤   W(k)*ORDERDECFRAC     │         ║
      ║                 │ increase order if                        │         ║
      ║                 │         W(k)   ≤ W(k-1)*ORDERINCFRAC     │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ OPT_RHO      &  │ safety factors for step control algorithm│    0.94 ║
      ║ OPT_RHO2        │ hnew=h*RHO*(RHO2*TOL/ERR)^(1/(k-1) )     │    0.65 ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝ 

  """
function odex(rhs::Function, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return odex_impl(rhs,t0,T,x0,opt,OdexArguments{Int64}())
end

"""
  odex with 32bit integers, see odex.
  """
function odex_i32(rhs::Function, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)

  return odex_impl(rhs,t0,T,x0,opt,OdexArguments{Int32}())
end


"""
        function odex_impl{FInt<:FortranInt}(rhs::Function, 
                t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
                args::OdexArguments{FInt})
  
  implementation of odex for FInt.
  """
function odex_impl{FInt<:FortranInt}(rhs::Function, 
        t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
        args::OdexArguments{FInt})
  
  (lio,l,l_g,l_solver,lprefix) = solver_start("odex",rhs,t0,T,x0,opt)
  
  (method_odex, method_contex) = getAllMethodPtrs(
     (FInt == Int64)? DL_ODEX : DL_ODEX_I32 )

  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)

  args.ITOL = [ scalarFlag?0:1 ]
  args.IOUT = [ FInt( output_mode == OUTPUTFCN_NEVER? 0:
                     (output_mode == OUTPUTFCN_DENSE?2:1) )]
  
  OPT = nothing; KM = FInt(0)
  try
    OPT=OPT_MAXEXCOLUMN; KM=convert(FInt,getOption(opt,OPT,9))
    @assert KM ≥ 3
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  # WORK memory
  args.LWORK = [ d*(KM+5)+5*KM+20+(2*KM*(KM+2)+5)*nrdense ]
  args.WORK = zeros(Float64,args.LWORK[1])

  # IWORK memory
  args.LIWORK = [ 2*KM+21+nrdense ]
  args.IWORK = zeros(FInt,args.LIWORK[1])

  OPT = nothing
  try
    # fill IWORK
    OPT=OPT_MAXSTEPS; args.IWORK[1] = convert(FInt,getOption(opt,OPT,10000))
    @assert 0 < args.IWORK[1]
    args.IWORK[2] = KM
    
    OPT = OPT_STEPSIZESEQUENCE
    args.IWORK[3] = convert(FInt,getOption(opt,OPT,
      output_mode == OUTPUTFCN_DENSE?4:1))
    @assert 1 ≤ args.IWORK[3] ≤ 5
    
    OPT=OPT_MAXSTABCHECKS; args.IWORK[4] = convert(FInt,getOption(opt,OPT,1))
    
    OPT = OPT_MAXSTABCHECKLINE
    args.IWORK[5] = convert(FInt,getOption(opt,OPT,1))
    
    OPT = OPT_DENSEOUTPUTWOEE
    supp_flag = convert(Bool,getOption(opt,OPT,false))
    args.IWORK[6] = FInt( supp_flag?1:0 )
    @assert !supp_flag || (supp_flag && output_mode == OUTPUTFCN_DENSE)
    
    OPT = OPT_INTERPOLDEGREE
    args.IWORK[7] = convert(FInt,getOption(opt,OPT,4))
    @assert 1 ≤ args.IWORK[7] ≤ 6
    
    args.IWORK[8] = nrdense
    
    # fill WORK
    OPT=OPT_EPS; args.WORK[1]=convert(Float64,getOption(opt,OPT,2.3e-16))
    @assert 1e-35 < args.WORK[1] < 1.0

    OPT=OPT_MAXSS; args.WORK[2]=convert(Float64,getOption(opt,OPT,T-t0))
    @assert 0 ≠ args.WORK[2]
    
    OPT=OPT_SSREDUCTION; args.WORK[3]=convert(Float64,getOption(opt,OPT,0.5))
    @assert args.WORK[1] < args.WORK[3] < 1.0
    
    OPT = OPT_SSSELECTPAR1; 
    args.WORK[4] = convert(Float64,getOption(opt,OPT,0.02))
    @assert 0 < args.WORK[4]

    OPT = OPT_SSSELECTPAR2
    args.WORK[5] = convert(Float64,getOption(opt,OPT,4.0))
    @assert 0 < args.WORK[5]

    OPT = OPT_ORDERDECFRAC
    args.WORK[6] = convert(Float64,getOption(opt,OPT,0.8))
    @assert 0 < args.WORK[6]

    OPT = OPT_ORDERINCFRAC
    args.WORK[7] = convert(Float64,getOption(opt,OPT,0.9))
    @assert 0 < args.WORK[7]

    OPT = OPT_RHO2
    args.WORK[8] = convert(Float64,getOption(opt,OPT,0.65))
    @assert 0 < args.WORK[8]
    
    OPT = OPT_RHO
    args.WORK[9] = convert(Float64,getOption(opt,OPT,0.94))
    @assert 0 < args.WORK[9]
    
    # H
    OPT = OPT_INITIALSS
    args.H = [ convert(Float64,getOption(opt,OPT,1e-4)) ]
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end
  
  args.RPAR = zeros(Float64,0)
  args.IDID = zeros(FInt,1)
  rhs_lprefix = "unsafe_HW1RHSCallback: "
  out_lprefix = "unsafe_odexSoloutCallback: "
  eval_lprefix = "eval_sol_fcn_closure: "

  cbi = OdexInternalCallInfos(lio,l,rhs,rhs_mode,rhs_lprefix,
      output_mode,output_fcn,
      Dict(),out_lprefix,eval_sol_fcn_noeval,eval_lprefix,
      NaN,NaN,Vector{Float64}(),
      Vector{FInt}(1),Vector{Float64}(1),
      Ptr{Float64}(C_NULL),Ptr{FInt}(C_NULL),
      Ptr{FInt}(C_NULL),Ptr{FInt}(C_NULL))

  if output_mode == OUTPUTFCN_DENSE
    cbi.eval_sol_fcn = create_odex_eval_sol_fcn_closure(cbi,d,method_contex)
  end

  args.FCN = unsafe_HW1RHSCallback_c(cbi, FInt(0))
  args.SOLOUT = output_mode ≠ OUTPUTFCN_NEVER?
        unsafe_odexSoloutCallback_c(cbi, FInt(0)):
     cfunction(dummy_func, Void, () )
  args.IPAR = cbi

  output_mode ≠ OUTPUTFCN_NEVER &&
    call_julia_output_fcn(cbi,OUTPUTFCN_CALL_INIT,
      args.t[1],args.tEnd[1],args.x,eval_sol_fcn_init) # ignore result

  if l_solver
    println(lio,lprefix,"call Fortran-odex $method_odex with")
    dump(lio,args);
  end

  ccall( method_odex, Void,
    (Ptr{FInt},  Ptr{Void},                    # N=d, Rightsidefunc
     Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, # t, x, tEnd
     Ptr{Float64},                             # h
     Ptr{Float64}, Ptr{Float64}, Ptr{FInt},    # RTOL, ATOL, ITOL
     Ptr{Void}, Ptr{FInt},                     # Soloutfunc, IOUT
     Ptr{Float64}, Ptr{FInt},                  # WORK, LWORK
     Ptr{FInt}, Ptr{FInt},                     # IWORK, LIWORK
     Ptr{Float64}, Ref{OdexInternalCallInfos}, # RPAR, IPAR,
     Ptr{FInt},                                # IDID
    ),
    args.N, args.FCN,
    args.t, args.x, args.tEnd,
    args.H,
    args.RTOL, args.ATOL, args.ITOL, 
    args.SOLOUT, args.IOUT,
    args.WORK, args.LWORK,
    args.IWORK, args.LIWORK,
    args.RPAR, args.IPAR, args.IDID,
  )

  if l_solver
    println(lio,lprefix,"Fortran-odex $method_odex returned")
    dump(lio,args);
  end

  output_mode ≠ OUTPUTFCN_NEVER &&
    call_julia_output_fcn(cbi,OUTPUTFCN_CALL_DONE,
      args.t[1],args.tEnd[1],args.x,eval_sol_fcn_done) # ignore result

  l_g && println(lio,lprefix,string("done IDID=",args.IDID[1]))
  stats = Dict{AbstractString,Any}(
    "step_predict"       => args.H[1],
    "no_rhs_calls"       => args.IWORK[17],
    "no_steps"           => args.IWORK[18],
    "no_steps_accepted"  => args.IWORK[19],
    "no_steps_rejected"  => args.IWORK[20]
          )
  return ( args.t[1], args.x, args.IDID[1], stats)
end

"""  
  ## Compile ODEX

  The Fortran source code can be found at:
  
       http://www.unige.ch/~hairer/software.html 
  
  See `help_odex_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile ODEX with `Float64` reals and
  `Int64` integers with `gfortran`:

       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o odex.o odex.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):
  
       gfortran -shared -fPIC -o odex.so odex.o
       gfortran -shared -fPIC -o odex.dylib odex.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile ODEX with `Float64` reals and
  `Int64` integers with `gfortran`:

       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o odex.o odex.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o odex.dll odex.o
  
  ### Using `gfortran` and 32bit integers (Linux and Mac)
  
  Here is an example how to compile ODEX with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fPIC  
                -fdefault-real-8 -fdefault-double-8 
                -o odex_i32.o   odex.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o odex_i32.so odex_i32.o
       gfortran -shared -fPIC -o odex_i32.dylib odex_i32.o
  
  ### Using `gfortran` and 32bit integers (Windows)
  
  Here is an example how to compile ODEX with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c
                -fdefault-real-8 -fdefault-double-8 
                -o odex_i32.o   odex.f
  
  In order to get create a shared library (from the object file above) use:

       gfortran -shared -o odex_i32.dll odex_i32.o
  
  """
function help_odex_compile()
  return Docs.doc(help_odex_compile)
end

function help_odex_license()
  return Docs.doc(help_odex_license)
end

@doc(@doc(hw_license),help_odex_license)

# Add informations about solver in global solverInfo-array.
push!(solverInfo,
  SolverInfo("odex",
    "GBS Extrapolation-Algorithm based on the explicit midpoint rule",
    tuple(:OPT_RTOL, :OPT_ATOL, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN, 
          :OPT_MAXSTEPS, :OPT_EPS, :OPT_MAXS, :OPT_INITIALSS,
          :OPT_MAXEXCOLUMN, :OPT_STEPSIZESEQUENCE, :OPT_MAXSTABCHECKS,
          :OPT_MAXSTABCHECKLINE, :OPT_DENSEOUTPUTWOEE, :OPT_INTERPOLDEGREE,
          :OPT_SSREDUCTION, :OPT_SSSELECTPAR1, :OPT_SSSELECTPAR2,
          :OPT_ORDERDECFRAC, :OPT_ORDERINCFRAC, :OPT_RHO, :OPT_RHO2,
          ),
    tuple(
      SolverVariant("odex_i64",
        "Odex with 64bit integers",
        DL_ODEX,
        tuple("odex", "contex")),
      SolverVariant("odex_i32",
        "Odex with 32bit integers",
        DL_ODEX_I32,
        tuple("odex", "contex")),
    ),
    help_odex_compile,
    help_odex_license,
  )
)


# vim:syn=julia:cc=79:fdm=indent:
