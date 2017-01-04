# Functions for ODE-Solver: Seulex

"""Name for Loading seulex solver (64bit integers)."""
const DL_SEULEX               = "seulex"

"""Name for Loading seulex solver (32bit integers)."""
const DL_SEULEX_I32           = "seulex_i32"

"""macro for import seulex solver."""
macro import_seulex()
  :(
    using ODEInterface: seulex, seulex_i32
  )
end

"""macro for import seulex dynamic lib names."""
macro import_DLseulex()
  :(
    using ODEInterface: DL_SEULEX, DL_SEULEX_I32
  )
end

"""macro for importing Seulex help."""
macro import_seulex_help()
  :(
    using ODEInterface: help_seulex_compile, help_seulex_license
  )
end

"""
  Type encapsulating all required data for Seulex-Solver-Callbacks.
  
  We have the typical calling stack:

       seulex       
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           ccall( SEULEX_  ... )
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW2RHSCallback                      │  ⎬ cb. rhs
              │    rhs                                    │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_seulexSoloutCallback                │  ⎪
              │    call_julia_output_fcn( ... STEP ...)   │  ⎪ cb. solout
              │        output_fcn ( ... STEP ...)         │  ⎬ with eval
              │            eval_sol_fcn                   │  ⎪
              │                ccall(CONTEX_ ... )        │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW1JacCallback:                     │  ⎬ cb. jacobian
              │    call_julia_jac_fcn(             )      │  ⎪
              └───────────────────────────────────────────┘  ⎭
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)
  """
type SeulexInternalCallInfos{FInt<:FortranInt, RHS_F<:Function,
        OUT_F<:Function, JAC_F<:Function} <: ODEinternalCallInfos
  logio        :: IO                    # where to log
  loglevel     :: UInt64                # log level
  # special structure
  M1           :: FInt
  M2           :: FInt
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
  cont_rc      :: Ptr{Float64}          # saved pointers for contex
  cont_lrc     :: Ptr{FInt}             # saved pointers for contex
  cont_ic      :: Ptr{FInt}             # saved pointers for contex
  cont_lic     :: Ptr{FInt}             # saved pointers for contex
  # Massmatrix
  massmatrix   :: AbstractArray{Float64}# saved mass matrix
  # Jacobimatrix
  jacobimatrix :: JAC_F                 # function for jacobi matrix
  jacobibandstruct                      # Bandstruktur oder nothing
  jac_lprefix  :: AbstractString        # saved log-prefix for jac
end

"""
       type SeulexArguments{FInt<:FortranInt} <: 
                AbstractArgumentsODESolver{FInt}
  
  Stores Arguments for Seulex solver.
  
  FInt is the Integer type used for the fortran compilation.
  """
type SeulexArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
  N       :: Vector{FInt}      # Dimension
  FCN     :: Ptr{Void}         # rhs callback
  IFCN    :: Vector{FInt}      # autonomous (0) or not (1)
  t       :: Vector{Float64}   # start time (and current)
  tEnd    :: Vector{Float64}   # end time
  x       :: Vector{Float64}   # initial value (and current state)
  H       :: Vector{Float64}   # initial step size
  RTOL    :: Vector{Float64}   # relative tolerance
  ATOL    :: Vector{Float64}   # absolute tolerance
  ITOL    :: Vector{FInt}      # switch for RTOL, ATOL
  JAC     :: Ptr{Void}         # jacobian callback
  IJAC    :: Vector{FInt}      # jacobian given as callback?
  MLJAC   :: Vector{FInt}      # lower bandwidth of jacobian
  MUJAC   :: Vector{FInt}      # upper bandwidth of jacobian
  MAS     :: Ptr{Void}
  IMAS    :: Vector{FInt}      # mass matrix given as callback?
  MLMAS   :: Vector{FInt}      # lower bandwidth of mass matrix
  MUMAS   :: Vector{FInt}      # upper bandwidth of mass matrix
  SOLOUT  :: Ptr{Void}         # solout callback
  IOUT    :: Vector{FInt}      # switch for SOLOUT
  WORK    :: Vector{Float64}   # double working array
  LWORK   :: Vector{FInt}      # length of WORK
  IWORK   :: Vector{FInt}      # integer working array
  LIWORK  :: Vector{FInt}      # length of IWORK
  RPAR    :: Vector{Float64}   # add. double-array
  IPAR    :: Ref{SeulexInternalCallInfos} # misuse IPAR
  IDID    :: Vector{FInt}      # Status code
    ## Allow uninitialized construction
  function SeulexArguments()
    return new()
  end
end

"""
        function unsafe_seulexSoloutCallback{FInt<:FortranInt,
                CI<:SeulexInternalCallInfos}(
                nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
                x_::Ptr{Float64}, rc_::Ptr{Float64}, lrc_::Ptr{FInt}, 
                ic_::Ptr{FInt}, lic_::Ptr{FInt}, n_::Ptr{FInt}, 
                rpar_::Ptr{Float64}, cbi::CI, irtrn_::Ptr{FInt})
  
  This is the solout given as callback to Fortran-seulex.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.

  This function saves the state informations of the solver in
  `SeulexInternalCallInfos`, where they can be found by
  the `eval_sol_fcn`, see `create_seulex_eval_sol_fcn_closure`.
  
  Then the user-supplied `output_fcn` is called (which in turn can use
  `eval_sol_fcn`, to evalutate the solution at intermediate points).
  
  The return value of the `output_fcn` is propagated to `SEULEX_`.
  
  For the typical calling sequence, see `SeulexInternalCallInfos`.
  """
function unsafe_seulexSoloutCallback{FInt<:FortranInt,
        CI<:SeulexInternalCallInfos}(
        nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
        x_::Ptr{Float64}, rc_::Ptr{Float64}, lrc_::Ptr{FInt}, 
        ic_::Ptr{FInt}, lic_::Ptr{FInt}, n_::Ptr{FInt}, 
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
  cbi.cont_rc = rc_; cbi.cont_lrc = lrc_; 
  cbi.cont_ic = ic_; cbi.cont_lic = lic_;
  cbi.output_data["nr"] = nr

  ret = call_julia_output_fcn(cbi,OUTPUTFCN_CALL_STEP,told,t,x,
                              cbi.eval_sol_fcn)

  if      ret == OUTPUTFCN_RET_STOP
    irtrn[1] = -1
  elseif  ret == OUTPUTFCN_RET_CONTINUE
    irtrn[1] = 0
  elseif  ret == OUTPUTFCN_RET_CONTINUE_XCHANGED
    throw(FeatureNotSupported(string("Sorry; seulex does not support ",
          "to change the solution inside the output function.")))
  else
    throw(InternalErrorODE(string("Unkown ret=",ret," of output function")))
  end
  
  return nothing
end

"""
        function unsafe_seulexSoloutCallback_c{FInt,CI}
                (cbi::CI, fint_flag::FInt)
  """
function unsafe_seulexSoloutCallback_c{FInt,CI}(cbi::CI, fint_flag::FInt)
  return cfunction(unsafe_seulexSoloutCallback, Void, (Ptr{FInt},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
    Ptr{Float64}, Ptr{FInt}, Ptr{FInt}, Ptr{FInt},
    Ptr{FInt}, Ptr{Float64}, Ref{CI}, Ptr{FInt}))
end

"""
        function create_seulex_eval_sol_fcn_closure{FInt<:FortranInt,
                CI<:SeulexInternalCallInfos}(
                cbi::CI, d::FInt, method_contex::Ptr{Void})
  
  generates a eval_sol_fcn for seulex.
  
  Why is a closure needed? We need a function `eval_sol_fcn`
  that calls `CONTEX_` (with `ccall`).
  But `CONTEX_` needs the informations for the current state. This
  informations were saved by `unsafe_seulexSoloutCallback` in the
  `SeulexInternalCallInfos`. `eval_sol_fcn` needs to get this informations.
  Here comes `create_seulex_eval_sol_fcn_closure` into play: this function
  takes the call informations and generates a `eval_sol_fcn` with this data.
  
  Why doesn't `unsafe_seulexSoloutCallback` generate a closure (then
  the current state needs not to be saved in `SeulexInternalCallInfos`)?
  Because then every call to `unsafe_seulexSoloutCallback` would have
  generated a closure function. That's a lot of overhead: 1 closure function
  for every solout call. With the strategy above, we have 1 closure function
  per ODE-solver-call, i.e. 1 closure function per ODE.

  For the typical calling sequence, see `SeulexInternalCallInfos`.
  """
function create_seulex_eval_sol_fcn_closure{FInt<:FortranInt,
        CI<:SeulexInternalCallInfos}(
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
          cbi.cont_i,cbi.cont_s, cbi.cont_rc, cbi.cont_lrc,
          cbi.cont_ic, cbi.cont_lic)
      end
    end
    
    l_eval && println(lio,lprefix,"contex returned ",result)
    return result
  end
  return eval_sol_fcn_closure
end

"""
       function seulex(rhs::Function, t0::Real, T::Real,
                       x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)
  

  `retcode` can have the following values:

        1: computation successful
        2: computation. successful, but interrupted by output function
       -1: computation unsuccessful


  main call for using Fortran seulex solver.
  
  This solver support problems with special structure, see
  `help_specialstructure`.
  
  In `opt` the following options are used:
  
      ╔═════════════════╤══════════════════════════════════════════╤═════════╗
      ║  Option OPT_…   │ Description                              │ Default ║
      ╠═════════════════╪══════════════════════════════════════════╪═════════╣
      ║ RHSAUTONOMOUS   │ Flag, if right-hand side is autonomous.  │   false ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ M1 & M2         │ parameter for special structure, see     │       0 ║
      ║                 │ above                                    │      M1 ║
      ║                 │ M1, M2 ≥ 0                               │         ║
      ║                 │ M1 +M2 ≤ length(x0)                      │         ║
      ║                 │ (M1==M2==0) || (M1≠0≠M2)                 │         ║
      ║                 │ M1 % M2 == 0 or M1==0                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
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
      ║ LAMBDADENSE     │ parameter λ of dense output              │       0 ║
      ║                 │ OPT_LAMBDADENSE ∈ {0,1}                  │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ EPS             │ the rounding unit                        │   1e-16 ║
      ║                 │ 0 < OPT_EPS < 1.0                        │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ TRANSJTOH       │ The solver transforms the jacobian       │   false ║
      ║                 │ matrix to Hessenberg form.               │         ║
      ║                 │ This option is not supported if the      │         ║
      ║                 │ system is "implicit" (i.e. a mass matrix │         ║
      ║                 │ is given) or if jacobian is banded.      │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTEPS        │ maximal number of allowed steps          │  100000 ║
      ║                 │ OPT_MAXSTEPS > 0                         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSS           │ maximal step size                        │  T - t0 ║
      ║                 │ OPT_MAXSS ≠ 0                            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ INITIALSS       │ initial step size guess                  │    1e-6 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXEXCOLUMN     │ the maximum number of columns in         │      12 ║
      ║                 │ the extrapolation table                  │         ║
      ║                 │ OPT_MAXEXCOLUMN ≥ 3                      │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ STEPSIZESEQUENCE│ switch for the step size sequence        │       2 ║
      ║                 │ 1: 1, 2, 3, 6, 8, 12, 16, 24, 32, 48, …  │         ║
      ║                 │ 2: 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, …  │         ║
      ║                 │ 3: 1, 2, 3, 4, 5,  6,  7,  8,  9, 10, …  │         ║
      ║                 │ 4: 2, 3, 4, 5, 6,  7,  8,  9, 10, 11, …  │         ║
      ║                 │ 1 ≤ OPT_STEPSIZESEQUENCE ≤ 4             │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SSSELECTPAR1 &  │ parameters for step size selection       │     0.1 ║
      ║ SSSELECTPAR2    │ the new step size for the k-th diagonal  │     4.0 ║
      ║                 │ entry is chosen subject to               │         ║
      ║                 │ FMIN/SSSELECTPAR2 ≤ hnewₖ/hold ≤ 1/FMIN  │         ║
      ║                 │ with FMIN = SSSELECTPAR1^(1/(k-1))       │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ ORDERDECFRAC &  │ parameters for the order selection       │     0.7 ║
      ║ ORDERINCFRAC    │ decrease order if                        │     0.9 ║
      ║                 │         W(k-1) ≤   W(k)*ORDERDECFRAC     │         ║
      ║                 │ increase order if                        │         ║
      ║                 │         W(k)   ≤ W(k-1)*ORDERINCFRAC     │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ JACRECOMPFACTOR │ decides whether the jacobian should be   │ min(    ║
      ║                 │ recomputed.                              │   1e-4, ║
      ║                 │ small (≈ 0.001): recompute often         │ RTOL[1])║
      ║                 │ large (≈ 0.1): recompute rarely          │         ║
      ║                 │ i.e. this number represents how costly   │         ║
      ║                 │ Jacobia evaluations are.                 │         ║
      ║                 │ OPT_JACRECOMPFACTOR ≠ 0                  │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ OPT_RHO      &  │ safety factors for step control algorithm│    0.93 ║
      ║ OPT_RHO2        │ hnew=h*RHO*(RHO2*TOL/ERR)^(1/(k-1) )     │    0.80 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MASSMATRIX      │ the mass matrix of the problem. If not   │ nothing ║
      ║                 │ given (nothing) then the identiy matrix  │         ║
      ║                 │ is used.                                 │         ║
      ║                 │ The size has to be (d-M1)×(d-M1).        │         ║
      ║                 │ It can be an full matrix or a banded     │         ║
      ║                 │ matrix (BandedMatrix).                   │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ JACOBIMATRIX    │ A function providing the Jacobian for    │ nothing ║
      ║                 │ ∂f/∂x or nothing. For nothing the solver │         ║
      ║                 │ uses finite differences to calculate the │         ║
      ║                 │ Jacobian.                                │         ║
      ║                 │ The function has to be of the form:      │         ║
      ║                 │   function (t,x,J) -> nothing       (A)  │         ║
      ║                 │ or for M1>0 & JACOBIBANDSTRUCT ≠ nothing │         ║
      ║                 │   function (t,x,J1,…,JK) -> nothing (B)  │         ║
      ║                 │ with K = 1+M1/M2 and (M1+M2==d)          │         ║
      ║                 │ see help_specialstructure                │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ JACOBIBANDSTRUCT│ A tuple (l,u) describing the banded      │ nothing ║
      ║                 │ structure of the Jacobian or nothing if  │         ║
      ║                 │ the Jacobian is full.                    │         ║
      ║                 │ see help_specialstructure                │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ WORKFORRHS      │ estimated works (complexity) for a call  │     1.0 ║
      ║ WORKFORJAC      │ to                                       │     5.0 ║
      ║ WORKFORDEC      │ WORKFORRHS: right-hand side f            │     1.0 ║
      ║ WORKFORSOL      │ WORKFORJAC: JACOBIMATRIX                 │     1.0 ║
      ║                 │ WORKFORDEC: LU-decomposition             │         ║
      ║                 │ WORKFORSOL: Forward- and Backward subst. │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝
  """
function seulex(rhs::Function, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return seulex_impl(rhs,t0,T,x0,opt,SeulexArguments{Int64}())
end

"""
  seulex with 32bit integers, see seulex.
  """
function seulex_i32(rhs::Function, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return seulex_impl(rhs,t0,T,x0,opt,SeulexArguments{Int32}())
end

function seulex_impl{FInt<:FortranInt}(rhs::Function, 
        t0::Real, T::Real, x0::Vector,
        opt::AbstractOptionsODE, args::SeulexArguments{FInt})
  
  (lio,l,l_g,l_solver,lprefix) = solver_start("seulex",rhs,t0,T,x0,opt)
  
  (method_seulex, method_contex) = getAllMethodPtrs(
     (FInt == Int64)? DL_SEULEX : DL_SEULEX_I32 )

  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)
  
  args.ITOL = [ scalarFlag?0:1 ]
  args.IOUT = [ FInt( output_mode == OUTPUTFCN_NEVER? 0:
                     (output_mode == OUTPUTFCN_DENSE?2:1) )]

  (M1,M2,NM1) = extractSpecialStructureOpt(d,opt)
  massmatrix = extractMassMatrix(M1,M2,NM1,args,opt)

  (jacobimatrix,jacobibandstruct,jac_lprefix) =
    extractJacobiOpt(d,M1,M2,NM1,args,opt)

  flag_implct = args.IMAS[1] ≠ 0
  flag_jband = args.MLJAC[1] < NM1

  OPT = nothing; KM = FInt(0)
  try
    OPT=OPT_MAXEXCOLUMN; KM=convert(FInt,getOption(opt,OPT,12))
    @assert KM ≥ 3
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end
  
  # WORK memory
  ljac = flag_jband? 1+args.MLJAC[1]+args.MUJAC[1]: NM1
  lmas = (!flag_implct)? 0 :
         (args.MLMAS[1] == NM1)? NM1 : 1+args.MLMAS[1]+args.MUMAS[1]
  le   = flag_jband? 1+2*args.MLJAC[1]+args.MUJAC[1] : NM1
  KM2  = 2+KM*ceil(FInt,(KM+3)/2)
  
  args.LWORK = [ (M1==0)? d*(ljac+lmas+le+KM+8)+4*KM+20+KM2*nrdense :
                          d*(ljac+KM+8)+NM1*(lmas+le)+4*KM+20+KM2*nrdense ]
  args.WORK = zeros(Float64,args.LWORK[1])
  
  # IWORK memory
  args.LIWORK = [ 2*d+KM+20+nrdense ]
  args.IWORK = zeros(FInt,args.LIWORK[1])
  
  try
    OPT=OPT_RHSAUTONOMOUS;
    rhsautonomous = convert(Bool,getOption(opt,OPT,false))
    args.IFCN = [ rhsautonomous? 0:1 ]

    OPT=OPT_TRANSJTOH; 
    transjtoh  = convert(Bool,getOption(opt,OPT,false))
    @assert (!transjtoh) || 
            ( (!flag_jband) && (!flag_implct ))  string(
            "Does not work, if the jacobian is banded or if the system ",
            " has a mass matrix (≠Id).")
    args.IWORK[1] = transjtoh? 1 : 0

    OPT=OPT_MAXSTEPS; args.IWORK[2] = convert(FInt,getOption(opt,OPT,100000))
    @assert 0 < args.IWORK[2]
    
    args.IWORK[3] = KM

    OPT = OPT_STEPSIZESEQUENCE
    args.IWORK[4] = convert(FInt,getOption(opt,OPT,2))
    @assert 1 ≤ args.IWORK[4] ≤ 4
    
    OPT = OPT_LAMBDADENSE
    args.IWORK[5] = convert(FInt,getOption(opt,OPT,0))
    @assert args.IWORK[5] ∈ (0,1)
    
    args.IWORK[6] = nrdense

    for k in 1:nrdense
      args.IWORK[20+k] = k
    end

    args.IWORK[9] = M1 
    args.IWORK[10] = M2

    OPT=OPT_EPS; args.WORK[1]=convert(Float64,getOption(opt,OPT,1e-16))
    @assert 0 < args.WORK[1] < 1.0

    OPT=OPT_MAXSS; args.WORK[2]=convert(Float64,getOption(opt,OPT,T-t0))
    @assert 0 ≠ args.WORK[2]

    OPT = OPT_JACRECOMPFACTOR
    args.WORK[3] = convert(Float64,getOption(opt,OPT,min(1e-4,args.RTOL[1])))
    @assert !(args.WORK[3]==0)

    OPT = OPT_SSSELECTPAR1; 
    args.WORK[4] = convert(Float64,getOption(opt,OPT,0.1))
    @assert 0 < args.WORK[4]

    OPT = OPT_SSSELECTPAR2
    args.WORK[5] = convert(Float64,getOption(opt,OPT,4.0))
    @assert 0 < args.WORK[5]

    OPT = OPT_ORDERDECFRAC
    args.WORK[6] = convert(Float64,getOption(opt,OPT,0.7))
    @assert 0 < args.WORK[6]

    OPT = OPT_ORDERINCFRAC
    args.WORK[7] = convert(Float64,getOption(opt,OPT,0.9))
    @assert 0 < args.WORK[7]

    OPT = OPT_RHO2
    args.WORK[8] = convert(Float64,getOption(opt,OPT,0.8))
    @assert 0 < args.WORK[8]
    
    OPT = OPT_RHO
    args.WORK[9] = convert(Float64,getOption(opt,OPT,0.93))
    @assert 0 < args.WORK[9]
    
    OPT = OPT_WORKFORRHS
    args.WORK[10] = convert(Float64,getOption(opt,OPT,1.0))
    @assert 0 < args.WORK[10]

    OPT = OPT_WORKFORJAC
    args.WORK[11] = convert(Float64,getOption(opt,OPT,5.0))
    @assert 0 < args.WORK[11]

    OPT = OPT_WORKFORDEC
    args.WORK[12] = convert(Float64,getOption(opt,OPT,1.0))
    @assert 0 < args.WORK[12]

    OPT = OPT_WORKFORSOL
    args.WORK[13] = convert(Float64,getOption(opt,OPT,1.0))
    @assert 0 < args.WORK[13]

    OPT = OPT_INITIALSS
    args.H = [ convert(Float64,getOption(opt,OPT,1e-6)) ]
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  args.RPAR = zeros(Float64,0)
  args.IDID = zeros(FInt,1)
  rhs_lprefix = "unsafe_HW2RHSCallback: "
  out_lprefix = "unsafe_seulexSoloutCallback: "
  eval_lprefix = "eval_sol_fcn_closure: "

  cbi = SeulexInternalCallInfos(lio,l,M1,M2,rhs,rhs_mode,rhs_lprefix,
      output_mode,output_fcn,Dict(),
      out_lprefix,eval_sol_fcn_noeval,eval_lprefix,
      NaN,NaN,Vector{Float64}(),Vector{FInt}(1),Vector{Float64}(1),
      Ptr{Float64}(C_NULL),
      Ptr{FInt}(C_NULL),Ptr{FInt}(C_NULL),Ptr{FInt}(C_NULL),
      massmatrix==nothing?zeros(0,0):massmatrix,
      jacobimatrix==nothing?dummy_func:jacobimatrix,
      jacobibandstruct,jac_lprefix)

  if output_mode == OUTPUTFCN_DENSE
    cbi.eval_sol_fcn =  create_seulex_eval_sol_fcn_closure(
        cbi,d,method_contex)
  end

  args.FCN = unsafe_HW2RHSCallback_c(cbi, FInt(0))
  args.SOLOUT = output_mode ≠ OUTPUTFCN_NEVER?
        unsafe_seulexSoloutCallback_c(cbi, FInt(0)):
        cfunction(dummy_func, Void, () )
  args.IPAR = cbi
  args.MAS = unsafe_HW1MassCallback_c(cbi, FInt(0))
  args.JAC = unsafe_HW1JacCallback_c(cbi, FInt(0))

  output_mode ≠ OUTPUTFCN_NEVER &&
    call_julia_output_fcn(cbi,OUTPUTFCN_CALL_INIT,
      args.t[1],args.tEnd[1],args.x,eval_sol_fcn_init) # ignore result

  if l_solver
    println(lio,lprefix,"call Fortran-seulex $method_seulex with")
    dump(lio,args);
  end

  ccall( method_seulex, Void,
    (Ptr{FInt},  Ptr{Void}, Ptr{FInt},           # N=d, Rightsidefunc, IFCN
     Ptr{Float64}, Ptr{Float64}, Ptr{Float64},   # t, x, tEnd
     Ptr{Float64},                               # h
     Ptr{Float64}, Ptr{Float64}, Ptr{FInt},      # RTOL, ATOL, ITOL
     Ptr{Void}, Ptr{FInt}, Ptr{FInt}, Ptr{FInt}, # JAC, IJAC, MLJAC, MUJAC
     Ptr{Void}, Ptr{FInt}, Ptr{FInt}, Ptr{FInt}, # MAS, IMAS, MLMAS, MUMAS
     Ptr{Void}, Ptr{FInt},                       # Soloutfunc, IOUT
     Ptr{Float64}, Ptr{FInt},                    # WORK, LWORK
     Ptr{FInt}, Ptr{FInt},                       # IWORK, LIWORK
     Ptr{Float64}, Ref{SeulexInternalCallInfos}, # RPAR, IPAR 
     Ptr{FInt},                                  # IDID
    ),
    args.N, args.FCN, args.IFCN,
    args.t, args.x, args.tEnd,
    args.H,
    args.RTOL, args.ATOL, args.ITOL, 
    args.JAC, args.IJAC, args.MLJAC, args.MUJAC,
    args.MAS, args.IMAS, args.MLMAS, args.MUMAS,
    args.SOLOUT, args.IOUT,
    args.WORK, args.LWORK,
    args.IWORK, args.LIWORK,
    args.RPAR, args.IPAR, args.IDID,
  )

  if l_solver
    println(lio,lprefix,"Fortran-seulex $method_seulex returned")
    dump(lio,args);
  end

  output_mode ≠ OUTPUTFCN_NEVER &&
    call_julia_output_fcn(cbi,OUTPUTFCN_CALL_DONE,
      args.t[1],args.tEnd[1],args.x,eval_sol_fcn_done) # ignore result

  l_g && println(lio,lprefix,string("done IDID=",args.IDID[1]))
  stats = Dict{AbstractString,Any}(
    "step_predict"       => args.H[1],
    "no_rhs_calls"       => args.IWORK[14],
    "no_jac_calls"       => args.IWORK[15],
    "no_steps"           => args.IWORK[16],
    "no_steps_accepted"  => args.IWORK[17],
    "no_steps_rejected"  => args.IWORK[18],
    "no_lu_decomp"       => args.IWORK[19],
    "no_fw_bw_subst"     => args.IWORK[20],
          )
  return ( args.t[1], args.x, args.IDID[1], stats)
end

"""  
  ## Compile SEULEX 

  The Fortran source code can be found at:
  
       http://www.unige.ch/~hairer/software.html 
  
  See `help_seulex_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile SEULEX with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack.o dc_lapack.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o lapack.o lapack.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o lapackc.o lapackc.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o seulex.o seulex.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o seulex.so 
                seulex.o dc_lapack.o lapack.o lapackc.o
       gfortran -shared -fPIC -o seulex.dylib
                seulex.o dc_lapack.o lapack.o lapackc.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile SEULEX with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack.o dc_lapack.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o lapack.o lapack.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o lapackc.o lapackc.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o seulex.o seulex.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o seulex.so 
                seulex.o dc_lapack.o lapack.o lapackc.o
  
  ### Using `gfortran` and 32bit integers (Linux and Mac)
  
  Here is an example how to compile SEULEX with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o lapackc_i32.o lapackc.f 
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o seulex_i32.o seulex.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):
  
       gfortran -shared -fPIC -o seulex_i32.so 
                 seulex_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
       gfortran -shared -fPIC -o seulex_i32.dylib
                 seulex_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
  
  ### Using `gfortran` and 32bit integers (Windows)
  
  Here is an example how to compile SEULEX with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o lapackc_i32.o lapackc.f 
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o seulex_i32.o seulex.f
  
  In order to get create a shared library (from the object file above) use:

       gfortran -shared -o seulex_i32.dll
                 seulex_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
  
  """
function help_seulex_compile()
  return Docs.doc(help_seulex_compile)
end

function help_seulex_license()
  return Docs.doc(help_seulex_license)
end

@doc(@doc(hw_license),help_seulex_license)

# Add informations about solvers in global solverInfo-array.
push!(solverInfo,
  SolverInfo("seulex",
    "Extrapolation method based on the linearly implicit Eueler method",
    tuple(:OPT_RTOL, :OPT_ATOL, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN, 
          :OPT_M1, :OPT_M2,
          :OPT_RHSAUTONOMOUS,
          :OPT_MASSMATRIX,
          :OPT_JACOBIMATRIX, :OPT_JACOBIBANDSTRUCT,
          :OPT_MAXEXCOLUMN,
          :OPT_TRANSJTOH, :OPT_MAXSTEPS, :OPT_STEPSIZESEQUENCE,
          :OPT_LAMBDADENSE, :OPT_EPS, :OPT_MAXSS, :OPT_JACRECOMPFACTOR,
          :OPT_SSSELECTPAR1, :OPT_SSSELECTPAR2,
          :OPT_ORDERDECFRAC, :OPT_ORDERINCFRAC,
          :OPT_RHO2, :OPT_RHO,
          :OPT_INITIALSS,
          :OPT_WORKFORRHS, :OPT_WORKFORJAC, :OPT_WORKFORDEC,
          :OPT_WORKFORSOL,
          ),
    tuple(
      SolverVariant("seulex_i64",
        "Seulex with 64bit integers",
        DL_SEULEX,
        tuple("seulex","contex")),
      SolverVariant("seulex_i32",
        "Seulex with 32bit integers",
        DL_SEULEX_I32,
        tuple("seulex","contex")),
    ),
    help_seulex_compile,
    help_seulex_license,
  )
)

# vim:syn=julia:cc=79:fdm=indent:
