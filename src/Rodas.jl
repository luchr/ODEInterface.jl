# Functions for ODE-Solver: Rodas

"""Name for Loading rodas solver (64bit integers)."""
const DL_RODAS                = "rodas"

"""Name for Loading rodas solver (32bit integers)."""
const DL_RODAS_I32            = "rodas_i32"

"""macro for import rodas solver."""
macro import_rodas()
  :(
    using ODEInterface: rodas, rodas_i32
  )
end

"""macro for import rodas dynamic lib names."""
macro import_DLrodas()
  :(
    using ODEInterface: DL_RODAS, DL_RODAS_I32
  )
end

"""macro for importing Rodas help."""
macro import_rodas_help()
  :(
    using ODEInterface: help_rodas_compile, help_rodas_license
  )
end

"""
  Type encapsulating all required data for Rodas-Solver-Callbacks.
  
  We have the typical calling stack:

       rodas
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           ccall( RODAS_  ... )
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW2RHSCallback                      │  ⎬ cb. rhs
              │    rhs                                    │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_rodasSoloutCallback:                │  ⎪
              │    call_julia_output_fcn( ... STEP ...)   │  ⎪ cb. solout
              │        output_fcn ( ... STEP ...)         │  ⎬ with eval
              │            eval_sol_fcn                   │  ⎪
              │                ccall(CONTRO_ ... )        │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW1JacCallback:                     │  ⎬ cb. jacobian
              │    call_julia_jac_fcn                     │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HWRhsTimeDerivCallback:             │  ⎬ cb. ∂f/∂t
              │    rhsdt                                  │  ⎪
              └───────────────────────────────────────────┘  ⎭
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)
  """
mutable struct RodasInternalCallInfos{FInt<:FortranInt, RHS_F,
        OUT_F, JAC_F, RHSDT_F} <: ODEinternalCallInfos
  logio        :: IO                    # where to log
  loglevel     :: UInt64                # log level
  # special structure
  M1           :: FInt
  M2           :: FInt
  # RHS
  rhs          :: RHS_F                 # right-hand-side 
  rhs_mode     :: RHS_CALL_MODE         # how to call rhs
  rhs_lprefix  :: AbstractString        # saved log-prefix for rhs
  # RHS time derivative
  rhsdt        :: RHSDT_F               # right-hand side time derivative
  rhsdt_prefix :: AbstractString        # saved log-prefix for rhsdt
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
  cont_i       :: Vector{FInt}          # argument to contro
  cont_s       :: Vector{Float64}       # argument to contro
  cont_cont    :: Ptr{Float64}          # saved pointers for contro
  cont_lrc     :: Ptr{FInt}             # saved pointers for contro
  # Massmatrix
  massmatrix   :: AbstractArray{Float64}# saved mass matrix
  # Jacobimatrix
  jacobimatrix :: JAC_F                 # function for jacobi matrix
  jacobibandstruct                      # (l,u) bandstructure or nothing
  jac_lprefix  :: AbstractString        # saved log-prefix for jac
end

"""
       mutable struct RodasArguments{FInt<:FortranInt} <: 
                AbstractArgumentsODESolver{FInt}
  
  Stores Arguments for Rodas solver.
  
  FInt is the Integer type used for the fortran compilation.
  """
mutable struct RodasArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
  N       :: Vector{FInt}      # Dimension
  FCN     :: Ptr{Cvoid}        # rhs callback
  IFCN    :: Vector{FInt}      # autonomous (0) or not (1)
  t       :: Vector{Float64}   # start time (and current)
  tEnd    :: Vector{Float64}   # end time
  x       :: Vector{Float64}   # initial value (and current state)
  H       :: Vector{Float64}   # initial step size
  RTOL    :: Vector{Float64}   # relative tolerance
  ATOL    :: Vector{Float64}   # absolute tolerance
  ITOL    :: Vector{FInt}      # switch for RTOL, ATOL
  JAC     :: Ptr{Cvoid}        # jacobian callback
  IJAC    :: Vector{FInt}      # jacobian given as callback?
  MLJAC   :: Vector{FInt}      # lower bandwidth of jacobian
  MUJAC   :: Vector{FInt}      # upper bandwidth of jacobian
  DFX     :: Ptr{Cvoid}        # dFCN/dt callback
  IDFX    :: Vector{FInt}      # DFX given as callback?
  MAS     :: Ptr{Cvoid}
  IMAS    :: Vector{FInt}      # mass matrix given as callback?
  MLMAS   :: Vector{FInt}      # lower bandwidth of mass matrix
  MUMAS   :: Vector{FInt}      # upper bandwidth of mass matrix
  SOLOUT  :: Ptr{Cvoid}        # solout callback
  IOUT    :: Vector{FInt}      # switch for SOLOUT
  WORK    :: Vector{Float64}   # double working array
  LWORK   :: Vector{FInt}      # length of WORK
  IWORK   :: Vector{FInt}      # integer working array
  LIWORK  :: Vector{FInt}      # length of IWORK
  RPAR    :: Vector{Float64}   # add. double-array
  IPAR    :: Ref{RodasInternalCallInfos} # misuse IPAR
  IDID    :: Vector{FInt}      # Status code
    ## Allow uninitialized construction
  function RodasArguments{FInt}(dummy::FInt) where FInt
    return new{FInt}()
  end
end


"""
       function unsafe_rodasSoloutCallback(
               nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
               x_::Ptr{Float64}, cont_::Ptr{Float64}, lrc_::Ptr{FInt}, 
               n_::Ptr{FInt}, rpar_::Ptr{Float64}, cbi::CI, 
               irtrn_::Ptr{FInt}) where {FInt<:FortranInt, 
                                         CI<:RodasInternalCallInfos}
  
  This is the solout given as callback to Fortran-rodas.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.

  This function saves the state informations of the solver in
  `RodasInternalCallInfos`, where they can be found by
  the `eval_sol_fcn`, see `create_rodas_eval_sol_fcn_closure`.
  
  Then the user-supplied `output_fcn` is called (which in turn can use
  `eval_sol_fcn`, to evalutate the solution at intermediate points).
  
  The return value of the `output_fcn` is propagated to `RODAS_`.
  
  For the typical calling sequence, see `RodasInternalCallInfos`.
  """
function unsafe_rodasSoloutCallback(
        nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
        x_::Ptr{Float64}, cont_::Ptr{Float64}, lrc_::Ptr{FInt}, 
        n_::Ptr{FInt}, rpar_::Ptr{Float64}, cbi::CI, 
        irtrn_::Ptr{FInt}) where {FInt<:FortranInt, 
                                  CI<:RodasInternalCallInfos}

  nr = unsafe_load(nr_); told = unsafe_load(told_); t = unsafe_load(t_)
  n = unsafe_load(n_)
  x = unsafe_wrap(Array, x_, (n,), own=false)
  irtrn = unsafe_wrap(Array, irtrn_, (1,), own=false)

  (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.out_lprefix)
  l_sol = l & LOG_SOLOUT>0

  l_sol && println(lio,lprefix,"called with nr=",nr," told=",told,
                               " t=",t," x=",x)

  cbi.tOld = told; cbi.tNew = t; cbi.xNew = x
  cbi.cont_cont = cont_; cbi.cont_lrc = lrc_
  cbi.output_data["nr"] = nr

  ret = call_julia_output_fcn(cbi,OUTPUTFCN_CALL_STEP,told,t,x,
                              cbi.eval_sol_fcn)

  if      ret == OUTPUTFCN_RET_STOP
    irtrn[1] = -1
  elseif  ret == OUTPUTFCN_RET_CONTINUE
    irtrn[1] = 0
  elseif  ret == OUTPUTFCN_RET_CONTINUE_XCHANGED
    throw(FeatureNotSupported(string("Sorry; rodas does not support ",
          "to change the solution inside the output function.")))
  else
    throw(InternalErrorODE(string("Unkown ret=",ret," of output function")))
  end

  return nothing
end

"""
       function unsafe_rodasSoloutCallback_c(cbi::CI, 
               fint_flag::FInt) where {FInt,CI}
  """
function unsafe_rodasSoloutCallback_c(cbi::CI, 
        fint_flag::FInt) where {FInt,CI}
  return @cfunction(unsafe_rodasSoloutCallback, Cvoid, (Ptr{FInt},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
    Ptr{Float64}, Ptr{FInt}, 
    Ptr{FInt}, Ptr{Float64}, Ref{CI}, Ptr{FInt}))
end

"""
       function create_rodas_eval_sol_fcn_closure(cbi::CI, d::FInt, 
               method_contro::Ptr{Cvoid}) where {FInt<:FortranInt, 
                                                CI<:RodasInternalCallInfos}
  
  generates a eval_sol_fcn for rodas.
  
  Why is a closure needed? We need a function `eval_sol_fcn`
  that calls `CONTRO_` (with `ccall`).
  But `CONTRO_` needs the informations for the current state. This
  informations were saved by `unsafe_rodasSoloutCallback` in the
  `RodasInternalCallInfos`. `eval_sol_fcn` needs to get this informations.
  Here comes `create_rodas_eval_sol_fcn_closure` into play: this function
  takes the call informations and generates a `eval_sol_fcn` with this data.
  
  Why doesn't `unsafe_rodasSoloutCallback` generate a closure (then
  the current state needs not to be saved in `RodasInternalCallInfos`)?
  Because then every call to `unsafe_rodasSoloutCallback` would have
  generated a closure function. That's a lot of overhead: 1 closure function
  for every solout call. With the strategy above, we have 1 closure function
  per ODE-solver-call, i.e. 1 closure function per ODE.

  For the typical calling sequence, see `RodasInternalCallInfos`.
  """
function create_rodas_eval_sol_fcn_closure(cbi::CI, d::FInt, 
        method_contro::Ptr{Cvoid}) where {FInt<:FortranInt, 
                                         CI<:RodasInternalCallInfos}
  
  function eval_sol_fcn_closure(s::Float64)
    (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.eval_lprefix)
    l_eval = l & LOG_EVALSOL>0

    l_eval && println(lio,lprefix,"called with s=",s)
    cbi.cont_s[1] = s
    result = Vector{Float64}(undef, d)
    if s == cbi.tNew
      result[:] = cbi.xNew
      l_eval && println(lio,lprefix,"not calling contro because s==tNew")
    else
      for k = 1:d
        cbi.cont_i[1] = k
        result[k] = ccall(method_contro,Float64,
          (Ptr{FInt}, Ptr{Float64}, Ptr{Float64}, Ptr{FInt}, ),
          cbi.cont_i,cbi.cont_s, cbi.cont_cont, cbi.cont_lrc)
      end
    end
    
    l_eval && println(lio,lprefix,"contro returned ",result)
    return result
  end
  return eval_sol_fcn_closure
end


"""
        function rodas(rhs, t0::Real, T::Real,
                       x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)
  

  `retcode` can have the following values:

        1: computation successful
        2: computation. successful, but interrupted by output function
       -1: computation unsuccessful


  main call for using Fortran rodas solver.

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
      ║ EPS             │ the rounding unit                        │   1e-16 ║
      ║                 │ 0 < OPT_EPS < 1.0                        │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ METHODCHOICE    │ Choice of coefficients:                  │       1 ║
      ║                 │ 1: Hairer, Wanner: Solving ODE II,       │         ║
      ║                 │    page 452                              │         ║
      ║                 │ 2: same as 1, with different params      │         ║
      ║                 │ 3: G. Steinbach (1993)                   │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTEPS        │ maximal number of allowed steps          │  100000 ║
      ║                 │ OPT_MAXSTEPS > 0                         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSS           │ maximal step size                        │  T - t0 ║
      ║                 │ OPT_MAXSS ≠ 0                            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ INITIALSS       │ initial step size guess                  │    1e-6 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ STEPSIZESTRATEGY│ Switch for step size strategy            │       1 ║
      ║                 │   1: mod. predictive controller          │         ║
      ║                 │      (Gustafsson)                        │         ║
      ║                 │   2: classical step size control         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ OPT_RHO         │ safety factor for step control algorithm │     0.9 ║
      ║                 │ 0.001 < OPT_RHO < 1.0                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SSMINSEL   &    │ parameters for step size selection       │     0.2 ║
      ║ SSMAXSEL        │ The new step size is chosen subject to   │     6.0 ║
      ║                 │ the restriction                          │         ║
      ║                 │ OPT_SSMINSEL ≤ hnew/hold ≤ OPT_SSMAXSEL  │         ║
      ║                 │ OPT_SSMINSEL ≤ 1, OPT_SSMAXSEL ≥ 1       │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ RHSTIMEDERIV    │ A function providing the time derivative │ nothing ║
      ║                 │ ∂f/∂t of the right-hand side or nothing. │         ║
      ║                 │ If the value given is nothing the solver │         ║
      ║                 │ uses finite differences to approximate   │         ║
      ║                 │ ∂f/∂t.                                   │         ║
      ║                 │ Obviously this options is only relevant  │         ║
      ║                 │ for non-autonomous problems.             │         ║
      ║                 │ The function has to be of the form:      │         ║
      ║                 │   function (t,x,dfdt) -> nothing         │         ║
      ║                 │ Even if the problem has special structure│         ║
      ║                 │ (M1>0, see help_specialstructure) x and  │         ║
      ║                 │ dfdt are always vectors with full length,│         ║
      ║                 │ i.e. length(x)==length(dfdt)==length(x0).│         ║
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
      ║ MASSMATRIX      │ the mass matrix of the problem. If not   │ nothing ║
      ║                 │ given (nothing) then the identiy matrix  │         ║
      ║                 │ is used.                                 │         ║
      ║                 │ The size has to be (d-M1)×(d-M1).        │         ║
      ║                 │ It can be an full matrix or a banded     │         ║
      ║                 │ matrix (BandedMatrix).                   │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝
  """
function rodas(rhs, t0::Real, T::Real,
               x0::Vector, opt::AbstractOptionsODE)
  return rodas_impl(rhs,t0,T,x0,opt,RodasArguments{Int64}(Int64(0)))
end

"""
  rodas with 32bit integers, see rodas.
  """
function rodas_i32(rhs, t0::Real, T::Real,
                   x0::Vector, opt::AbstractOptionsODE)
  return rodas_impl(rhs,t0,T,x0,opt,RodasArguments{Int32}(Int32(0)))
end

function rodas_impl(rhs, 
        t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
        args::RodasArguments{FInt}) where FInt<:FortranInt

  (lio,l,l_g,l_solver,lprefix) = solver_start("rodas",rhs,t0,T,x0,opt)

  (method_rodas, method_contro) = getAllMethodPtrs(
     (FInt == Int64) ? DL_RODAS : DL_RODAS_I32 )

  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)

  args.ITOL = [ scalarFlag ? 0 : 1 ]
  args.IOUT = [ FInt( output_mode == OUTPUTFCN_NEVER ? 0 :
                     (output_mode == OUTPUTFCN_DENSE ? 2 : 1) )]

  (M1,M2,NM1) = extractSpecialStructureOpt(d,opt)
  massmatrix = extractMassMatrix(M1,M2,NM1,args,opt)

  (jacobimatrix,jacobibandstruct,jac_lprefix) =
    extractJacobiOpt(d,M1,M2,NM1,args,opt)

  (rhsdt,rhsdt_prefix) = extractRhsTimeDerivOpt(args,opt)

  flag_implct = args.IMAS[1] ≠ 0
  flag_jband = args.MLJAC[1] < NM1

  # WORK memory
  ljac = flag_jband ? 1+args.MLJAC[1]+args.MUJAC[1] : NM1
  lmas = (!flag_implct) ? 0 :
         (args.MLMAS[1] == NM1) ? NM1 : 1+args.MLMAS[1]+args.MUMAS[1]
  le   = flag_jband ? 1+2*args.MLJAC[1]+args.MUJAC[1] : NM1

  args.LWORK = [ (M1==0) ? d*(ljac+lmas+le+14) + 20 :
                           d*(ljac+14) + NM1*(lmas+le) + 20 ]
  args.WORK = zeros(Float64,args.LWORK[1])

  # IWORK memory
  args.LIWORK = [ d+20 ]
  args.IWORK = zeros(FInt,args.LIWORK[1])

  OPT = nothing
  try
    OPT=OPT_RHSAUTONOMOUS;
    rhsautonomous = convert(Bool,getOption(opt,OPT,false))
    args.IFCN = [ FInt(rhsautonomous ? 0 : 1) ]

    OPT=OPT_MAXSTEPS; args.IWORK[1] = convert(FInt,getOption(opt,OPT,100000))
    @assert 0 < args.IWORK[1]

    OPT=OPT_METHODCHOICE; args.IWORK[2] = convert(FInt,getOption(opt,OPT,1))
    @assert 1 ≤ args.IWORK[2] ≤ 3

    OPT = OPT_STEPSIZESTRATEGY; 
    args.IWORK[3] = convert(FInt,getOption(opt,OPT,1))
    @assert args.IWORK[3] ∈ (1,2,)

    args.IWORK[9] = M1 
    args.IWORK[10] = M2

    OPT=OPT_EPS; args.WORK[1]=convert(Float64,getOption(opt,OPT,1e-16))
    @assert 0 < args.WORK[1] < 1.0

    OPT=OPT_MAXSS; args.WORK[2]=convert(Float64,getOption(opt,OPT,T-t0))
    @assert 0 ≠ args.WORK[2]

    OPT=OPT_SSMINSEL; args.WORK[3]=convert(Float64,getOption(opt,OPT,0.2))
    @assert args.WORK[3] ≤ 1
    OPT=OPT_SSMAXSEL; args.WORK[4]=convert(Float64,getOption(opt,OPT,6.0))
    @assert args.WORK[4] ≥ 1

    OPT = OPT_RHO
    args.WORK[5] = convert(Float64,getOption(opt,OPT,0.9))
    @assert 0.001 < args.WORK[5] < 1.0

    OPT = OPT_INITIALSS
    args.H = [ convert(Float64,getOption(opt,OPT,1e-6)) ]
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  args.RPAR = zeros(Float64,0)
  args.IDID = zeros(FInt,1)
  rhs_lprefix = "unsafe_HW2RHSCallback: "
  out_lprefix = "unsafe_rodasSoloutCallback: "
  eval_lprefix = "eval_sol_fcn_closure: "

  cbi = RodasInternalCallInfos(lio,l,M1,M2,rhs,rhs_mode,rhs_lprefix,
      rhsdt === nothing ? dummy_func : rhsdt, rhsdt_prefix,
      output_mode,output_fcn,Dict(),
      out_lprefix,eval_sol_fcn_noeval,eval_lprefix,
      NaN,NaN,Vector{Float64}(),Vector{FInt}(undef, 1),
                                Vector{Float64}(undef, 1),
      Ptr{Float64}(C_NULL),Ptr{FInt}(C_NULL),
      massmatrix === nothing ? zeros(0,0) : massmatrix,
      jacobimatrix === nothing ? dummy_func : jacobimatrix,
      jacobibandstruct,jac_lprefix)

  if output_mode == OUTPUTFCN_DENSE
    cbi.eval_sol_fcn = create_rodas_eval_sol_fcn_closure(cbi,d,method_contro)
  end

  args.FCN = unsafe_HW2RHSCallback_c(cbi, FInt(0))

  args.SOLOUT = output_mode ≠ OUTPUTFCN_NEVER ?
        unsafe_rodasSoloutCallback_c(cbi, FInt(0)) :
        @cfunction(dummy_func, Cvoid, ())
  args.IPAR = cbi
  args.MAS = unsafe_HW1MassCallback_c(cbi, FInt(0))
  args.JAC = unsafe_HW1JacCallback_c(cbi, FInt(0))
  args.DFX = unsafe_HWRhsTimeDerivCallback_c(cbi, FInt(0))

  output_mode ≠ OUTPUTFCN_NEVER &&
    call_julia_output_fcn(cbi,OUTPUTFCN_CALL_INIT,
      args.t[1],args.tEnd[1],args.x,eval_sol_fcn_init) # ignore result

  if l_solver
    println(lio,lprefix,"call Fortran-rodas $method_rodas with")
    dump(lio,args);
  end

  ccall( method_rodas, Cvoid,
    (Ptr{FInt},  Ptr{Cvoid}, Ptr{FInt},          # N=d, Rightsidefunc, IFCN
     Ptr{Float64}, Ptr{Float64}, Ptr{Float64},   # t, x, tEnd
     Ptr{Float64},                               # h
     Ptr{Float64}, Ptr{Float64}, Ptr{FInt},      # RTOL, ATOL, ITOL
     Ptr{Cvoid}, Ptr{FInt}, Ptr{FInt}, Ptr{FInt},# JAC, IJAC, MLJAC, MUJAC
     Ptr{Cvoid}, Ptr{FInt},                      # DFX, IDFX
     Ptr{Cvoid}, Ptr{FInt}, Ptr{FInt}, Ptr{FInt},# MAS, IMAS, MLMAS, MUMAS
     Ptr{Cvoid}, Ptr{FInt},                      # Soloutfunc, IOUT
     Ptr{Float64}, Ptr{FInt},                    # WORK, LWORK
     Ptr{FInt}, Ptr{FInt},                       # IWORK, LIWORK
     Ptr{Float64}, Ref{RodasInternalCallInfos},  # RPAR, IPAR
     Ptr{FInt},                                  # IDID
    ),
    args.N, args.FCN, args.IFCN,
    args.t, args.x, args.tEnd,
    args.H,
    args.RTOL, args.ATOL, args.ITOL, 
    args.JAC, args.IJAC, args.MLJAC, args.MUJAC,
    args.DFX, args.IDFX,
    args.MAS, args.IMAS, args.MLMAS, args.MUMAS,
    args.SOLOUT, args.IOUT,
    args.WORK, args.LWORK,
    args.IWORK, args.LIWORK,
    args.RPAR, args.IPAR, args.IDID,
  )

  if l_solver
    println(lio,lprefix,"Fortran-rodas $method_rodas returned")
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
  ## Compile RODAS

  The julia ODEInterface tries to compile and link the solvers
  automatically at the build-time of this module. The following
  calls need only be done, if one uses a different compiler and/or if
  one wants to change/add some compiler options.

  The Fortran source code can be found at:
  
       http://www.unige.ch/~hairer/software.html 
  
  See `help_rodas_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile RODAS with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack.o dc_lapack.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o lapack.o lapack.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o rodas.o rodas.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o rodas.so 
                rodas.o dc_lapack.o lapack.o
       gfortran -shared -fPIC -o rodas.dylib
                rodas.o dc_lapack.o lapack.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile RODAS with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack.o dc_lapack.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o lapack.o lapack.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o rodas.o rodas.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o rodas.so 
                rodas.o dc_lapack.o lapack.o
  
  ### Using `gfortran` and 32bit integers (Linux and Mac)
  
  Here is an example how to compile RODAS with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o rodas_i32.o rodas.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):
  
       gfortran -shared -fPIC -o rodas_i32.so 
                 rodas_i32.o dc_lapack_i32.o lapack_i32.o
       gfortran -shared -fPIC -o rodas_i32.dylib
                 rodas_i32.o dc_lapack_i32.o lapack_i32.o
  
  ### Using `gfortran` and 32bit integers (Windows)
  
  Here is an example how to compile RODAS with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o rodas_i32.o rodas.f
  
  In order to get create a shared library (from the object file above) use:

       gfortran -shared -o rodas_i32.dll
                 rodas_i32.o dc_lapack_i32.o lapack_i32.o
  
  """
function help_rodas_compile()
  return Docs.doc(help_rodas_compile)
end

function help_rodas_license()
  return Docs.doc(help_rodas_license)
end

@doc(@doc(hw_license),help_rodas_license)

# Add informations about solvers in global solverInfo-array.
push!(solverInfo,
  SolverInfo("rodas",
    "Rosenbrock method of order 4(3)",
    tuple(:OPT_RTOL, :OPT_ATOL, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN, 
          :OPT_M1, :OPT_M2,
          :OPT_RHSAUTONOMOUS,
          :OPT_MASSMATRIX,
          :OPT_JACOBIMATRIX, :OPT_JACOBIBANDSTRUCT,
          :OPT_RHSTIMEDERIV,
          :OPT_MAXSTEPS, :OPT_METHODCHOICE, :OPT_STEPSIZESTRATEGY,
          :OPT_RHO,
          :OPT_EPS, :OPT_MAXSS,
          :OPT_SSMINSEL, :OPT_SSMAXSEL, :OPT_INITIALSS,
          ),
    tuple(
      SolverVariant("rodas_i64",
        "Rodas with 64bit integers",
        DL_RODAS,
        tuple("rodas","contro")),
      SolverVariant("rodas_i32",
        "Rodas with 32bit integers",
        DL_RODAS_I32,
        tuple("rodas","contro")),
    ),
    help_rodas_compile,
    help_rodas_license,
  )
)

# vim:syn=julia:cc=79:fdm=indent:
