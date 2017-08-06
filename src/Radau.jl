# Functions for ODE-Solver: Radau5 and Radau

"""Name for Loading radau5 solver (64bit integers)."""
const DL_RADAU5               = "radau5"

"""Name for Loading radau5 solver (32bit integers)."""
const DL_RADAU5_I32           = "radau5_i32"

"""Name for Loading radau solver (64bit integers)."""
const DL_RADAU                = "radau"

"""Name for Loading radau solver (32bit integers)."""
const DL_RADAU_I32            = "radau_i32"

"""macro for import Radau5 solver."""
macro import_radau5()
  :(
    using ODEInterface: radau5, radau5_i32
  )
end

"""macro for import Radau solver."""
macro import_radau()
  :(
    using ODEInterface: radau, radau_i32
  )
end

"""macro for import Radau5 dynamic lib names."""
macro import_DLradau5()
  :(
    using ODEInterface: DL_RADAU5, DL_RADAU5_I32
  )
end

"""macro for import Radau dynamic lib names."""
macro import_DLradau()
  :(
    using ODEInterface: DL_RADAU, DL_RADAU_I32
  )
end

"""macro for importing Radau5 help."""
macro import_radau5_help()
  :(
    using ODEInterface: help_radau5_compile, help_radau5_license
  )
end

"""macro for importing Radau help."""
macro import_radau_help()
  :(
    using ODEInterface: help_radau_compile, help_radau_license
  )
end

"""
  Type encapsulating all required data for Radau5/Radau-Solver-Callbacks.
  
  We have the typical calling stack:

       radau5/radau
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           ccall( RADAU5_/RADAU_ ... )
               unsafe_HW1MassCallback  
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW2RHSCallback                      │  ⎬ cb. rhs
              │    rhs                                    │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_radauSoloutCallback                 │  ⎪
              │    call_julia_output_fcn( ... STEP ...)   │  ⎪ cb. solout
              │        output_fcn ( ... STEP ...)         │  ⎬ with eval
              │            eval_sol_fcn                   │  ⎪
              │                ccall(CONTR5_/CONTRA_ ... )│  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_HW1JacCallback:                     │  ⎬ cb. jacobian
              │    call_julia_jac_fcn(             )      │  ⎪
              └───────────────────────────────────────────┘  ⎭
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)
  """
type RadauInternalCallInfos{FInt<:FortranInt, RHS_F, 
        OUT_F, JAC_F} <: ODEinternalCallInfos
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
  cont_cont    :: Ptr{Float64}          # saved pointers for contr5
  cont_lrc     :: Ptr{FInt}             # saved pointers for contr5
  # Massmatrix
  massmatrix   :: AbstractArray{Float64}# saved mass matrix
  # Jacobimatrix
  jacobimatrix :: JAC_F                 # function for jacobi matrix
  jacobibandstruct                      # (l,u) bandstructure or nothing
  jac_lprefix  :: AbstractString        # saved log-prefix for jac
end

"""
       type RadauArguments{FInt<:FortranInt} <: 
                AbstractArgumentsODESolver{FInt}
  
  Stores Arguments for Radau5 and Radau solver.
  
  FInt is the Integer type used for the fortran compilation.
  """
type RadauArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
  N       :: Vector{FInt}      # Dimension
  FCN     :: Ptr{Void}         # rhs callback
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
  IPAR    :: Ref{RadauInternalCallInfos} # misuse IPAR
  IDID    :: Vector{FInt}      # Status code
    ## Allow uninitialized construction
  @WHEREFUNC(FInt,
  function RadauArguments{FInt}(dummy::FInt)
    return new{FInt}()
  end
  )
end


"""
        function unsafe_radauSoloutCallback{FInt<:FortranInt}(
                nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
                x_::Ptr{Float64}, cont_::Ptr{Float64}, lrc_::Ptr{FInt}, 
                n_::Ptr{FInt}, rpar_::Ptr{Float64}, 
                ipar_::Ptr{FInt}, irtrn_::Ptr{FInt})
  
  This is the solout given as callback to Fortran radau5/radau.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.

  This function saves the state informations of the solver in
  `RadauInternalCallInfos`, where they can be found by
  the `eval_sol_fcn`, see `create_radau_eval_sol_fcn_closure`.
  
  Then the user-supplied `output_fcn` is called (which in turn can use
  `eval_sol_fcn`, to evalutate the solution at intermediate points).
  
  The return value of the `output_fcn` is propagated to `RADAU5_`/`RADAU_`.
  
  For the typical calling sequence, see `RadauInternalCallInfos`.
  """
function unsafe_radauSoloutCallback{FInt<:FortranInt,
        CI<:RadauInternalCallInfos}(
        nr_::Ptr{FInt}, told_::Ptr{Float64}, t_::Ptr{Float64}, 
        x_::Ptr{Float64}, cont_::Ptr{Float64}, lrc_::Ptr{FInt}, 
        n_::Ptr{FInt}, rpar_::Ptr{Float64}, 
        cbi::CI, irtrn_::Ptr{FInt})
  
  nr = unsafe_load(nr_); told = unsafe_load(told_); t = unsafe_load(t_);
  n = unsafe_load(n_)
  x = unsafe_wrap(Array, x_,(n,),false)
  irtrn = unsafe_wrap(Array, irtrn_,(1,),false)

  (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.out_lprefix)
  l_sol = l & LOG_SOLOUT>0

  l_sol && println(lio,lprefix,"called with nr=",nr," told=",told,
                               " t=",t," x=",x)
  
  cbi.tOld = told; cbi.tNew = t; cbi.xNew = x;
  cbi.cont_cont = cont_; cbi.cont_lrc = lrc_;
  cbi.output_data["nr"] = nr
  
  ret = call_julia_output_fcn(cbi,OUTPUTFCN_CALL_STEP,told,t,x,
                              cbi.eval_sol_fcn)

  if      ret == OUTPUTFCN_RET_STOP
    irtrn[1] = -1
  elseif  ret == OUTPUTFCN_RET_CONTINUE
    irtrn[1] = 0
  elseif  ret == OUTPUTFCN_RET_CONTINUE_XCHANGED
    throw(FeatureNotSupported(string("Sorry; radau5/radau does not support ",
          "to change the solution inside the output function.")))
  else
    throw(InternalErrorODE(string("Unkown ret=",ret," of output function")))
  end
  
  return nothing
end

"""
        function unsafe_radauSoloutCallback_c{FInt,CI}(cbi::CI, 
                fint_flag::FInt)
  """
function unsafe_radauSoloutCallback_c{FInt,CI}(cbi::CI, fint_flag::FInt)
  return cfunction(unsafe_radauSoloutCallback, Void, (Ptr{FInt},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    Ptr{Float64}, Ptr{FInt}, Ptr{FInt},
    Ptr{Float64}, Ref{CI}, Ptr{FInt}))
end

"""
        function create_radau_eval_sol_fcn_closure{FInt<:FortranInt,
                CI<:RadauInternalCallInfos}(
                cbi::CI, d::FInt, method_cont::Ptr{Void})
  
  generates a eval_sol_fcn for radau and radau5.
  
  Why is a closure needed? We need a function `eval_sol_fcn`
  that calls `CONTR5_` OR `CONTRA_` (with `ccall`).
  But `CONTR5_`/`CONTRA_` need the informations for the current state. This
  informations were saved by `unsafe_radauSoloutCallback` in the
  `RadauInternalCallInfos`. `eval_sol_fcn` needs to get this informations.
  Here comes `create_radau_eval_sol_fcn_closure` into play: this function
  takes call informations and generates a `eval_sol_fcn` with this data.
  
  Why doesn't `unsafe_radauSoloutCallback` generate a closure (then
  the current state needs not to be saved in `RadauInternalCallInfos`)?
  Because then every call to `unsafe_radauSoloutCallback` would have
  generated a closure function. That's a lot of overhead: 1 closure function
  for every solout call. With the strategy above, we have 1 closure function
  per ODE-solver-call, i.e. 1 closure function per ODE.

  For the typical calling sequence, see `RadauInternalCallInfos`.
  """
function create_radau_eval_sol_fcn_closure{FInt<:FortranInt,
        CI<:RadauInternalCallInfos}(
        cbi::CI, d::FInt, method_cont::Ptr{Void})
  
  function eval_sol_fcn_closure(s::Float64)
    (lio,l,lprefix)=(cbi.logio,cbi.loglevel,cbi.eval_lprefix)
    l_eval = l & LOG_EVALSOL>0

    l_eval && println(lio,lprefix,"called with s=",s)
    cbi.cont_s[1] = s
    result = Vector{Float64}(d)
    if s == cbi.tNew
      result[:] = cbi.xNew
      l_eval && println(lio,lprefix,"not calling cont because s==tNew")
    else
      for k = 1:d
        cbi.cont_i[1] = k
        result[k] = ccall(method_cont,Float64,
          (Ptr{FInt}, Ptr{Float64}, Ptr{Float64}, Ptr{FInt},),
          cbi.cont_i,cbi.cont_s, cbi.cont_cont, cbi.cont_lrc,
          )
      end
    end
    
    l_eval && println(lio,lprefix,"contr5 returned ",result)
    return result
  end
  return eval_sol_fcn_closure
end


"""
  extracts options specific to radau5 and to radau.
  Fills in `args`: `IWORK[1,2,4,5,6,7,8]`, `WORK[1,2,3,5,6]`,
  `RPAR`, `IDID`, `FCN`, `SOLOUT`
  """
function extractCommonRadauOpt{FInt<:FortranInt}(
        d::FInt,T,t0, args::RadauArguments,opt::AbstractOptionsODE)
  OPT = nothing
  try
    # fill IWORK
    OPT=OPT_TRANSJTOH; 
    transjtoh  = convert(Bool,getOption(opt,OPT,false))
    @assert (!transjtoh) || 
            ( (!flag_jband) && (!flag_implct ))  string(
            "Does not work, if the jacobian is banded or if the system ",
            " has a mass matrix (≠Id).")
    args.IWORK[1] = transjtoh? 1 : 0

    OPT=OPT_MAXSTEPS; args.IWORK[2] = convert(FInt,getOption(opt,OPT,100000))
    @assert 0 < args.IWORK[2]

    OPT=OPT_NEWTONSTARTZERO; zflag = convert(Bool,getOption(opt,OPT,false))
    args.IWORK[4] = zflag? 1 : 0
    
    OPT=OPT_DIMOFIND1VAR; args.IWORK[5] = convert(FInt,getOption(opt,OPT,d))
    @assert 0 < args.IWORK[5]
    OPT=OPT_DIMOFIND2VAR; args.IWORK[6] = convert(FInt,getOption(opt,OPT,0))
    @assert 0 ≤ args.IWORK[6]
    OPT=OPT_DIMOFIND3VAR; args.IWORK[7] = convert(FInt,getOption(opt,OPT,0))
    @assert 0 ≤ args.IWORK[7]

    OPT = string(OPT_DIMOFIND1VAR,", ",OPT_DIMOFIND2VAR,", ",
                 OPT_DIMOFIND3VAR)
    @assert(d == args.IWORK[5]+args.IWORK[6]+args.IWORK[7],string(
      "Sum of dim1, dim2 and dim3 variables must be the dimension of the ",
      "system."))
    
    OPT = OPT_STEPSIZESTRATEGY; 
    args.IWORK[8] = convert(FInt,getOption(opt,OPT,1))
    @assert args.IWORK[8] ∈ (1,2,)

    # fill WORK
    OPT=OPT_EPS; args.WORK[1]=convert(Float64,getOption(opt,OPT,1e-16))
    @assert 1e-19 < args.WORK[1] < 1.0

    OPT = OPT_RHO
    args.WORK[2] = convert(Float64,getOption(opt,OPT,0.9))
    @assert 0.001 < args.WORK[2] < 1.0
    
    OPT = OPT_JACRECOMPFACTOR
    args.WORK[3] = convert(Float64,getOption(opt,OPT,0.001))
    @assert !(args.WORK[3]==0)

    OPT = OPT_FREEZESSLEFT
    args.WORK[5] = convert(Float64,getOption(opt,OPT,1.0))
    @assert args.WORK[5] ≤ 1.0

    OPT = OPT_FREEZESSRIGHT
    args.WORK[6] = convert(Float64,getOption(opt,OPT,1.2))
    @assert args.WORK[6] ≥ 1.0
    
    OPT=OPT_MAXSS; args.WORK[7]=convert(Float64,getOption(opt,OPT,T-t0))
    @assert 0 ≠ args.WORK[7]
    
    OPT=OPT_SSMINSEL; args.WORK[8]=convert(Float64,getOption(opt,OPT,0.2))
    @assert args.WORK[8] ≤ 1
    OPT=OPT_SSMAXSEL; args.WORK[9]=convert(Float64,getOption(opt,OPT,8.0))
    @assert args.WORK[9] ≥ 1

    # H
    OPT = OPT_INITIALSS
    args.H = [ convert(Float64,getOption(opt,OPT,1e-6)) ]
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end
  args.RPAR = zeros(Float64,0)
  args.IDID = zeros(FInt,1)

  rhs_lprefix = "unsafe_HW2RHSCallback: "
  out_lprefix = "unsafe_radauSoloutCallback: "
  eval_lprefix = "eval_sol_fcn_closure: "
  
  return (rhs_lprefix,out_lprefix,eval_lprefix)
end

"""
  calls the radau5 or radau solver after all solver arguments are prepared.
  """
function doRadauSolverCall{FInt<:FortranInt}(
        lio,l,l_g,l_solver,lprefix, d::FInt,M1::FInt,M2::FInt, 
        rhs,rhs_mode,rhs_lprefix, output_mode,output_fcn,out_lprefix,
        eval_lprefix,massmatrix, jacobimatrix,jacobibandstruct,
        jac_lprefix,args,method_solver,method_cont)

  cbi = RadauInternalCallInfos(lio,l,M1,M2,rhs,rhs_mode,rhs_lprefix,
      output_mode,output_fcn,
      Dict(),out_lprefix,eval_sol_fcn_noeval,eval_lprefix,
      NaN,NaN,Vector{Float64}(),
      Vector{FInt}(1),Vector{Float64}(1),
      Ptr{Float64}(C_NULL),Ptr{FInt}(C_NULL),
      massmatrix==nothing?zeros(0,0):massmatrix,
      jacobimatrix==nothing?dummy_func:jacobimatrix,
      jacobibandstruct,jac_lprefix
      )

  if output_mode == OUTPUTFCN_DENSE
    cbi.eval_sol_fcn = create_radau_eval_sol_fcn_closure(cbi,d,method_cont)
  end

  args.FCN = unsafe_HW2RHSCallback_c(cbi, FInt(0))
  args.SOLOUT = output_mode ≠ OUTPUTFCN_NEVER?
        unsafe_radauSoloutCallback_c(cbi, FInt(0)):
        cfunction(dummy_func, Void, () )
  args.IPAR = cbi
  args.MAS = unsafe_HW1MassCallback_c(cbi, FInt(0))
  args.JAC = unsafe_HW1JacCallback_c(cbi, FInt(0))

  output_mode ≠ OUTPUTFCN_NEVER &&
    call_julia_output_fcn(cbi,OUTPUTFCN_CALL_INIT,
      args.t[1],args.tEnd[1],args.x,eval_sol_fcn_init) # ignore result

  if l_solver
    println(lio,lprefix,"call Fortran radau5/radau $method_solver with")
    dump(lio,args);
  end

  ccall( method_solver, Void,
    (Ptr{FInt},  Ptr{Void},                     # N=d, Rightsidefunc
     Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # t, x, tEnd
     Ptr{Float64},                              # h
     Ptr{Float64}, Ptr{Float64}, Ptr{FInt},     # RTOL, ATOL, ITOL
     Ptr{Void}, Ptr{FInt}, Ptr{FInt}, Ptr{FInt},# JAC, IJAC, MLJAC, MUJAC
     Ptr{Void}, Ptr{FInt}, Ptr{FInt}, Ptr{FInt},# MAS, IJAC, MLJAC, MUJAC
     Ptr{Void}, Ptr{FInt},                      # Soloutfunc, IOUT
     Ptr{Float64}, Ptr{FInt},                   # WORK, LWORK
     Ptr{FInt}, Ptr{FInt},                      # IWORK, LIWORK
     Ptr{Float64}, Ref{RadauInternalCallInfos}, # RPAR, IPAR
     Ptr{FInt},                                 # IDID
    ),
    args.N, args.FCN,
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
    println(lio,lprefix,"Fortran radau5/radau $method_solver returned")
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

# Same documentation than for radau, is copied after
# the definition of the radau function.
function radau5(rhs, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return radau5_impl(rhs,t0,T,x0,opt,RadauArguments{Int64}(Int64(0)))
end

"""
  radau5 with 32bit integers, see radau5.
  """
function radau5_i32(rhs, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return radau5_impl(rhs,t0,T,x0,opt,RadauArguments{Int32}(Int32(0)))
end

"""
        function radau5_impl{FInt<:FortranInt}(rhs, 
                t0::Real, T::Real, x0::Vector,
                opt::AbstractOptionsODE, args::RadauArguments{FInt})
  
  implementation of radau5 for FInt.
  """
function radau5_impl{FInt<:FortranInt}(rhs, 
        t0::Real, T::Real, x0::Vector,
        opt::AbstractOptionsODE, args::RadauArguments{FInt})

  (lio,l,l_g,l_solver,lprefix) = solver_start("radau5",rhs,t0,T,x0,opt)
  
  (method_radau5, method_contr5) = getAllMethodPtrs(
     (FInt == Int64)? DL_RADAU5 : DL_RADAU5_I32 )
  
  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)
  
  args.ITOL = [ scalarFlag?0:1 ]
  args.IOUT = [ output_mode == OUTPUTFCN_NEVER? 0: 1 ]

  (M1,M2,NM1) = extractSpecialStructureOpt(d,opt)
  massmatrix = extractMassMatrix(M1,M2,NM1,args,opt)
  
  (jacobimatrix,jacobibandstruct,jac_lprefix) =
    extractJacobiOpt(d,M1,M2,NM1,args,opt)
  
  flag_implct = args.IMAS[1] ≠ 0
  flag_jband = args.MLJAC[1] < NM1
  
  # WORK memory
  ljac = flag_jband? 1+args.MLJAC[1]+args.MUJAC[1]: NM1
  lmas = (!flag_implct)? 0 :
         (args.MLMAS[1] == NM1)? NM1 : 1+args.MLMAS[1]+args.MUMAS[1]
  le   = flag_jband? 1+2*args.MLJAC[1]+args.MUJAC[1] : NM1

  args.LWORK = [ (M1==0)? d*(ljac+lmas+3*le+12)+20 :
                        d*(ljac+12)+NM1*(lmas+3*le)+20 ]
  args.WORK = zeros(Float64,args.LWORK[1])

  # IWORK memoery
  args.LIWORK = [ 3*d + 20 ]
  args.IWORK = zeros(FInt,args.LIWORK[1])

  args.IWORK[9] = M1; args.IWORK[10] = M2
  (rhs_lprefix,out_lprefix,eval_lprefix) = 
    extractCommonRadauOpt(d,T,t0,args,opt)

  OPT = nothing
  try
    OPT=OPT_MAXNEWTONITER; args.IWORK[3] = convert(FInt,getOption(opt,OPT,7))
    @assert 0 < args.IWORK[3]

    OPT = OPT_NEWTONSTOPCRIT
    args.WORK[4] = convert(Float64,getOption(opt,OPT, 
                           max(10*args.WORK[1]/args.RTOL[1],
                               min(0.03,sqrt(args.RTOL[1])))))
    @assert args.WORK[4]>args.WORK[1]/args.RTOL[1]
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  return doRadauSolverCall(lio,l,l_g,l_solver,lprefix,d,M1,M2,
         rhs,rhs_mode,rhs_lprefix,output_mode,output_fcn,
         out_lprefix,eval_lprefix,massmatrix,
         jacobimatrix,jacobibandstruct,jac_lprefix,args,
         method_radau5,method_contr5)
end

"""
       function radau(rhs, t0::Real, T::Real,
                       x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)
       
       function radau5(rhs, t0::Real, T::Real,
                       x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)
  
  `retcode` can have the following values:

        1: computation successful
        2: computation. successful, but interrupted by output function
       -1: input is not consistent
       -2: larger OPT_MAXSTEPS is needed
       -3: step size becomes too small
       -4: matrix is repeatedly singular
  
  main call for using Fortran radau or radau5 solver.
  
  This solver support problems with special structure, see
  `help_specialstructure`.
  
  Remark:
  Because radau and radau5 are collocation methods, there is no difference 
  in the computational costs for OUTPUTFCN_WODENSE and OUTPUTFCN_DENSE.
  
  In `opt` the following options are used:
  
      ╔═════════════════╤══════════════════════════════════════════╤═════════╗
      ║  Option OPT_…   │ Description                              │ Default ║
      ╠═════════════════╪══════════════════════════════════════════╪═════════╣
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
      ║                 │ 1e-19 < OPT_EPS < 1.0                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ TRANSJTOH       │ The solver transforms the jacobian       │   false ║
      ║                 │ matrix to Hessenberg form.               │         ║
      ║                 │ This option is not supported if the      │         ║
      ║                 │ system is "implicit" (i.e. a mass matrix │         ║
      ║                 │ is given) or if jacobian is banded.      │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXNEWTONITER   │ maximum number of Newton iterations for  │       7 ║
      ║                 │ the solution of the implicit system in   │         ║
      ║                 │ each step.                               │         ║
      ║                 │ for radau: MAXNEWTONITER + (NS-3)*2.5    │         ║
      ║                 │   where NS is number of current stages   │         ║
      ║                 │ for radau5:     OPT_MAXNEWTONITER > 0    │         ║
      ║                 │ for radau: 50 > OPT_MAXNEWTONITER > 0    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ NEWTONSTARTZERO │ if `false`, the extrapolated collocation │   false ║
      ║                 │ solution is taken as starting vector for │         ║
      ║                 │ Newton's method. If `true` zero starting │         ║
      ║                 │ values are used. The latter is           │         ║
      ║                 │ recommended if Newton's method has       │         ║
      ║                 │ difficulties with convergence.           │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ NEWTONSTOPCRIT  │ only for radau5:                         │ see     ║
      ║                 │ Stopping criterion for Newton's method.  │   left  ║
      ║                 │ Smaller values make the code slower, but │         ║
      ║                 │ safer.                                   │         ║
      ║                 │ Default:                                 │         ║
      ║                 │  max(10*OPT_EPS/OPT_RTOL[1],             │         ║
      ║                 │       min(0.03,sqrt(OPT_RTOL[1])))       │         ║
      ║                 │ OPT_NEWTONSTOPCRIT > OPT_EPS/OPT_RTOL[1] │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ DIMOFIND1VAR  & │ For differential-algebraic systems of    │ len(x0) ║
      ║ DIMOFIND2VAR  & │ index > 1. The right-hand side should be │       0 ║
      ║ DIMOFIND3VAR    │ written such that the index 1,2,3        │       0 ║
      ║                 │ variables appear in this order.          │         ║
      ║                 │ DIMOFINDzVAR: number of index z vars.    │         ║
      ║                 │ ∑ DIMOFINDzVAR == length(x0)             │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTEPS        │ maximal number of allowed steps          │  100000 ║
      ║                 │ OPT_MAXSTEPS > 0                         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSS           │ maximal step size                        │  T - t0 ║
      ║                 │ OPT_MAXSS ≠ 0                            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ INITIALSS       │ initial step size guess                  │    1e-6 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MINSTAGES     & │ only for radau:                          │       3 ║
      ║ MAXSTAGES       │ minimal and maximal number of stages.    │       7 ║
      ║                 │ The order is given by: 2⋅stages-1        │         ║
      ║                 │ MINSTAGES,MAXSTAGES ∈ (1,3,5,7)          │         ║
      ║                 │ MINSTAGES ≤ MAXSTAGES                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ INITSTAGES      │ only for radau:                          │MINSTAGES║
      ║                 │ number of stages to start with.          │         ║
      ║                 │ MINSTAGES ≤ INITSTAGES ≤ MAXSTAGES       │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ STEPSIZESTRATEGY│ Switch for step size strategy            │       1 ║
      ║                 │   1: mod. predictive controller          │         ║
      ║                 │      (Gustafsson)                        │         ║
      ║                 │   2: classical step size control         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ OPT_RHO         │ safety factor for step control algorithm │     0.9 ║
      ║                 │ 0.001 < OPT_RHO < 1.0                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ JACRECOMPFACTOR │ decides whether the jacobian should be   │   0.001 ║
      ║                 │ recomputed.                              │         ║
      ║                 │ <0: recompute after every accepted step  │         ║
      ║                 │ small (≈ 0.001): recompute often         │         ║
      ║                 │ large (≈ 0.1): recompute rarely          │         ║
      ║                 │ i.e. this number represents how costly   │         ║
      ║                 │ Jacobia evaluations are.                 │         ║
      ║                 │ OPT_JACRECOMPFACTOR ≠ 0                  │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ FREEZESSLEFT  & │ Step size freezing: If                   │     1.0 ║
      ║ FREEZESSRIGHT   │ FREEZESSLEFT < hnew/hold < FREEZESSRIGHT │     1.2 ║
      ║                 │ then the step size is not changed. This  │         ║
      ║                 │ saves, together with a large             │         ║
      ║                 │ JACRECOMPFACTOR, LU-decompositions and   │         ║
      ║                 │ computing time for large systems.        │         ║
      ║                 │ small systems:                           │         ║
      ║                 │    FREEZESSLEFT  ≈ 1.0                   │         ║
      ║                 │    FREEZESSRIGHT ≈ 1.2                   │         ║
      ║                 │ large full systems:                      │         ║
      ║                 │    FREEZESSLEFT  ≈ 0.99                  │         ║
      ║                 │    FREEZESSRIGHT ≈ 2.0                   │         ║
      ║                 │                                          │         ║
      ║                 │ OPT_FREEZESSLEFT  ≤ 1.0                  │         ║
      ║                 │ OPT_FREEZESSRIGHT ≥ 1.0                  │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SSMINSEL   &    │ parameters for step size selection       │     0.2 ║
      ║ SSMAXSEL        │ The new step size is chosen subject to   │     8.0 ║
      ║                 │ the restriction                          │         ║
      ║                 │ OPT_SSMINSEL ≤ hnew/hold ≤ OPT_SSMAXSEL  │         ║
      ║                 │ OPT_SSMINSEL ≤ 1, OPT_SSMAXSEL ≥ 1       │         ║
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
      ║ ORDERDECFACTOR &│ only for radau:                          │     0.8 ║
      ║ ORDERINCFACTOR  │ Order is decreased, if the contractivity │   0.002 ║
      ║                 │ factor is smaller than ORDERDECFACTOR.   │         ║
      ║                 │ Order is increased, if the contractivity │         ║
      ║                 │ factor is larger than ORDERINCFACTOR.    │         ║
      ║                 │ ORDERDECFACTOR > ORDERINCFACTOR > 0      │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ ORDERDECSTEPFAC1│ only for radau:                          │     1.2 ║
      ║ ORDERDECSTEPFAC2│ the order is only decreased if the       │     0.8 ║
      ║                 │ stepsize ratio satisfies                 │         ║
      ║                 │  ORDERDECSTEPFAC2 ≤ hnew/hold ≤          │         ║
      ║                 │               ORDERDECSTEPFAC1           │         ║
      ║                 │ 0 < ORDERDECSTEPFAC2 < ORDERDECSTEPFAC1  │         ║
      ║                 │                                          │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝
  """
function radau(rhs, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return radau_impl(rhs,t0,T,x0,opt,RadauArguments{Int64}(Int64(0)))
end

@doc(@doc(radau),radau5)

"""
  radau with 32bit integers, see radau.
  """
function radau_i32(rhs, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return radau_impl(rhs,t0,T,x0,opt,RadauArguments{Int32}(Int32(0)))
end

"""
        function radau_impl{FInt<:FortranInt}(rhs, 
                t0::Real, T::Real, x0::Vector,
                opt::AbstractOptionsODE, args::RadauArguments{FInt})
  
  implementation of radau for FInt.
  """
function radau_impl{FInt<:FortranInt}(rhs, 
        t0::Real, T::Real, x0::Vector,
        opt::AbstractOptionsODE, args::RadauArguments{FInt})

  (lio,l,l_g,l_solver,lprefix) = solver_start("radau5",rhs,t0,T,x0,opt)
  
  (method_radau, method_contra) = getAllMethodPtrs(
     (FInt == Int64)? DL_RADAU : DL_RADAU_I32 )
  
  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)
  
  args.ITOL = [ scalarFlag?0:1 ]
  args.IOUT = [ output_mode == OUTPUTFCN_NEVER? 0: 1 ]

  (M1,M2,NM1) = extractSpecialStructureOpt(d,opt)
  massmatrix = extractMassMatrix(M1,M2,NM1,args,opt)
  
  (jacobimatrix,jacobibandstruct,jac_lprefix) =
    extractJacobiOpt(d,M1,M2,NM1,args,opt)
  
  flag_implct = args.IMAS[1] ≠ 0
  flag_jband = args.MLJAC[1] < NM1

  OPT = nothing
  NSMAX = 0
  try
    OPT = OPT_MAXSTAGES; NSMAX = convert(FInt,getOption(opt,OPT,7))
    @assert NSMAX ∈ (1,3,5,7)
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  # WORK memory
  ljac = flag_jband? 1+args.MLJAC[1]+args.MUJAC[1]: NM1
  lmas = (!flag_implct)? 0 :
         (args.MLMAS[1] == NM1)? NM1 : 1+args.MLMAS[1]+args.MUMAS[1]
  le   = flag_jband? 1+2*args.MLJAC[1]+args.MUJAC[1] : NM1
  
  args.LWORK = [ (M1==0)? d*(ljac+lmas+NSMAX*le+3*NSMAX+3)+20 :
                          d*(ljac+3*NSMAX+3)+NM1*(lmas+NSMAX*le)+20 ]
  args.WORK = zeros(Float64,args.LWORK[1])
   
  # IWORK memoery
  args.LIWORK = [ (2+(NSMAX-1)÷2)*d+20  ]
  args.IWORK = zeros(FInt,args.LIWORK[1])

  args.IWORK[9] = M1; args.IWORK[10] = M2
  (rhs_lprefix,out_lprefix,eval_lprefix) = 
    extractCommonRadauOpt(d,T,t0,args,opt)

  try
    OPT=OPT_MAXNEWTONITER; args.IWORK[3] = convert(FInt,getOption(opt,OPT,7))
    @assert 0 < args.IWORK[3] < 50
    
    OPT = OPT_MINSTAGES; NSMIN = convert(FInt,getOption(opt,OPT,3))
    @assert NSMIN ∈ (1,3,5,7) && NSMIN ≤ NSMAX
    args.IWORK[11] = NSMIN
    
    OPT = OPT_INITSTAGES; NSINIT = convert(FInt,getOption(opt,OPT,NSMIN))
    @assert NSINIT ∈ (1,3,5,7) && NSMIN ≤ NSINIT ≤ NSMAX
    args.IWORK[13] = NSINIT
    
    OPT = OPT_ORDERINCFACTOR
    orderincfactor = convert(Float64,getOption(opt,OPT,0.002))
    OPT = OPT_ORDERDECFACTOR
    orderdecfactor = convert(Float64,getOption(opt,OPT,0.8))
    OPT = string(OPT_ORDERDECFACTOR," & ",OPT_ORDERINCFACTOR)
    @assert 0 < orderincfactor < orderdecfactor
    args.WORK[10] = orderincfactor
    args.WORK[11] = orderdecfactor
    
    OPT = OPT_ORDERDECSTEPFAC1
    orderdecstepfac1 = convert(Float64,getOption(opt,OPT,1.2))
    OPT = OPT_ORDERDECSTEPFAC2
    orderdecstepfac2 = convert(Float64,getOption(opt,OPT,0.8))
    OPT = string(OPT_ORDERDECSTEPFAC1," & ",OPT_ORDERDECSTEPFAC2)
    @assert 0 < orderdecstepfac2 < orderdecstepfac1
    args.WORK[12] = orderdecstepfac1 
    args.WORK[13] = orderdecstepfac2
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end
  
  return doRadauSolverCall(lio,l,l_g,l_solver,lprefix,d,M1,M2,
         rhs,rhs_mode,rhs_lprefix,output_mode,output_fcn,
         out_lprefix,eval_lprefix,massmatrix,
         jacobimatrix,jacobibandstruct,jac_lprefix,args,
         method_radau,method_contra)
end

"""  
  ## Compile RADAU5 

  The julia ODEInterface tries to compile and link the solvers
  automatically at the build-time of this module. The following
  calls need only be done, if one uses a different compiler and/or if
  one wants to change/add some compiler options.

  The Fortran source code can be found at:
  
       http://www.unige.ch/~hairer/software.html 
  
  See `help_radau5_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile RADAU5 with `Float64` reals and
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
                -o radau5.o radau5.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o radau5.so 
                radau5.o dc_lapack.o lapack.o lapackc.o
       gfortran -shared -fPIC -o radau5.dylib
                radau5.o dc_lapack.o lapack.o lapackc.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile RADAU5 with `Float64` reals and
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
                -o radau5.o radau5.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o radau5.so 
                radau5.o dc_lapack.o lapack.o lapackc.o
  
  ### Using `gfortran` and 32bit integers (Linux and Mac)
  
  Here is an example how to compile RADAU5 with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o lapackc_i32.o lapackc.f 
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o radau5_i32.o radau5.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):
  
       gfortran -shared -fPIC -o radau5_i32.so 
                 radau5_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
       gfortran -shared -fPIC -o radau5_i32.dylib
                 radau5_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
  
  ### Using `gfortran` and 32bit integers (Windows)
  
  Here is an example how to compile RADAU5 with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o lapackc_i32.o lapackc.f 
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o radau5_i32.o radau5.f
  
  In order to get create a shared library (from the object file above) use:

       gfortran -shared -o radau5_i32.dll
                 radau5_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
  
  """
function help_radau5_compile()
  return Docs.doc(help_radau5_compile)
end

function help_radau5_license()
  return Docs.doc(help_radau5_license)
end

@doc(@doc(hw_license),help_radau5_license)

"""  
  ## Compile RADAU 

  The julia ODEInterface tries to compile and link the solvers
  automatically at the build-time of this module. The following
  calls need only be done, if one uses a different compiler and/or if
  one wants to change/add some compiler options.

  The Fortran source code can be found at:
  
       http://www.unige.ch/~hairer/software.html 
  
  See `help_radau_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile RADAU with `Float64` reals and
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
                -o radau.o radau.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o radau.so 
                radau.o dc_lapack.o lapack.o lapackc.o
       gfortran -shared -fPIC -o radau.dylib
                radau.o dc_lapack.o lapack.o lapackc.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile RADAU with `Float64` reals and
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
                -o radau.o radau.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o radau.so 
                radau.o dc_lapack.o lapack.o lapackc.o
  
  ### Using `gfortran` and 32bit integers (Linux and Mac)
  
  Here is an example how to compile RADAU with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o lapackc_i32.o lapackc.f 
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o radau_i32.o radau.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):
  
       gfortran -shared -fPIC -o radau_i32.so 
                 radau_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
       gfortran -shared -fPIC -o radau_i32.dylib
                 radau_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
  
  ### Using `gfortran` and 32bit integers (Windows)
  
  Here is an example how to compile RADAU with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o dc_lapack_i32.o dc_lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o lapack_i32.o lapack.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o lapackc_i32.o lapackc.f 
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o radau_i32.o radau.f
  
  In order to get create a shared library (from the object file above) use:

       gfortran -shared -o radau_i32.dll
                 radau_i32.o dc_lapack_i32.o lapack_i32.o lapackc_i32.o
  
  """
function help_radau_compile()
  return Docs.doc(help_radau_compile)
end

function help_radau_license()
  return Docs.doc(help_radau_license)
end

@doc(@doc(hw_license),help_radau_license)

# Add informations about solvers in global solverInfo-array.
push!(solverInfo,
  SolverInfo("radau5",
    "Implicit Runge-Kutta method (Radau IIA) of order 5",
    tuple(:OPT_RTOL, :OPT_ATOL, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN, 
          :OPT_M1, :OPT_M2,
          :OPT_MASSMATRIX,
          :OPT_TRANSJTOH, :OPT_MAXSTEPS, :OPT_MAXNEWTONITER,
          :OPT_NEWTONSTARTZERO, 
          :OPT_DIMOFIND1VAR, :OPT_DIMOFIND2VAR, :OPT_DIMOFIND3VAR,
          :OPT_STEPSIZESTRATEGY, 
          :OPT_RHO, :OPT_JACRECOMPFACTOR, :OPT_NEWTONSTOPCRIT,
          :OPT_FREEZESSLEFT, :OPT_FREEZESSRIGHT,
          :OPT_SSMINSEL, :OPT_SSMAXSEL, :OPT_INITIALSS,
          :OPT_JACOBIMATRIX, :OPT_JACOBIBANDSTRUCT
          ),
    tuple(
      SolverVariant("radau5_i64",
        "Radau5 with 64bit integers",
        DL_RADAU5,
        tuple("radau5","contr5")),
      SolverVariant("radau5_i32",
        "Radau5 with 32bit integers",
        DL_RADAU5_I32,
        tuple("radau5","contr5")),
    ),
    help_radau5_compile,
    help_radau5_license,
  )
)

push!(solverInfo,
  SolverInfo("radau",
    "Implicit Runge-Kutta method (Radau IIA) of variable order ∈ (5,9,13)",
    tuple(:OPT_RTOL, :OPT_ATOL, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN, 
          :OPT_M1, :OPT_M2,
          :OPT_MASSMATRIX,
          :OPT_TRANSJTOH, :OPT_MAXSTEPS, :OPT_MAXNEWTONITER,
          :OPT_NEWTONSTARTZERO, 
          :OPT_DIMOFIND1VAR, :OPT_DIMOFIND2VAR, :OPT_DIMOFIND3VAR,
          :OPT_STEPSIZESTRATEGY, 
          :OPT_RHO, :OPT_JACRECOMPFACTOR, :OPT_NEWTONSTOPCRIT,
          :OPT_FREEZESSLEFT, :OPT_FREEZESSRIGHT,
          :OPT_SSMINSEL, :OPT_SSMAXSEL, :OPT_INITIALSS,
          :OPT_JACOBIMATRIX, :OPT_JACOBIBANDSTRUCT,
          :OPT_MAXSTAGES, :OPT_MINSTAGES, :OPT_INITSTAGES,
          :OPT_ORDERDECFACTOR, :OPT_ORDERINCFACTOR,
          :OPT_ORDERDECSTEPFAC1, :OPT_ORDERDECSTEPFAC2
          ),
    tuple(
      SolverVariant("radau_i64",
        "Radau with 64bit integers",
        DL_RADAU,
        tuple("radau","contra")),
      SolverVariant("radau5_i32",
        "Radau with 32bit integers",
        DL_RADAU_I32,
        tuple("radau","contra")),
    ),
    help_radau_compile,
    help_radau_license,
  )
)

# vim:syn=julia:cc=79:fdm=indent:
