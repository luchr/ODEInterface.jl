# Functions for ODE-Solver: ddeabm

"""Name for Loading ddeabm solver (64bit integers)."""
const DL_DDEABM = "ddeabm"

"""Name for Loading ddeabm solver (32bit integers)."""
const DL_DDEABM_I32 = "ddeabm_i32"

"""macro for import Ddeabm solver."""
macro import_ddeabm()
  :(
    using ODEInterface: ddeabm, ddeabm_i32
  )
end

"""macro for import Ddeabm dynamic lib names."""
macro import_DLddeabm()
  :(
    using ODEInterface: DL_DDEABM, DL_DDEABM_I32
  )
end

"""macro for importing Ddeabm help."""
macro import_ddeabm_help()
  :(
    using ODEInterface: help_ddeabm_compile, help_ddeabm_license
  )
end

"""
  Type encapsulating all required data for Ddeabm-Solver-Callbacks.

  We have the typical calling stack:

       ddeabm
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           loop with iterations (intermediate-mode and/or OPT_OUTPUTATTIMES)
           │ ccall( DDEABM_ ... ) # after 1st call: continuation call
           │    ┌───────────────────────────────────────────┐  ⎫
           │    │unsafe_SLATEC1RHSCallback                  │  ⎬ cb. rhs
           │    │    rhs                                    │  ⎪
           │    └───────────────────────────────────────────┘  ⎭
           └ call_julia_output_fcn( ... STEP ...)
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)

  see also help of `ODEInterface.SLATEC_continuation_call`.
  """
mutable struct DdeabmInternalCallInfos{FInt<:FortranInt,
        RHS_F<:Function, OUT_F<:Function} <: ODEinternalCallInfos
  logio        :: IO                    # where to log
  loglevel     :: UInt64                # log level
  # RHS:
  N            :: FInt                  # Dimension of Problem
  rhs          :: RHS_F                 # right-hand-side 
  rhs_mode     :: RHS_CALL_MODE         # how to call rhs
  rhs_lprefix  :: AbstractString        # saved log-prefix for rhs
  rhs_count    :: Int                   # count: number of calls to rhs
  # SOLOUT & output function
  output_mode  :: OUTPUTFCN_MODE        # what mode for output function
  output_fcn   :: OUT_F                 # the output function to call
  output_data  :: Dict                  # extra_data for output_fcn
end

"""
       mutable struct DdeabmArguments{FInt} <: AbstractArgumentsODESolver{FInt}

  Stores Arguments for Ddeabm solver.

  FInt is the Integer type used for the fortran compilation.
  """
mutable struct DdeabmArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
  FCN     :: Ptr{Void}         # rhs callback
  N       :: Vector{FInt}      # Dimension: NEQ
  t       :: Vector{Float64}   # start time (and current)
  x       :: Vector{Float64}   # initial value (and current state)
  tEnd    :: Vector{Float64}   # (current) end time
  INFO    :: Vector{FInt}      # INFO-Array
  RTOL    :: Vector{Float64}   # relative tolerance
  ATOL    :: Vector{Float64}   # absolute tolerance
  IDID    :: Vector{FInt}      # Status code
  RWORK   :: Vector{Float64}   # double working array
  LRW     :: Vector{FInt}      # length of WORK
  IWORK   :: Vector{FInt}      # integer working array
  LIW     :: Vector{FInt}      # length of IWORK
  RPAR    :: Vector{Float64}   # add. double-array
  IPAR    :: Ref{DdeabmInternalCallInfos} # misuse IPAR
    ## Allow uninitialized construction
  function DdeabmArguments{FInt}(dummy::FInt) where FInt
    return new{FInt}()
  end
end

"""
       function ddeabm(rhs::Function, t0::Real, T::Real,
                       x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)

  `retcode` can have the following values:

        1: computation successful
        2: computation. successful, but interrupted by output function
       <0: error

  main call for using Fortran-ddeabm solver. In `opt` the following
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
      ║                 │   call OPT_OUTPUTFCN either              │         ║
      ║                 │   (a) either for all intermediate steps  │         ║
      ║                 │       choosen by the solver or           │         ║
      ║                 │   (b) at the times given in the option   │         ║
      ║                 │       OPT_OUTPUTATTIMES                  │         ║
      ║                 │ OUTPUTFCN_DENSE                          │         ║
      ║                 │   is *not* supported!                    │         ║
      ║                 │   but see OUTPUTATTIMES for an           │         ║
      ║                 │   alternative approach                   │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ OUTPUTATTIMES   │ If OPT_OUTPUTMODE is OUTPUTFCN_WODENSE   │ nothing ║
      ║                 │ then one can specify with this vector    │         ║
      ║                 │ the time points where the OPT_OUTPUTFCN  │         ║
      ║                 │ should be called.                        │         ║
      ║                 │ All values of OPT_OUTPUTATTIMES *must*   │         ║
      ║                 │ be sorted (ascending, if T>t0, and       │         ║
      ║                 │ descending, if T<t0) and they must be    │         ║
      ║                 │ between t0 and T.                        │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ TSTOP           │ tell the solver, that it is not          │     NaN ║
      ║                 │ permissable to integrate past the point  │         ║
      ║                 │ TSTOP. If TSTOP is NaN then the solver   │         ║
      ║                 │ may integrate past T and interpolate the │         ║
      ║                 │ result at T. Sometimes there are         │         ║
      ║                 │ right-hand sides, where this is not      │         ║
      ║                 │ possible.                                │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTEPS        │ maximal number of allowed steps          │  100000 ║
      ║                 │ (allowed intermediate steps)             │         ║
      ║                 │ between t0, T and the values given       │         ║
      ║                 │ in OPT_OUTPUTATTIMES.                    │         ║
      ║                 │ The value will be rounded up to a        │         ║
      ║                 │ multiple of 500.                         │         ║
      ║                 │ OPT_MAXSTEPS > 0                         │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝ 

  """
function ddeabm(rhs::Function, t0::Real, T::Real,
                  x0::Vector, opt::AbstractOptionsODE)
  return ddeabm_impl(rhs, t0, T, x0, opt, DdeabmArguments{Int64}(Int64(0)))
end

"""
  ddeabm with 32bit integers, see ddeabm.
  """
function ddeabm_i32(rhs::Function, t0::Real, T::Real,
                  x0::Vector, opt::AbstractOptionsODE)
  return ddeabm_impl(rhs, t0, T, x0, opt, DdeabmArguments{Int32}(Int32(0)))
end

"""
  MAXNUM value in ddeabm.
  """
const ddeabm_maxnum = 500

"""
       function ddeabm_impl{FInt<:FortranInt}(rhs::Function, 
                t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
                args::DdeabmArguments{FInt})
  
  implementation of ddeabm-call for FInt.

  This solver has the conecpt of continuation calls (CC), see
  help for `ODEInterface.SLATEC_continuation_call`.

  This CC are used in this method to get the solution at user-given
  t-values (see `OPT_OUTPUTATTIMES`) in order to support some alternative
  method for dense output.

  The solver also has an intermediate-output mode. If activated then
  the numerical integration is stoped at every (solver-chosen) intermediate
  step and can be continued with a CC.

  Different cases for this method:

      ╔═══════════════════╤═══════════════════╤════════════════════╤═════════╗
      ║ OPT_OUTPUTMODE    │ OPT_OUTPUTATTIMES │ intermediate-mode? │  case   ║
      ╠═══════════════════╪═══════════════════╪════════════════════╪═════════╣
      ║ OUTPUTFCN_NEVER   │ <ignored>         │ no   INFO(3)=0     │  C1     ║
      ║ OUTPUTFCN_WODENSE │ nothing or empty  │ yes  INFO(3)=1     │  C2     ║
      ║ OUTPUTFCN_WODENSE │ given vector      │ no   INFO(3)=0     │  C3     ║
      ╠═══════════════════╪═══════════════════╧════════════════════╧═════════╣
      ║ OUTPUTFCN_DENSE   │        N O T    S U P P O R T E D !              ║
      ╚═══════════════════╧══════════════════════════════════════════════════╝ 

  """
function ddeabm_impl{FInt<:FortranInt}(rhs::Function, 
        t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
        args::DdeabmArguments{FInt})

  (lio,l,l_g,l_solver,lprefix) = solver_start("ddeabm",rhs,t0,T,x0,opt)

  method_ddeabm, = getAllMethodPtrs(
    (FInt == Int64) ? DL_DDEABM : DL_DDEABM_I32)

  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)

  check_slatec_output_mode(output_mode)

  # RWORK
  args.LRW = [ FInt(130+21*d) ]
  args.RWORK = zeros(Float64, args.LRW[1])

  # IWORK
  args.LIW = [ FInt(51) ]
  args.IWORK = zeros(FInt, args.LIW[1])

  # INFO
  args.INFO = zeros(FInt, 15)

  args.INFO[1] = 0
  args.INFO[2] = scalarFlag ? 0 : 1

  tstop = NaN
  maxsteps = 0
  try
    maxsteps = convert(Int, getOption(opt, OPT_MAXSTEPS, 100000))
    @assert 0<maxsteps
  catch e
    throw(ArgumentErrorODE("OPT_MAXSTEPS is not valid", :opt, e))
  end

  try
    tstop = convert(Float64, getOption(opt, OPT_TSTOP, NaN))
  catch e
    throw(ArgumentErrorODE("Cannot convert OPT_TSTOP to Float64", :opt, e))
  end
  args.INFO[4] = isnan(tstop) ? 0 : 1
  if !isnan(tstop)
    args.RWORK[1] = tstop
  end

  t_values = extractSlatecOutputAtTimes(t0, T, opt)
  output_attimes_given = length(t_values)>0
  push!(t_values, args.tEnd[1])

  args.RPAR = zeros(Float64, 0)
  args.IDID = zeros(FInt, 1)
  rhs_lprefix = "unsafe_SLATEC1RHSCallback: "

  case = output_mode == OUTPUTFCN_NEVER ? 1 : (output_attimes_given ? 3 : 2)

  args.INFO[3] = (case == 2) ? 1 : 0
  retcode = -100

  cbi = DdeabmInternalCallInfos(lio, l, d, rhs, rhs_mode, rhs_lprefix, 
        0, output_mode, output_fcn, Dict())

  args.FCN = unsafe_SLATEC1RHSCallback_c(cbi, FInt(0))
  args.IPAR = cbi

  if output_mode ≠ OUTPUTFCN_NEVER
    call_julia_output_fcn(cbi, OUTPUTFCN_CALL_INIT,
      args.t[1], args.tEnd[1], args.x, eval_sol_fcn_init) # ignore result
    call_julia_output_fcn(cbi, OUTPUTFCN_CALL_STEP,
      args.t[1], args.t[1], args.x, eval_sol_fcn_noeval)  # 1 step
  end

  if l_solver
    println(lio,lprefix,"call Fortran-ddeabm $method_ddeabm with")
    dump(lio,args);
  end

  maxsteps_seen = 0
  while (true) 
    told = args.t[1]
    args.tEnd[1] = t_values[1]
    ccall( method_ddeabm, Void,
      (Ptr{Void}, Ptr{FInt},                        # Rightsidefunc, N=NEQ=d
       Ptr{Float64}, Ptr{Float64}, Ptr{Float64},    # t, x, tEnd
       Ptr{FInt}, Ptr{Float64}, Ptr{Float64},       # INFO, RTOL, ATOL
       Ptr{FInt},                                   # IDID
       Ptr{Float64}, Ptr{FInt},                     # RWORK, LRW
       Ptr{FInt}, Ptr{FInt},                        # IWORK, LIW
       Ptr{Float64}, Ref{DdeabmInternalCallInfos},  # RPAR, IPAR=CBI
      ),
       args.FCN, args.N,
       args.t, args.x, args.tEnd,
       args.INFO, args.RTOL, args.ATOL,
       args.IDID,
       args.RWORK, args.LRW,
       args.IWORK, args.LIW,
       args.RPAR, args.IPAR)
    if l_solver
      println(lio,lprefix,"call Fortran-ddeabm $method_ddeabm returned")
      dump(lio,args);
    end
    if args.IDID[1] < -1     # -1 handled below
      retcode = Int(args.IDID[1])
      break
    end
    if args.IDID[1] ≥ 0 && output_mode ≠ OUTPUTFCN_NEVER
      out_result = call_julia_output_fcn(cbi, OUTPUTFCN_CALL_STEP,
        told, args.t[1], args.x, eval_sol_fcn_noeval)
      if out_result == OUTPUTFCN_RET_CONTINUE_XCHANGED
        throw(FeatureNotSupported(string("Sorry; ddeabm does not support ",
              "to change the solution inside the output function.")))
      elseif out_result == OUTPUTFCN_RET_STOP
        retcode = 2
        break
      end
    end
    if args.IDID[1] == -1
      # => MAXNUM steps done
      maxsteps_seen += ddeabm_maxnum
      if maxsteps_seen ≥ maxsteps
        retcode = -1
        break
      end
    elseif args.IDID[1] == 1
      # => intermediate-output mode => case C2
      # args.tEnd[1] yet not reached => continue
      args.INFO[1] = 1
    elseif args.IDID[1] ∈ (2,3,)
      # => args.tEnd[1] reached
      shift!(t_values)     # next t-value to compute
      maxsteps_seen = 0
      if length(t_values) == 0
        retcode = 1   # reached T
        break
      end
    else
      throw(InternalErrorODE(
        string("unknown IDID=",args.IDID[1]," result of ddeabm")))
    end
    args.INFO[1] = 1  # continuation call
  end

  output_mode ≠ OUTPUTFCN_NEVER &&
    call_julia_output_fcn(cbi,OUTPUTFCN_CALL_DONE,
      args.t[1], args.tEnd[1], args.x, eval_sol_fcn_done) # ignore result

  l_g && println(lio,lprefix,string("done IDID=",args.IDID[1]))
  stats = Dict{AbstractString,Any}(
    "no_rhs_calls"       => cbi.rhs_count,
    "ddeabm_idid"        => args.IDID[1],
        )

  return (args.t[1], args.x, retcode, stats)
end

"""
  ## Compile DDEABM

  The julia ODEInterface tries to compile and link the solvers
  automatically at the build-time of this module. The following
  calls need only be done, if one uses a different compiler and/or if
  one wants to change/add some compiler options.

  The Fortran source code can be found at:

       http://www.netlib.org/slatec/src/

  See `help_ddeabm_license` for the licsense information.

  ### Using `gfortran` and 64bit integers (Linux, Mac and Windows)
  
  Here is an example how to compile DDEABM with `Float64` reals and
  `Int64` integers with `gfortran`:

       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 -o slatec.o slatec.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 -o ddeabm.o ddeabm.f

  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac, 3rd for Window):

       gfortran -shared -fPIC -o ddeabm.so ddeabm.o slatec.o
       gfortran -shared -fPIC -o ddeabm.dylib ddeabm.o slatec.o
       gfortran -shared       -o ddeabm.dll ddeabm.o slatec.o
  """
function help_ddeabm_compile()
  return Docs.doc(help_ddeabm_compile)
end

function help_ddeabm_license()
  return Docs.doc(help_ddeabm_license)
end

@doc(@doc(SLATEC_license), help_ddeabm_license)


# Add informations about solver in global solverInfo-array.
push!(solverInfo,
  SolverInfo("ddeabm",
    "Adams-Bashforth-Moulton Predictor-Corrector method (orders: 1-12)",
    tuple(:OPT_RTOL, :OPT_ATOL, :OPT_RHS_CALLMODE, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN,
          :OPT_OUTPUTATTIMES,
          :OPT_TSTOP, :OPT_MAXSTEPS,
    ),
    tuple(
      SolverVariant("ddeabm_i64",
        "Ddeabm with 64bit integers",
        DL_DDEABM,
        tuple("ddeabm",)),
      SolverVariant("ddeabm_i32",
        "Ddeabm with 32bit integers",
        DL_DDEABM_I32,
        tuple("ddeabm",)),
    ),
    help_ddeabm_compile,
    help_ddeabm_license,
  )
)


# vim:syn=julia:cc=79:fdm=indent:
