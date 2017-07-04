# Functions for ODE-Solver: ddebdf

"""Name for Loading ddebdf solver (64bit integers)."""
const DL_DDEBDF = "ddebdf"

"""Name for Loading ddebdf solver (32bit integers)."""
const DL_DDEBDF_I32 = "ddebdf_i32"

"""macro for import Ddebdf solver."""
macro import_ddebdf()
  :(
    using ODEInterface: ddebdf, ddebdf_i32
  )
end

"""macro for import Ddebdf dynamic lib names."""
macro import_DLddebdf()
  :(
    using ODEInterface: DL_DDEBDF, DL_DDEBDF_I32
  )
end

"""macro for importing Ddebdf help."""
macro import_ddebdf_help()
  :(
    using ODEInterface: help_ddebdf_compile, help_ddebdf_license
  )
end

"""
  Type encapsulating all required data for Ddebdf-Solver-Callbacks.

  We have the typical calling stack:

       ddebdf
           call_julia_output_fcn(  ... INIT ... )
               output_fcn ( ... INIT ...)
           loop with iterations (intermediate-mode and/or OPT_OUTPUTATTIMES)
           │ ccall( DDEBDF_ ... ) # after 1st call: continuation call
           │    ┌───────────────────────────────────────────┐  ⎫
           │    │unsafe_SLATEC1RHSCallback                  │  ⎬ cb. rhs
           │    │    rhs                                    │  ⎪
           │    └───────────────────────────────────────────┘  ⎭
           │    ┌───────────────────────────────────────────┐  ⎫
           │    │unsafe_SLATEC1JacCallback                  │  ⎬ cb. jacobian
           │    │    call_julia_jac_fcn(             )      │  ⎪
           │    └───────────────────────────────────────────┘  ⎭
           └ call_julia_output_fcn( ... STEP ...)
           call_julia_output_fcn(  ... DONE ... )
               output_fcn ( ... DONE ...)

  see also help of `ODEInterface.SLATEC_continuation_call`.
  """
type DdebdfInternalCallInfos{FInt<:FortranInt,
        RHS_F<:Function, OUT_F<:Function,
        JAC_F<:Function, VIEW_F<:Function} <: ODEinternalCallInfos
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
  # Jacobimatrix
  jacobimatrix :: JAC_F                 # function for jacobi matrix
  jacobibandstruct                      # (l,u) bandstructure or nothing
  jac_lprefix  :: AbstractString        # saved log-prefix for jac
  #
  submatrix    :: VIEW_F                # method for generating view/sub
end

"""
       type DdebdfArguments{FInt} <: AbstractArgumentsODESolver{FInt}

  Stores Arguments for Ddebdf solver.

  FInt is the Integer type used for the fortran compilation.
  """
type DdebdfArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
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
  IPAR    :: Ref{DdebdfInternalCallInfos} # misuse IPAR
  JAC     :: Ptr{Void}         # Jacobian callback
    ## Allow uninitialized construction
  @WHEREFUNC(FInt,
  function DdebdfArguments{FInt}(dummy::FInt)
    return new{FInt}()
  end
  )
end

"""
       function ddebdf(rhs::Function, t0::Real, T::Real,
                       x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)

  `retcode` can have the following values:

        1: computation successful
        2: computation. successful, but interrupted by output function
       <0: error

  main call for using Fortran-ddebdf solver. In `opt` the following
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
      ║ JACOBIMATRIX    │ A function providing the Jacobian for    │ nothing ║
      ║                 │ ∂f/∂x or nothing. For nothing the solver │         ║
      ║                 │ uses finite differences to calculate the │         ║
      ║                 │ Jacobian.                                │         ║
      ║                 │ The function has to be of the form:      │         ║
      ║                 │   function (t,x,J) -> nothing            │         ║
      ║                 │ Depending on OPT_JACOBIBANDSTRUCT the    │         ║
      ║                 │ argument J will then by a full or a      │         ║
      ║                 │ banded matrix, where the user-given      │         ║
      ║                 │ function has to fill in the entries.     │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ JACOBIBANDSTRUCT│ A tuple (l,u) describing the banded      │ nothing ║
      ║                 │ structure of the Jacobian or nothing if  │         ║
      ║                 │ the Jacobian is full.                    │         ║
      ║                 │ Even if the option JACOBIMATRIX is empty,│         ║
      ║                 │ the solver will perform much better if   │         ║
      ║                 │ the Jacobian matrix is banded and the    │         ║
      ║                 │ code is told this.                       │         ║
      ║                 │ see also help of BandedMatrix            │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝ 

  """
function ddebdf(rhs::Function, t0::Real, T::Real,
                  x0::Vector, opt::AbstractOptionsODE)
  return ddebdf_impl(rhs, t0, T, x0, opt, DdebdfArguments{Int64}(Int64(0)))
end

"""
  ddebdf with 32bit integers, see ddebdf.
  """
function ddebdf_i32(rhs::Function, t0::Real, T::Real,
                  x0::Vector, opt::AbstractOptionsODE)
  return ddebdf_impl(rhs, t0, T, x0, opt, DdebdfArguments{Int32}(Int32(0)))
end


"""
       function ddebdf_impl{FInt<:FortranInt}(rhs::Function, 
                t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
                args::DdebdfArguments{FInt})
  
  implementation of ddebdf-call for FInt.

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
function ddebdf_impl{FInt<:FortranInt}(rhs::Function, 
        t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE, 
        args::DdebdfArguments{FInt})

  (lio,l,l_g,l_solver,lprefix) = solver_start("ddeabm",rhs,t0,T,x0,opt)

  method_ddebdf, = getAllMethodPtrs(
    (FInt == Int64) ? DL_DDEBDF : DL_DDEBDF_I32)

  (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn) =
    solver_extract_commonOpt(t0,T,x0,opt,args)

  check_slatec_output_mode(output_mode)

  (jacobimatrix, jacobibandstruct) = extractSlatecJacobiOpt(d, opt)

  # INFO
  args.INFO = zeros(FInt, 15)

  # RWORK
  args.LRW = [ FInt( 
    jacobibandstruct == nothing ? 
      250+10*d+d*d : 250+10*d+(2*jacobibandstruct[1]+jacobibandstruct[2]+1)*d
    ) ]
  args.RWORK = zeros(Float64, args.LRW[1])

  # IWORK
  args.LIW = [ FInt(56+d) ]
  args.IWORK = zeros(FInt, args.LIW[1])

  args.INFO[1] = 0
  args.INFO[2] = scalarFlag ? 0 : 1

  tstop = NaN
  try
    tstop = convert(Float64, getOption(opt, OPT_TSTOP, NaN))
  catch e
    throw(ArgumentErrorODE("Cannot convert OPT_TSTOP to Float64", :opt, e))
  end
  args.INFO[4] = isnan(tstop) ? 0 : 1
  if !isnan(tstop)
    args.RWORK[1] = tstop
  end

  args.INFO[5] = jacobimatrix == nothing ? 0 : 1
  if jacobibandstruct ≠ nothing
    args.INFO[6] = 1
    args.IWORK[1] = jacobibandstruct[1]
    args.IWORK[2] = jacobibandstruct[2]
  else
    args.INFO[6] = 0
  end

  t_values = extractSlatecOutputAtTimes(t0, T, opt)
  output_attimes_given = length(t_values)>0
  push!(t_values, args.tEnd[1])

  args.RPAR = zeros(Float64, 0)
  args.IDID = zeros(FInt, 1)
  rhs_lprefix = "unsafe_SLATEC1RHSCallback: "
  jac_lprefix = "unsafe_SLATEC1JacCallback: "

  case = output_mode == OUTPUTFCN_NEVER ? 1 : (output_attimes_given ? 3 : 2)

  args.INFO[3] = (case == 2) ? 1 : 0
  retcode = -100

  cbi = DdebdfInternalCallInfos(lio, l, d, rhs, rhs_mode, rhs_lprefix,
        0, output_mode, output_fcn, Dict(), 
        jacobimatrix == nothing ? dummy_func : jacobimatrix,
        jacobibandstruct, jac_lprefix, get_view_function())

  args.IPAR = cbi
  args.FCN = unsafe_SLATEC1RHSCallback_c(cbi, FInt(0))
  args.JAC = unsafe_SLATEC1JacCallback_c(cbi, FInt(0))

  if output_mode ≠ OUTPUTFCN_NEVER
    call_julia_output_fcn(cbi, OUTPUTFCN_CALL_INIT,
      args.t[1], args.tEnd[1], args.x, eval_sol_fcn_init) # ignore result
    call_julia_output_fcn(cbi, OUTPUTFCN_CALL_STEP,
      args.t[1], args.t[1], args.x, eval_sol_fcn_noeval)  # 1 step
  end

  if l_solver
    println(lio,lprefix,"call Fortran-ddeabm $method_ddebdf with")
    dump(lio,args);
  end

  while (true)
    told = args.t[1]
    args.tEnd[1] = t_values[1]
    ccall( method_ddebdf, Void,
      (Ptr{Void}, Ptr{FInt},                        # Rightsidefunc, N=NEQ=d
       Ptr{Float64}, Ptr{Float64}, Ptr{Float64},    # t, x, tEnd
       Ptr{FInt}, Ptr{Float64}, Ptr{Float64},       # INFO, RTOL, ATOL,
       Ptr{FInt},                                   # IDID
       Ptr{Float64}, Ptr{FInt},                     # RWORK, LRW
       Ptr{FInt}, Ptr{FInt},                        # IWORK, LIW
       Ptr{Float64}, Ref{DdebdfInternalCallInfos},  # RPAR, IPAR=CBI
       Ptr{Void},                                   # DJAC
      ),
       args.FCN, args.N,
       args.t, args.x, args.tEnd,
       args.INFO, args.RTOL, args.ATOL,
       args.IDID,
       args.RWORK, args.LRW,
       args.IWORK, args.LIW,
       args.RPAR, args.IPAR,
       args.JAC)
    if l_solver
      println(lio,lprefix,"call Fortran-ddeabm $method_ddebdf returned")
      dump(lio,args);
    end
    if args.IDID[1] < 0
      retcode = Int(args.IDID[1])
      break
    end
    if output_mode ≠ OUTPUTFCN_NEVER
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
    if args.IDID[1] == 1
      # => intermediate-output mode => case C2
      # args.tEnd[1] yet not reached => continue
      args.INFO[1] = 1
    elseif args.IDID[1] ∈ (2,3,)
      # => args.tEnd[1] reached
      shift!(t_values)     # next t-value to compute
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
    "ddebdf_idid"        => args.IDID[1],
        )

  return (args.t[1], args.x, retcode, stats)
end

"""
  ## Compile DDEBDF

  The julia ODEInterface tries to compile and link the solvers
  automatically at the build-time of this module. The following
  calls need only be done, if one uses a different compiler and/or if
  one wants to change/add some compiler options.

  The Fortran source code can be found at:

       http://www.netlib.org/slatec/src/

  See `help_ddebdf_license` for the licsense information.

  ### Using `gfortran` and 64bit integers (Linux, Mac and Windows)
  
  Here is an example how to compile DDEBDF with `Float64` reals and
  `Int64` integers with `gfortran`:

       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 -o d1mach.o d1mach.f

  Do the same for the other components: d1mach, daxpy, dcfod, dgesl, dhstrt,
  dhvnrm, dgbfa, dgefa, dintyd, dlsod, dpjac, dslvs, dscal, dstod,  dvnrms,
  fdump, i1mach, idamax, j4save, xercnt, xerhlt, xermsg, xerprn, xersve, xgetua

  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac, 3rd for Window):

       gfortran -shared -fPIC -o ddebdf.so ddebdf.o d1mach.o daxpy.o dcfod.o
                dgefa.o dgesl.o dhstrt.o dhvnrm.o dintyd.o dlsod.o dpjac.o 
                dslvs.o dscal.o dstod.o dvnrms.o fdump.o i1mach.o idamax.o
                j4save.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o xgetua.o
       gfortran -shared -fPIC -o ddebdf.dylib ddebdf.o d1mach.o daxpy.o dcfod.o
                dgefa.o dgesl.o dhstrt.o dhvnrm.o dintyd.o dlsod.o dpjac.o 
                dslvs.o dscal.o dstod.o dvnrms.o fdump.o i1mach.o idamax.o
                j4save.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o xgetua.o
       gfortran -shared       -o ddebdf.dll ddebdf.o d1mach.o daxpy.o dcfod.o
                dgefa.o dgesl.o dhstrt.o dhvnrm.o dintyd.o dlsod.o dpjac.o 
                dslvs.o dscal.o dstod.o dvnrms.o fdump.o i1mach.o idamax.o
                j4save.o xercnt.o xerhlt.o xermsg.o xerprn.o xersve.o xgetua.o
  """
function help_ddebdf_compile()
  return Docs.doc(help_ddebdf_compile)
end

function help_ddebdf_license()
  return Docs.doc(help_ddebdf_license)
end

@doc(@doc(SLATEC_license), help_ddebdf_license)

# Add informations about solver in global solverInfo-array.
push!(solverInfo,
  SolverInfo("ddebdf",
    "Backward Differentiation Formula (orders: 1-5)",
    tuple(:OPT_RTOL, :OPT_ATOL, :OPT_RHS_CALLMODE, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN,
          :OPT_OUTPUTATTIMES,
          :OPT_TSTOP,
          :OPT_JACOBIMATRIX, :OPT_JACOBIBANDSTRUCT,
    ),
    tuple(
      SolverVariant("ddebdf_i64",
        "Ddebdf with 64bit integers",
        DL_DDEBDF,
        tuple("ddebdf",)),
      SolverVariant("ddebdf_i32",
        "Ddebdf with 32bit integers",
        DL_DDEBDF_I32,
        tuple("ddebdf",)),
    ),
    help_ddebdf_compile,
    help_ddebdf_license,
  )
)


# vim:syn=julia:cc=79:fdm=indent:
