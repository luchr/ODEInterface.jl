# Functions for ODE-Solver: Dop853

"""Name for Loading dop853 solver (64bit integers)."""
const DL_DOP853               = "dop853"

"""Name for Loading dop853 solver (32bit integers)."""
const DL_DOP853_I32           = "dop853_i32"

"""macro for import Dop853 solver."""
macro import_dop853()
  :(
    using ODEInterface: dop853, dop853_i32
  )
end

"""macro for import Dop853 dynamic lib names."""
macro import_DLdop853()
  :(
    using ODEInterface: DL_DOP853, DL_DOP853_I32
  )
end

"""macro for importing Dopri5 help."""
macro import_dop853_help()
  :(
    using ODEInterface: help_dop853_compile, help_dop853_license
  )
end

"""
      function dop853(rhs::Function, t0::Real, T::Real,
                      x0::Vector, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)

  `retcode` can have the following values:

        1: computation successful
        2: computation. successful, but interrupted by output function
       -1: input is not consistent
       -2: larger OPT_MAXSTEPS is needed
       -3: step size becomes too small
       -4: problem is probably stiff (interrupted)
  
  main call for using Fortran-dopri5 solver. In `opt` the following
  options are used:

      ╔═════════════════╤══════════════════════════════════════════╤═════════╗
      ║  Option OPT_…   │ Description                              │ Default ║
      ╠═════════════════╪══════════════════════════════════════════╪═════════╣
      ║ RTOL     &      │ relative and absolute error tolerances   │    1e-3 ║
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
      ║ MAXSTEPS        │ maximal number of allowed steps          │  100000 ║
      ║                 │ OPT_MAXSTEPS > 0                         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ STEST           │ stiffness test                           │    1000 ║
      ║                 │ done after every step number k*OPT_STEST │         ║
      ║                 │ OPT_STEST < 0 for turning test off       │         ║
      ║                 │ OPT_STEST ≠ 0                            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ EPS             │ the rounding unit                        │ 2.3e-16 ║
      ║                 │ 1e-35 < OPT_EPS < 1.0                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ RHO             │ safety factor in step size predcition    │     0.9 ║
      ║                 │ 1e-4  < OPT_RHO < 1.0                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SSMINSEL   &    │ parameters for step size selection       │   0.333 ║
      ║ SSMAXSEL        │ The new step size is chosen subject to   │     6.0 ║
      ║                 │ the restriction                          │         ║
      ║                 │ OPT_SSMINSEL ≤ hnew/hold ≤ OPT_SSMAXSEL  │         ║
      ║                 │ OPT_SSMINSEL, OPT_SSMAXSEL > 0           │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SSBETA          │ β for stabilized step size control       │     0.0 ║
      ║                 │ OPT_SSBETA ≤ 0.2                         │         ║
      ║                 │ if OPT_SSBETA < 0 then OPT_SSBETA = 0    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSS           │ maximal step size                        │  T - t0 ║
      ║                 │ OPT_MAXSS ≠ 0                            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ INITIALSS       │ initial step size                        │     0.0 ║
      ║                 │ if OPT_INITIALSS == 0 then a initial     │         ║
      ║                 │ guess is computed                        │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝ 
  
  """
function dop853(rhs::Function, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return dop853_impl(rhs,t0,T,x0,opt,DopriArguments{Int64}())
end

"""
  dop853 with 32bit integers, see dop853
  """
function dop853_i32(rhs::Function, t0::Real, T::Real,
                x0::Vector, opt::AbstractOptionsODE)
  return dop853_impl(rhs,t0,T,x0,opt,DopriArguments{Int32}())
end

"""
       function dop853_impl{FInt}(rhs::Function, t0::Real, T::Real, x0::Vector,
                        opt::AbstractOptionsODE,arg::DopriArguments{FInt})
  
  implementation of dop853 for FInt ∈ (Int32,Int4)
  """
function dop853_impl{FInt}(rhs::Function, t0::Real, T::Real, x0::Vector, 
                 opt::AbstractOptionsODE,args::DopriArguments{FInt})

  (lio,l,l_g,l_solver,lprefix,cid,cid_str) = 
    solver_start("dop853",rhs,t0,T,x0,opt)

  FInt ∉ (Int32,Int64) &&
    throw(ArgumentErrorODE("only FInt ∈ (Int32,Int4) allowed"))
  
  (method_dop853, method_contd8) = getAllMethodPtrs(
     (FInt == Int64)? DL_DOP853 : DL_DOP853_I32 )

  (d,nrdense,rhs_mode,output_mode,output_fcn) = 
    dopri_extract_commonOpt(t0,T,x0,opt,args)

  # WORK memory
  args.LWORK = [ 11*d+8*nrdense+21 ]
  args.WORK = zeros(Float64,args.LWORK[1])

  # IWORK memory
  args.LIWORK = [ nrdense+21 ]
  args.IWORK = zeros(FInt,args.LIWORK[1])

  # fill IWORK
  OPT = nothing
  try
    OPT=OPT_MAXSTEPS; args.IWORK[1] = convert(FInt,getOption(opt,OPT,100000))
    @assert 0 < args.IWORK[1]
    args.IWORK[2]=1        # choice for coefficients
    args.IWORK[3]=-1       # don't print anything
    OPT=OPT_STEST; args.IWORK[4]=convert(FInt,getOption(opt,OPT_STEST,1000))
    @assert 0 ≠ args.IWORK[4]
    args.IWORK[5]=nrdense
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  # fill WORK
  try
    OPT=OPT_EPS; args.WORK[1]=convert(Float64,getOption(opt,OPT,2.3e-16))
    @assert 1e-35 < args.WORK[1] < 1.0
    OPT=OPT_RHO; args.WORK[2]=convert(Float64,getOption(opt,OPT,0.9))
    @assert 1e-4  < args.WORK[2] < 1.0
    OPT=OPT_SSMINSEL; args.WORK[3]=convert(Float64,getOption(opt,OPT,0.333))
    @assert 0.0 < args.WORK[3]
    OPT=OPT_SSMAXSEL; args.WORK[4]=convert(Float64,getOption(opt,OPT,6.0))
    @assert 0.0 < args.WORK[4]
    OPT=OPT_SSBETA; args.WORK[5]=convert(Float64,getOption(opt,OPT,0.0))
    @assert args.WORK[5] ≤ 0.2
    OPT=OPT_MAXSS; args.WORK[6]=convert(Float64,getOption(opt,OPT,T-t0))
    @assert 0 ≠ args.WORK[6]
    OPT=OPT_INITIALSS; args.WORK[7]=convert(Float64,getOption(opt,OPT,0.0))
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  args.RPAR=zeros(Float64,0)
  args.IDID=zeros(FInt,1)
  args.FCN    = (FInt == Int64)? unsafe_HW1RHSCallback_c : 
                                 unsafe_HW1RHSCallbacki32_c
  rhs_lprefix = string(cid_str,"unsafe_HW1RHSCallback: ")
  args.SOLOUT = (FInt == Int64)? unsafe_dopriSoloutCallback_c:
                                 unsafe_dopriSoloutCallbacki32_c
  out_lprefix = string(cid_str,"unsafe_dopriSoloutCallback: ")
  eval_lprefix = string(cid_str,"eval_sol_fcn_closure: ")

  try
    eval_sol_fcn = 
      (output_mode == OUTPUTFCN_DENSE)?
        create_dopri_eval_sol_fcn_closure(cid[1],d,method_contd8):
        eval_sol_fcn_noeval

    GlobalCallInfoDict[cid[1]] = 
      DopriInternalCallInfos{FInt}(cid,lio,l,rhs,rhs_mode,rhs_lprefix, 
        output_mode,output_fcn,
        Dict(), out_lprefix, eval_sol_fcn,eval_lprefix,
        NaN,NaN,Vector{Float64}(),
        Vector{FInt}(1),Vector{Float64}(1),C_NULL,C_NULL,C_NULL)

    args.IPAR = Vector{FInt}(2)   # enough for cid[1] even in 32bit case
    packUInt64ToVector!(args.IPAR,cid[1])

    output_mode ≠ OUTPUTFCN_NEVER &&
      call_julia_output_fcn(cid[1],OUTPUTFCN_CALL_INIT,
        args.t[1],args.tEnd[1],args.x,eval_sol_fcn_init) # ignore result

    if l_solver
      println(lio,lprefix,"call Fortran-dop853 $method_dop853 with")
      dump(lio,args);
    end

    ccall( method_dop853, Void,
      (Ptr{FInt},  Ptr{Void},                    # N=d, Rightsidefunc
       Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, # t, x, tEnd
       Ptr{Float64}, Ptr{Float64}, Ptr{FInt},    # RTOL, ATOL, ITOL
       Ptr{Void}, Ptr{FInt},                     # Soloutfunc, IOUT
       Ptr{Float64}, Ptr{FInt},                  # WORK, LWORK
       Ptr{FInt}, Ptr{FInt},                     # IWORK, LIWORK
       Ptr{Float64}, Ptr{FInt}, Ptr{FInt},       # RPAR, IPAR, IDID
      ),
      args.N, args.FCN, 
      args.t, args.x, args.tEnd,
      args.RTOL, args.ATOL, args.ITOL, 
      args.SOLOUT, args.IOUT,
      args.WORK, args.LWORK,
      args.IWORK, args.LIWORK,
      args.RPAR, args.IPAR, args.IDID,
    )

    if l_solver
      println(lio,lprefix,"Fortran-dop853 $method_dop853 returned")
      dump(lio,args);
    end

    output_mode ≠ OUTPUTFCN_NEVER &&
      call_julia_output_fcn(cid[1],OUTPUTFCN_CALL_DONE,
        args.t[1],args.tEnd[1],args.x,eval_sol_fcn_done) # ignore result

  finally
    delete!(GlobalCallInfoDict,cid[1])
  end
  l_g && println(lio,lprefix,string("done IDID=",args.IDID[1]))
  stats = Dict{ASCIIString,Any}(
    "no_rhs_calls"       => args.IWORK[17],
    "no_steps"           => args.IWORK[18],
    "no_steps_accepted"  => args.IWORK[19],
    "no_steps_rejected"  => args.IWORK[20]
          )
  return ( args.t[1], args.x, args.IDID[1], stats)

end

"""  
  ## Compile DOP853 

  The Fortran source code can be found at:
  
       http://www.unige.ch/~hairer/software.html 
  
  See `help_dop853_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile DOP853 with `Float64` reals and
  `Int64` integers with `gfortran`:

       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o dop853.o dop853.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):
  
       gfortran -shared -fPIC -o dop853.so dop853.o
       gfortran -shared -fPIC -o dop853.dylib dop853.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile DOP853 with `Float64` reals and
  `Int64` integers with `gfortran`:

       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o dop853.o dop853.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o dop853.dll dop853.o
  
  ### Using `gfortran` and 32bit integers (Linux and Mac)
  
  Here is an example how to compile DOP853 with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fPIC  
                -fdefault-real-8 -fdefault-double-8 
                -o dop853_i32.o   dop853.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o dop853_i32.so dop853_i32.o
       gfortran -shared -fPIC -o dop853_i32.dylib dop853_i32.o
  
  ### Using `gfortran` and 32bit integers (Windows)
  
  Here is an example how to compile DOP853 with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c
                -fdefault-real-8 -fdefault-double-8 
                -o dop853_i32.o   dop853.f
  
  In order to get create a shared library (from the object file above) use:

       gfortran -shared -o dop853_i32.dll dop853_i32.o
  
  """
function help_dop853_compile()
  return Docs.doc(help_dop853_compile)
end

function help_dop853_license()
  return Docs.doc(help_dop853_license)
end

@doc(@doc(hw_license),help_dop853_license)

# Add informations about solver in global solverInfo-array.
push!(solverInfo,
  SolverInfo("dop853",
    "Runge-Kutta method of order 8(5,3) due to Dormand & Prince",
    tuple(:OPT_RTOL, :OPT_ATOL, 
          :OPT_OUTPUTMODE, :OPT_OUTPUTFCN, 
          :OPT_MAXSTEPS, :OPT_STEST, :OPT_EPS, :OPT_RHO, 
          :OPT_SSMINSEL, :OPT_SSMAXSEL, :OPT_SSBETA, 
          :OPT_MAXSS, :OPT_INITIALSS),
    tuple(
      SolverVariant("dop853_i64",
        "Dop853 with 64bit integers",
        DL_DOP853,
        tuple("dop853", "contd8")),
      SolverVariant("dop853_i32",
        "Dop853 with 32bit integers",
        DL_DOP853_I32,
        tuple("dop853", "contd8")),
    ),
    help_dop853_compile,
    help_dop853_license,
  )
)


# vim:syn=julia:cc=79:fdm=indent:
