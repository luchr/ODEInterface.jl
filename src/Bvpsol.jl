# Functions for Boundary-Value Solver: Bvpsol

"""Name for Loading bvpsol solver (64bit integers)."""
const DL_BVPSOL               = "bvpsol"

"""Name for Loading bvpsol solver (32bit integers)."""
const DL_BVPSOL_I32           = "bvpsol_i32"

"""macro for import bvp solver."""
macro import_bvpsol()
  :(
    using ODEInterface: bvpsol, bvpsol_i32
  )
end

"""macro for import bvpsol dynamic lib names."""
macro import_DLbvpsol()
  :(
    using ODEInterface: DL_BVPSOL, DL_BVPSOL_I32
  )
end

"""macro for importing Bvpsol help."""
macro import_bvpsol_help()
  :(
    using ODEInterface: help_bvpsol_compile, help_bvpsol_license
  )
end

@enum(ODE_SOLVER_USAGE,
      ODE_SOLVER_INTERNAL, ODE_SOLVER_JULIA)

@doc """
  The ode-solver to use for the initial value problems.
  """ -> ODE_SOLVER_USAGE

@doc """
  Use the DIFEX1-solver builtin in bvpsol.
  """ -> ODE_SOLVER_INTERNAL

@doc """
  Use the julia-function for solving the initial value problems.
  """ -> ODE_SOLVER_JULIA


"""
  bvpsol does not support "pass-through" arguments for FCN and BC.
  Hence we can only support one bvpsol-call at a time.
  """
bvpsol_global_cbi = nothing

"""
  Type encapsulating all required data for Bvpsol-Callbacks.
  
  Unfortunately bvpsol.f does not support passthrough arguments.
  
  We have the typical calling stack:

       bvpsol       
           ccall( BVPSOL_  ... )
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_bvpsolrhs                           │  ⎬ cb. rhs
              │    rhs                                    │  ⎪
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_bvpssolbc                           │  ⎬ cb. boundary
              │    bc                                     │  ⎪     conditions
              └───────────────────────────────────────────┘  ⎭
              ┌───────────────────────────────────────────┐  ⎫
              │unsafe_bvpsolivp                           │  ⎬ cb. solving
              │    odesolver(rhs,t,tEnd,x,opt)            │  ⎪ IVP
              └───────────────────────────────────────────┘  ⎭
  """
type BvpsolInternalCallInfos{FInt<:FortranInt, RHS_F<:Function, 
        BC_F<:Function, ODESOL_F<:Function} <: ODEinternalCallInfos
  logio        :: IO                    # where to log
  loglevel     :: UInt64                # log level
  # RHS:
  rhs          :: RHS_F                 # right-hand-side 
  rhs_mode     :: RHS_CALL_MODE         # how to call rhs
  rhs_lprefix  :: AbstractString        # saved log-prefix for rhs
  # BC:
  bc           :: BC_F                  # julia: functions for b-conditions
  bc_lprefix   :: AbstractString        # saved log-prefix for bc
  # problem specific
  N            :: FInt                  # Dimension of the problem
  # ODE-Solver
  odesol_usage :: ODE_SOLVER_USAGE      # what IVP-Solver to use?
  odesol_julia :: ODESOL_F              # a Julia-Function
  odeopt       :: OptionsODE            # Options for IVP-Solver
  ivp_lprefix  :: AbstractString        # saved log-prefix for ivp-call
end

type BvpsolArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
  FCN     :: Ptr{Void}         # rhs callback
  BC      :: Ptr{Void}         # boundary conditions
  IVPSOL  :: Ptr{Void}         # Initial Value Problem Solver
  N       :: Vector{FInt}      # Dimension
  M       :: Vector{FInt}      # Number of shooting nodes
  T       :: Vector{Float64}   # shooting nodes
  X       :: Matrix{Float64}   # start data at shooting nodes
  EPS     :: Vector{Float64}   # required rel. accuracy for solution
  IOPT    :: Vector{FInt}      # Option array
  INFO    :: Vector{FInt}      # Return code
  IRW     :: Vector{FInt}      # length of workspace RW
  RW      :: Vector{Float64}   # real workspace
  IIW     :: Vector{FInt}      # length of workspace IW
  IW      :: Vector{FInt}      # integer workspace

    ## Allow uninitialized construction
    (::Type{BvpsolArguments{T}}){T}() = new{T}()
end

"""
        function unsafe_bvpsolrhs{FInt<:FortranInt}(n_::Ptr{FInt}, 
                t_::Ptr{Float64}, x_::Ptr{Float64}, f_::Ptr{Float64})
  
  This is the right-hand side given as callback to bvpsol.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  
  uses hw1rhs
  """
function unsafe_bvpsolrhs{FInt<:FortranInt}(n_::Ptr{FInt}, 
        t_::Ptr{Float64}, x_::Ptr{Float64}, f_::Ptr{Float64})
  
  n = unsafe_load(n_); t = unsafe_load(t_)
  x = unsafe_wrap(Array, x_,(n,),false)
  f = unsafe_wrap(Array, f_,(n,),false)
  cbi = bvpsol_global_cbi :: BvpsolInternalCallInfos
  hw1rhs(n,t,x,f,cbi)
  return nothing
end

"""
        function unsafe_bvpsolrhs_c{FInt}(fint_flag::FInt)
  """
function unsafe_bvpsolrhs_c{FInt}(fint_flag::FInt)
  return cfunction(unsafe_bvpsolrhs, Void, (Ptr{FInt},Ptr{Float64},
    Ptr{Float64},Ptr{Float64}))
end

"""
        function bvpsolbc{CI}(xa,xb,r,cbi::CI)
  
  This function calls `bc` saved in `BvpsolInternalCallInfos`.
  """
function bvpsolbc{CI}(xa,xb,r,cbi::CI)
  lprefix = cbi.bc_lprefix
  
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_bc = l & LOG_BC > 0
  
  l_bc && println(lio,lprefix,"called with xa=",xa," xb=",xb)
  cbi.bc(xa,xb,r)
  l_bc && println(lio,lprefix,"bc result=",r)
  return nothing
end

"""
       function unsafe_bvpsolbc(xa_::Ptr{Float64}, xb_::Ptr{Float64}, 
         r_::Ptr{Float64}) -> nothing
  
  This is the callback for the boundary conditions given to bvpsol.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  
  uses bvpsolbc
  """
function unsafe_bvpsolbc(xa_::Ptr{Float64}, xb_::Ptr{Float64}, 
  r_::Ptr{Float64})

  cbi = bvpsol_global_cbi :: BvpsolInternalCallInfos
  n = cbi.N
  xa = unsafe_wrap(Array, xa_,(n,),false)
  xb = unsafe_wrap(Array, xb_,(n,),false)
  r  = unsafe_wrap(Array, r_ ,(n,),false)
  bvpsolbc(xa,xb,r,cbi)
  return nothing
end

"""
        function unsafe_bvpsolbc_c()
  """
function unsafe_bvpsolbc_c()
  return cfunction(unsafe_bvpsolbc, Void, 
        (Ptr{Float64},Ptr{Float64},Ptr{Float64}))
end

function bvpsolivp{FInt,CI}(t::Vector{Float64},
        x::Vector{Float64}, tend,tol,hmax,h::Vector{Float64},
        kflag::Vector{FInt}, cbi::CI)

  @assert cbi.odesol_usage == ODE_SOLVER_JULIA

  (lio,l)=(cbi.logio,cbi.loglevel)
  lprefix = cbi.ivp_lprefix
  l_ivp = l & LOG_BVPIVPSOL > 0
  
  opt = cbi.odeopt
  setOptions!(opt, 
    OPT_RTOL => tol, OPT_MAXSS => hmax, 
    OPT_INITIALSS => h[1], "KFLAG" => kflag[1])
  
  l_ivp && println(lio,lprefix,"calling ivp solver with rhs=",cbi.rhs,
    " tStart=",t[1]," tEnd=",tend," x=",x," opt=",opt)
  
  (ret_t, ret_x, ret_code, ret_stats) = 
    cbi.odesol_julia(cbi.rhs,t[1],tend,x,opt)
  l_ivp && println(lio,lprefix,"ivp solver returned ret_t=",ret_t,
    " ret_x=",ret_x," ret_code=",ret_code," ret_stats=",ret_stats)
  t[1]=ret_t
  x[:]=ret_x
  kflag[1] = ret_code
  if haskey(ret_stats,"step_predict")
    h[1] = ret_stats["step_predict"]
  else
    throw(OutputErrorODE(
      "IVP-solver had no field 'step_predict' in stats",cbi.odesol_julia))
  end
  if ret_t ≠ tend
    h[1] = 0
  end
  l_ivp && println(lio,lprefix,"returning h=",h[1]," kflag=",kflag[1])
  
  return nothing
end

"""
        function unsafe_bvpsolivp{FInt<:FortranInt}(n_::Ptr{FInt}, 
                fcn_::Ptr{Void}, t_::Ptr{Float64}, x_::Ptr{Float64}, 
                tend_::Ptr{Float64}, tol_::Ptr{Float64}, hmax_::Ptr{Float64}, 
                h_::Ptr{Float64}, kflag_::Ptr{FInt})

  This is the callback for bvpsol to solve initial value problems.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  
  uses bvpsolivp
  """
function unsafe_bvpsolivp{FInt<:FortranInt}(n_::Ptr{FInt}, 
        fcn_::Ptr{Void}, t_::Ptr{Float64}, x_::Ptr{Float64}, 
        tend_::Ptr{Float64}, tol_::Ptr{Float64}, hmax_::Ptr{Float64}, 
        h_::Ptr{Float64}, kflag_::Ptr{FInt})

  cbi = bvpsol_global_cbi::BvpsolInternalCallInfos
  n = cbi.N
  t = unsafe_wrap(Array, t_,(1,),false)
  tend = unsafe_load(tend_)
  x = unsafe_wrap(Array, x_,(n,),false)
  tol = unsafe_load(tol_); hmax = unsafe_load(hmax_); 
  h = unsafe_wrap(Array, h_,(1,),false)
  kflag = unsafe_wrap(Array, kflag_,(1,),false)
  bvpsolivp(t,x,tend,tol,hmax,h,kflag,cbi)
  return nothing
end

function unsafe_bvpsolivp_c{FInt}(fint_flag::FInt)
  return cfunction(unsafe_bvpsolivp, Void, 
    (Ptr{FInt},Ptr{Void},Ptr{Float64},
    Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Float64},
    Ptr{FInt}))
end


function bvpsol_ivp_dummy(rhs,t,tend,x,opt)
  throw(InternalErrorODE("bvpsol_ivp_dummy was called"))
end

"""
       function bvpsol(rhs::Function, bc::Function,
         t::Vector, x::Matrix, odesolver, opt::AbstractOptionsODE)
           -> (t,x,retcode,stats)
  
  The `bc` has to be a function of the following form:

       function bc(xa,xb,r) -> nothing
  
  It has to calculate the residual for the boundary conditions and save
  them in `r`.
  
  `t` is a Vector with all the multiple-shooting nodes.
  
  `x` gives the initial guess for all multiple-shooting nodes. Hence
  `size(x,2)==length(t)`.
  
  `odesolver`: Either `nothing`: then the internal solver of `bvpsol` is
  used. Or `odesolver` is a ode-solver (like `dopri5`, `dop853`, `seulex`, 
  etc.).

  `retcode` can have the following values:

        >0: computation successful: number of iterations
        -1:        Iteration stops at stationary point for OPT_SOLMETHOD==0
                   Gaussian elimination failed due to singular 
                   Jacobian for OPT_SOLMETHOD==1
        -2: Iteration stops after OPT_MAXSTEPS 
        -3: Integrator failed to complete the trajectory
        -4: Gauss Newton method failed to converge
        -5: Given initial values inconsistent with separable linear bc
        -6:        Iterative refinement faild to converge for OPT_SOLMETHOD==0
                   Termination since multiple shooting condition or
                   condition of Jacobian is too bad for OPT_SOLMETHOD==1
        -7: wrong EPS (should not happen; checked by ODEInterface module)
        -8: Condensing algorithm for linear block system fails, try
            OPT_SOLMETHOD==1
        -9: Sparse linear solver failed
       -10: Real or integer work-space exhausted
       -11: Rank reduction failed - resulting rank is zero
 
  In `opt` the following options are used:
  
      ╔═════════════════╤══════════════════════════════════════════╤═════════╗
      ║  Option OPT_…   │ Description                              │ Default ║
      ╠═════════════════╪══════════════════════════════════════════╪═════════╣
      ║ RTOL            │ relative accuracy for soltuion           │    1e-6 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSTEPS        │ maximum permitted number of iteration    │      40 ║
      ║                 │ steps                                    │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ BVPCLASS        │ boundary value problem classification:   │       2 ║
      ║                 │ 0: linear                                │         ║
      ║                 │ 1: nonlinear with good initial data      │         ║
      ║                 │ 2: highly nonlinear & bad initial data   │         ║
      ║                 │ 3: highly nonlinear & bad initial data & │         ║
      ║                 │    initial rank reduction to separable   │         ║
      ║                 │    linear boundary conditions            │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SOLMETHOD       │ switch for solution method               │       0 ║
      ║                 │ 0: use local linear solver with          │         ║
      ║                 │    condensing algorithm                  │         ║
      ║                 │ 1: use global sparse linear solver       │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ IVPOPT          │ An OptionsODE-object with the options    │ empty   ║
      ║                 │ for the solver of the initial value      │ options ║
      ║                 │ problem.                                 │         ║
      ║                 │ In this OptionsODE-object bvpsol changes │         ║
      ║                 │ OPT_MAXSS, OPT_INITIALSS, OPT_RTOL       │         ║
      ║                 │ to give the IVP-solver solution hints.   │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ RHS_CALLMODE    │ see help_callsolvers()                   │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝
  """
function bvpsol(rhs::Function, bc::Function,
  t::Vector, x::Matrix, odesolver, opt::AbstractOptionsODE)
  return bvpsol_impl(rhs,bc,t,x,odesolver,opt,BvpsolArguments{Int64}())
end

"""
  bvpsol with 32bit integers, see bvpsol.
  """
function bvpsol_i32(rhs::Function, bc::Function,
  t::Vector, x::Matrix, odesolver, opt::AbstractOptionsODE)
  return bvpsol_impl(rhs,bc,t,x,odesolver,opt,BvpsolArguments{Int32}()) 
end

function bvpsol_impl{FInt<:FortranInt}(rhs::Function, bc::Function,
  t::Vector, x::Matrix, odesolver, 
  opt::AbstractOptionsODE, args::BvpsolArguments{FInt})

  (lio,l,l_g,l_solver,lprefix) = solver_init("bvpsol",opt)

  if l_g
    println(lio,lprefix,
      "called with rhs=",rhs," t=",t," x=",x," and opt")
    show(lio,opt); println(lio)
  end

  (method_bvpsol, method_bldfx1) = getAllMethodPtrs(
     (FInt == Int64)? DL_BVPSOL : DL_BVPSOL_I32 )
  
  d = FInt(0)
  try
    d = convert(FInt,size(x,1))
  catch e
    throw(ArgumentErrorODE("Cannot convert size(x,1) to $FInt",:x,e))
  end
  0==d && throw(ArgumentErrorODE("x has no rows",:x))
  N=d; args.N = [d]
  
  m = FInt(0)
  try
    m = convert(FInt,size(x,2))
  catch e
    throw(ArgumentErrorODE("Cannot convert size(x,2) to $FInt",:x,e))
  end
  m<2 && throw(ArgumentErrorODE("x has less than 2 columns",:x))
  try
    @assert length(t) == m
  catch e
    throw(ArgumentErrorODE("size(x,2) ≠ length(t)",:t,e))
  end
  M = m; args.M = [m]
  
  try
    args.T = getVectorCheckLength(t,Float64,m);
  catch e
    throw(ArgumentErrorODE("cannot convert t to Vector{Float64}",:t,e))
  end
  
  try
    args.X = getMatrixCheckSize(x,Float64,d,m)
  catch e
    throw(ArgumentErrorODE("cannot convert x to Matrix{Float64}",:x,e))
  end

  rhs_mode = solver_extract_rhsMode(opt)
  rhs_lprefix = "unsafe_bvpsolrhs: "
  bc_lprefix = "unsafe_bvpsolbc: "
  
  try
    args.EPS = [ convert(Float64,getOption(opt,OPT_RTOL,1e-6)) ]
    @assert 0 < args.EPS[1] <= 1.0e-2
  catch e
    throw(ArgumentErrorODE("Cannot convert OPT_RTOL to Float64",:opt,e))
  end
  
  args.IOPT = Vector{FInt}(6)
  ivpopt = nothing
  ivp_lprefix = "unsafe_bvpsolivp: "
  OPT = nothing
  try
    OPT = OPT_MAXSTEPS; args.IOPT[1] = convert(FInt,getOption(opt,OPT,40))
    @assert args.IOPT[1] > 0
    
    OPT = OPT_BVPCLASS; args.IOPT[2] = convert(FInt,getOption(opt,OPT,2))
    @assert 0 ≤ args.IOPT[2] ≤ 3
    
    OPT = OPT_SOLMETHOD; args.IOPT[3] = convert(FInt,getOption(opt,OPT,0))
    @assert 0 ≤ args.IOPT[3] ≤ 1
    
    args.IOPT[4] = -1; args.IOPT[5] = 0;
    
    OPT = OPT_IVPOPT; ivpopt = getOption(opt,OPT,nothing)
    if ivpopt == nothing
      ivpopt = OptionsODE("Autogenerated by bvpsol")
    end
    @assert isa(ivpopt,OptionsODE)
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end
  
  if     args.IOPT[3]==0
    args.IIW = [ FInt(4*N + 2*N*N) ]
    args.IRW = [ FInt(N*N*(M+5) + 10*M*N + 10*N +M) ]
  elseif args.IOPT[3]==1
    NZ = N*N*(M+1)*N*(M-1)
    LICN = 2*NZ
    h1 = FInt(ceil(1.5*NZ))
    h2 = NZ+4*M*N
    LIRN = min(max(h1,h2),LICN)
    args.IIW = [ FInt(2*N*N + 3*N + 16*M*N + 2*NZ + LICN + LIRN) ]
    args.IRW = [ FInt(N*N*(M+1) + 12*N*M + 4*N + M -1 + LICN) ]
  end
  args.RW = zeros(Float64,args.IRW[1])
  args.IW = zeros(FInt,args.IIW[1])
  args.FCN = unsafe_bvpsolrhs_c(FInt(0))
  args.BC = unsafe_bvpsolbc_c()
  args.INFO = zeros(FInt,1)

  try
    @assert (odesolver == nothing) || isa(odesolver,Function)
  catch e
    throw(ArgumentErrorODE("odesolver not valid",:odesolver,e))
  end
  if odesolver==nothing && d>1024
    throw(ArgumentErrorODE(string("internal solver of bvpsol cannot handle ",
      "more than 1024 equations")))
  end
  if odesolver == nothing
    args.IVPSOL = method_bldfx1
  else
    args.IVPSOL = unsafe_bvpsolivp_c(FInt(0))
  end

  if bvpsol_global_cbi ≠ nothing
    throw(ArgumentErrorODE(string("The Fortran solver bvpsol does not ",
      "support 'pass-through' arguments. Hence this julia module does not ",
      "support concurrent/nested bvpsol calls. Sorry.")))
  end
  
  try
    global bvpsol_global_cbi = BvpsolInternalCallInfos(lio,l,rhs,rhs_mode,
        rhs_lprefix, bc,bc_lprefix,N, 
        (odesolver==nothing)?ODE_SOLVER_INTERNAL : ODE_SOLVER_JULIA,
        (odesolver==nothing)?bvpsol_ivp_dummy : odesolver,
        ivpopt, ivp_lprefix)
    
    if l_solver
      println(lio,lprefix,"call Fortran-bvpsol $method_bvpsol with")
      dump(lio,args)
    end
    
    ccall( method_bvpsol, Void,
      (Ptr{Void}, Ptr{Void}, Ptr{Void},           # Rightsidefunc, BC, IVPSOL
       Ptr{FInt}, Ptr{FInt},                      # N, M 
       Ptr{Float64}, Ptr{Float64}, Ptr{Float64},  # t, x, EPS
       Ptr{FInt}, Ptr{FInt},                      # IOPT, INFO
       Ptr{FInt}, Ptr{Float64},                   # IRW, RW
       Ptr{FInt}, Ptr{FInt},                      # IIW, IW
      ),
      args.FCN, args.BC, args.IVPSOL,
      args.N, args.M,
      args.T, args.X, args.EPS,
      args.IOPT, args.INFO,
      args.IRW, args.RW,
      args.IIW, args.IW,
    )
    
    if l_solver
      println(lio,lprefix,"Fortran-bvpsol $method_bvpsol returned")
      dump(lio,args)
    end
  finally
    global bvpsol_global_cbi = nothing
  end
  l_g && println(lio,lprefix,string("done INFO=",args.INFO[1]))
  stats = Dict{AbstractString,Any}()
  return (args.T, args.X, args.INFO[1], stats)
end


"""  
  ## Compile BVPSOL

  The Fortran source code can be found at:
  
       http://elib.zib.de/pub/elib/codelib/bvpsol/
  
  See `help_bvpsol_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile BVPSOL with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o bvpsol.o bvpsol.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o linalg_bvpsol.o linalg_bvpsol.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o zibconst.o zibconst.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o ma28_bvpsol.o ma28_bvpsol.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o bvpsol.so 
                bvpsol.o linalg_bvpsol.o zibconst.o ma28_bvpsol.o
       gfortran -shared -fPIC -o bvpsol.dylib
                bvpsol.o linalg_bvpsol.o zibconst.o ma28_bvpsol.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile BVPSOL with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o bvpsol.o bvpsol.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o linalg_bvpsol.o linalg_bvpsol.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o zibconst.o zibconst.f
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o ma28_bvpsol.o ma28_bvpsol.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o bvpsol.so 
                bvpsol.o linalg_bvpsol.o zibconst.o ma28_bvpsol.o
  
  ### Using `gfortran` and 32bit integers (Linux and Mac)
  
  Here is an example how to compile BVPSOL with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o bvpsol_i32.o bvpsol.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o linalg_bvpsol_i32.o linalg_bvpsol.f
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o zibonst_i32.o zibconst.f 
       gfortran -c -fPIC -fdefault-real-8 -fdefault-double-8 
                -o ma28_bvpsol_i32.o ma28_bvpsol.f 
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):
  
       gfortran -shared -fPIC -o bvpsol_i32.so 
                 bvpsol_i32.o linalg_bvpsol_i32.o zibconst_i32.o ma28_bvpsol_i32.o
       gfortran -shared -fPIC -o bvpsol_i32.dylib
                 bvpsol_i32.o linalg_bvpsol_i32.o zibconst_i32.o ma28_bvpsol_i32.o
  
  ### Using `gfortran` and 32bit integers (Windows)
  
  Here is an example how to compile BVPSOL with `Float64` reals and
  `Int32` integers with `gfortran`:
  
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o bvpsol_i32.o bvpsol.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o linalg_bvpsol_i32.o linalg_bvpsol.f
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o zibconst_i32.o zibconst.f 
       gfortran -c -fdefault-real-8 -fdefault-double-8 
                -o ma28_bvpsol_i32.o ma28_bvpsol.f 
  
  In order to get create a shared library (from the object file above) use:

       gfortran -shared -o bvpsol_i32.dll
                 bvpsol_i32.o linalg_bvplsol_i32.o zibconst_i32.o
                 ma28_bvpsol_i32.o
  
  """
function help_bvpsol_compile()
  return Docs.doc(help_bvpsol_compile)
end

"""
  # Licence
  You may use or modify this code for your own non commercial
  purposes for an unlimited time.
  This code can be distributed in unmodified form as part of the
  complete ODEInterface packet. Apart from that, you should not
  deliver it without a special permission of the Zuse Institute
  Berlin (ZIB).
  In case you intend to use the code commercially, we oblige you
  to sign an according licence agreement with ZIB.
  
  # Warranty
  This code has been tested up to a certain level. Defects and
  weaknesses, which may be included in the code, do not establish
  any warranties by ZIB. ZIB does not take over any liabilities
  which may follow from aquisition or application of this code.
  """
function help_bvpsol_license()
  return Docs.doc(help_bvpsol_license)
end


# Add informations about solvers in global solverInfo-array.
push!(solverInfo,
  SolverInfo("bvpsol",
    "Boundary Value Problem Solver for highly nonlinear problems",
    tuple( :OPT_RTOL, :OPT_MAXSTEPS, :OPT_BVPCLASS, :OPT_SOLMETHOD,
           :OPT_RHS_CALLMODE, :OPT_IVPOPT
          ),
    tuple(
      SolverVariant("bvpsol_i64",
        "Bvpsol with 64bit integers",
        DL_BVPSOL,
        tuple("bvpsol","bldfx1")),
      SolverVariant("bvpsol_i32",
        "Bvplsolwith 32bit integers",
        DL_BVPSOL_I32,
        tuple("bvpsol","bldfx1")),
    ),
    help_bvpsol_compile,
    help_bvpsol_license,
  )
)

# vim:syn=julia:cc=79:fdm=indent:
