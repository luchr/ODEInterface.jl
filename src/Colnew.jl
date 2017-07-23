# Functions for Boundary-Value Solver: Colsys/Colnew

"""Name for Loading colnew solver (64bit integers)."""
const DL_COLNEW               = "colnew"

"""Name for Loading colnew solver (32bit integers)."""
const DL_COLNEW_I32           = "colnew_i32"

"""macro for import colnew solver."""
macro import_colnew()
  :(
    using ODEInterface: colnew, colnew_i32, evalSolution
  )
end

"""macro for import colnew dynamic lib names."""
macro import_DLcolnew()
  :(
    using ODEInterface: DL_COLNEW, DL_COLNEW_I32
  )
end

"""macro for importing Colnew help."""
macro import_colnew_help()
  :(
    using ODEInterface: help_colnew_compile, help_colnew_license
  )
end

"""
  colnew does not support "pass-through" arguments for FSUB, DFSUB,
  GSUB, DGSUB, GUES.
  Hence we can only support one colnew-call at a time.
  """
colnew_global_cbi = nothing

"""
  Type encapsulating all required data for colnew-Callbacks.
  
  Unfortunately colnew.f does not support passthrough arguments.
  
  We have the typical calling stack:

  """
type ColnewInternalCallInfos{FInt<:FortranInt, 
        RHS_F<:Function, DRHS_F<:Function,
        BC_F<:Function, DBC_F<:Function,
        GUESS_F<:Function} <: ODEinternalCallInfos
  logio        :: IO                    # where to log
  loglevel     :: UInt64                # log level
  n            :: FInt                  # number of ODEs
  d            :: FInt                  # dimension of (1st-order) system
  # RHS:
  rhs          :: RHS_F                 # right-hand-side of ODEs
  rhs_lprefix  :: AbstractString        # saved log-prefix for rhs
  rhs_count    :: Int                   # count: number of calls to rhs
  # DRHS:
  Drhs         :: DRHS_F                # Jacobian of rhs
  Drhs_lprefix :: AbstractString        # saved log-prefix for Drhs
  Drhs_count   :: Int                   # count: number of calls to Drhs
  # BC:
  bc           :: BC_F                  # side conditions (at ζ)
  bc_lprefix   :: AbstractString        # saved log-prefix for bc
  bc_count     :: Int                   # count: number of calls to bc
  # DBC:
  Dbc          :: DBC_F                 # Jacobian of side conditions
  Dbc_lprefix  :: AbstractString        # saved log-prefix for Dbc
  Dbc_count    :: Int                   # count: number of calls to Dbc
  # guess:
  guess        :: GUESS_F               # Function for guess
  guess_lprefix:: AbstractString        # saved log-prefix for guess
end

type ColnewArguments{FInt<:FortranInt} <: AbstractArgumentsODESolver{FInt}
  NCOMP   :: Vector{FInt}      # number of diff. equations
  M       :: Vector{FInt}      # order of the j-th diff. eq.
  ALEFT   :: Vector{Float64}   # left end of (time-)interval
  ARIGHT  :: Vector{Float64}   # right end of (time-)interval
  ZETA    :: Vector{Float64}   # j-th side condition point
  IPAR    :: Vector{FInt}      # (integer) parameters
  LTOL    :: Vector{FInt}      # assign TOL[j] to component of Z(U)
  TOL     :: Vector{Float64}   # tolerances
  FIXPNT  :: Vector{Float64}   # points to be included in all meshes
  ISPACE  :: Vector{FInt}      # integer work space
  FSPACE  :: Vector{Float64}   # floating-point work space
  IFLAG   :: Vector{FInt}      # return-code of colnew
    ## Allow uninitialized construction
  @WHEREFUNC(FInt,
  function ColnewArguments{FInt}(dummy::FInt)
    return new{FInt}()
  end
  )
end

type ColnewSolution{FInt<:FortranInt} <: AbstractODESolution{FInt}
  method_appsln :: Ptr{Void}   # Ptr to Fortran-method for solution eval
  t_a     :: Float64           # a
  t_b     :: Float64           # b
  n       :: FInt              # number of ODEs
  d       :: FInt              # dimension of (1st-order) system
  ISPACE  :: Vector{FInt}      # part of integer work-space
  FSPACE  :: Vector{Float64}   # solution coefficients extracted from FSPACE
end


"""
  This function calls `rhs` saved in ColnewInternalCallInfos.
  """
function colnew_rhs{CI}(t, z, f, cbi::CI)
  lprefix = cbi.rhs_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_rhs = l & LOG_RHS > 0

  cbi.rhs_count += 1

  l_rhs && println(lio, lprefix, "called with t=", t, " z=", z)
  cbi.rhs(t, z, f)
  l_rhs && println(lio, lprefix, "rhs result=", f)
  return nothing
end

"""
        function unsafe_colnew_rhs(t_::Ptr{Float64}, z_::Ptr{Float64}, 
          f_::Ptr{Float64})

  This is the right-hand side given as callback to colnew.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_colnew_rhs(t_::Ptr{Float64}, z_::Ptr{Float64}, 
  f_::Ptr{Float64})

  cbi = colnew_global_cbi :: ColnewInternalCallInfos
  n = cbi.n; d = cbi.d
  
  t = unsafe_load(t_)
  z = unsafe_wrap(Array, z_, (d,), false)
  f = unsafe_wrap(Array, f_, (n,), false)

  colnew_rhs(t, z, f, cbi)
  return nothing
end

function unsafe_colnew_rhs_c{FInt}(fint_flag::FInt)
  return cfunction(unsafe_colnew_rhs, Void, 
    (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}))
end

"""
  This function calls `Drhs` saved in ColnewInternalCallInfos.
  """
function colnew_Drhs{CI}(t, z, df, cbi::CI)
  lprefix = cbi.Drhs_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_Drhs = l & LOG_JAC > 0

  cbi.Drhs_count += 1

  l_Drhs && println(lio, lprefix, "called with t=", t, " z=", z)
  cbi.Drhs(t, z, df)
  l_Drhs && println(lio, lprefix, "Drhs result=", df)
  return nothing
end

"""
        function unsafe_colnew_Drhs(t_::Ptr{Float64}, z_::Ptr{Float64},
          df_::Ptr{Float64})

  This is the Jacobian for the right-hand side given as callback to colnew.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_colnew_Drhs(t_::Ptr{Float64}, z_::Ptr{Float64},
  df_::Ptr{Float64})

  cbi = colnew_global_cbi :: ColnewInternalCallInfos
  n = cbi.n; d = cbi.d

  t = unsafe_load(t_)
  z = unsafe_wrap(Array, z_, (d,), false)
  df = unsafe_wrap(Array, df_, (n,d,), false)

  colnew_Drhs(t, z, df, cbi)
  return nothing
end

function unsafe_colnew_Drhs_c{FInt}(fint_flag::FInt)
  return cfunction(unsafe_colnew_Drhs, Void,
    (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}))
end

"""
  This function calls `bc` saved in ColnewInternalCallInfos.
  """
function colnew_bc{CI}(i, z, bc, cbi::CI)
  lprefix = cbi.bc_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_bc = l & LOG_BC > 0

  cbi.bc_count += 1

  l_bc && println(lio, lprefix, "called with i=", i, " z=", z)
  cbi.bc(i, z, bc)
  l_bc && println(lio, lprefix, "bc result=", bc)
  return nothing
end

"""
        function unsafe_colnew_bc{FInt<:FortranInt}(i_::Ptr{FInt},
          z_::Ptr{Float64}, bc_::Ptr{Float64})

  This is the side-/boundary-conditions given as callback
  to colnew.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_colnew_bc{FInt<:FortranInt}(i_::Ptr{FInt},
  z_::Ptr{Float64}, bc_::Ptr{Float64})

  cbi = colnew_global_cbi :: ColnewInternalCallInfos
  d = cbi.d

  i = unsafe_load(i_)
  z = unsafe_wrap(Array, z_, (d,), false)
  bc = unsafe_wrap(Array, bc_, (1,), false)

  colnew_bc(i, z, bc, cbi)
  return nothing
end

function unsafe_colnew_bc_c{FInt}(fint_flag::FInt)
  return cfunction(unsafe_colnew_bc, Void,
    (Ptr{FInt}, Ptr{Float64}, Ptr{Float64}))
end

"""
  This function calls `Dbc` saved in ColnewInternalCallInfos.
  """
function colnew_Dbc{CI}(i, z, dbc, cbi::CI)
  lprefix = cbi.Dbc_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_Dbc = l & LOG_JACBC > 0

  cbi.Dbc_count += 1

  l_Dbc && println(lio, lprefix, "called with i=", i, " z=", z)
  cbi.Dbc(i, z, dbc)
  l_Dbc && println(lio, lprefix, "Dbc result=", dbc)
  return nothing
end

"""
        function unsafe_colnew_Dbc{FInt<:FortranInt}(i_::Ptr{FInt},
          z_::Ptr{Float64}, Dbc_::Ptr{Float64})

  This is the jacobian for the side-/boundary-conditions given as callback
  to colnew.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_colnew_Dbc{FInt<:FortranInt}(i_::Ptr{FInt},
  z_::Ptr{Float64}, Dbc_::Ptr{Float64})

  cbi = colnew_global_cbi :: ColnewInternalCallInfos
  d = cbi.d

  i = unsafe_load(i_)
  z = unsafe_wrap(Array, z_, (d,), false)
  Dbc = unsafe_wrap(Array, Dbc_, (d,), false)

  colnew_Dbc(i, z, Dbc, cbi)
  return nothing
end

function unsafe_colnew_Dbc_c{FInt}(fint_flag::FInt)
  return cfunction(unsafe_colnew_Dbc, Void,
    (Ptr{FInt}, Ptr{Float64}, Ptr{Float64}))
end

"""
  This function calls `guess` saved in ColnewInternalCallInfos.
  """
function colnew_guess{CI}(t, z, dxm, cbi::CI)
  lprefix = cbi.guess_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_guess = l & LOG_GUESS > 0

  l_guess && println(lio, lprefix, "called with t=", t)
  cbi.guess(t, z, dxm)
  l_guess && println(lio, lprefix, "guess z-result=", z, "dxm-result=", dxm)
  return nothing
end

"""
  This is the guess function given as callback to colnew.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_colnew_guess(t_::Ptr{Float64}, z_::Ptr{Float64},
  dmx_::Ptr{Float64})

  cbi = colnew_global_cbi :: ColnewInternalCallInfos
  n = cbi.n; d = cbi.d

  t = unsafe_load(t_)
  z = unsafe_wrap(Array, z_, (d,), false)
  dmx = unsafe_wrap(Array, dmx_, (n,), false)

  colnew_guess(t, z, dmx, cbi)
  return nothing
end

function unsafe_colnew_guess_c{FInt}(fint_flag::FInt)
  return cfunction(unsafe_colnew_guess, Void,
    (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}))
end


"""
        function colnew(interval::Vector, orders::Vector, ζ::Vector,
          rhs::Function, Drhs::Function,
          bc::Function, Dbc::Function, guess, opt::AbstractOptionsODE)

  Solve multi-point boundary value problem with colnew.

  ζ∊ℝᵈ with a ≤ ζ(1)=ζ₁ ≤ ζ(2)=ζ₂ ≤ ⋯ ≤ ζ(d) ≤ b are the (time-)points
  were side/boundary conditions are given:

        ├─────┼─────┼─────────┼─....───┼─────┤
       t=a  t=ζ(1) t=ζ(2)    t=ζ(3)  t=ζ(d)  t=b

  for the n ODEs
         ∂xᵢ
        ──────  = xᵢ⁽ᵐ⁽ⁱ⁾⁾ = fᵢ(t, z(x(t))          (i=1,2,…,n)   [*]
        ∂tᵐ⁽ⁱ⁾

  where the i-th ODE has order m(i). [x(t)∊ℝⁿ].

  z is the transformation to first order:
  z(x(t))∊ℝᵈ is the "first-order" state one gets if the n ODEs [*]
  are transformed to a first-order system:

       z(x(t)) = ( x₁(t), x₁'(t), x₁''(t), …, x₁⁽ᵐ⁽¹⁾⁻¹⁾,
                   x₂(t), x₂'(t), x₂''(t), …, x₂⁽ᵐ⁽²⁾⁻¹⁾,
                   ⋯                                    ,
                   xₙ(t), xₙ'(t), xₙ''(t), …, xₙ⁽ᵐ⁽ⁿ⁾⁻¹⁾  )

  Hence one has the requirement: ∑m(i) = d.

  The boundary-/side-conditions at the points ζ(j) are given in the form

       bcⱼ(ζⱼ, z(x(ζⱼ))) = 0                         (j=1,2,…,d)


  Restrictions (in the colnew code):
  * at maximum 20 ODEs: n ≤ 20
  * at maximum 40 dimensions: d ≤ 40
  * The orders m(i) have to satisfy: 1 ≤ m(i) ≤ 4   for all i=1,2,…,n.

  All (Julia-)callback-functions (like rhs, etc.) use the in-situ call-mode,
  i.e. they have to write the result in an preallocated vector.

  ## rhs

  `rhs` must be a function of the form

      function rhs(t, z, f)

  with the input data: t (scalar) time and z∈ℝᵈ (z=z(x(t))).
  The values of the right-hand side have to be saved in f: f∈ℝⁿ! 
  Only the non-trivial parts of the right-hand side must be calculated.

  ## Drhs

  `Drhs` must be a function of the form

      function Drhs(t, z, df)

  with the input data: t (scalar) time and z∈ℝᵈ (z=z(x(t))).
  The values of the jacobian of the right-hand side have to be saved
  in df: df∈ℝⁿˣᵈ!

                 ∂fᵢ
      df(i,j) = ─────      (i=1,…,n;  j=1,…,d)
                 ∂zⱼ

  ## bc

  `bc` must be a function of the form

      function bc(i, z, bc)

  with the input data: integer index i and z∈ℝᵈ (z=z(x(t))).
  The scalar(!) value of the i-th side-condition (at time ζ(i)) has to
  be saved in bc, which is a vector of length 1.

  ## Dbc

  `Dbc` must be a function of the form

      function Dbc(i, z, dbc)

  with the input data: integer index i and z∈ℝᵈ (z=z(x(t))).
  The  values of the derivative of the i-th side-condition 
  (at time ζ(i)) has to be saved in dbc:

                ∂bcᵢ
      dbc(j) = ─────      (j=1,…,d)
                ∂zⱼ

  ## guess

  `guess` must be function of the form

      function guess(t, z, dmx)

  with the input data t∈[a,b]. Guesses are needed for the following
  values: z=z(x(t))∈ℝᵈ and

                ∂xᵢ
      dmx(i) = ────────      (i=1,…,n)
                ∂tᵐ⁽ⁱ⁾

  In `opt` the following options are used:
  
      ╔═════════════════╤══════════════════════════════════════════╤═════════╗
      ║  Option OPT_…   │ Description                              │ Default ║
      ╠═════════════════╪══════════════════════════════════════════╪═════════╣
      ║ BVPCLASS        │ boundary value problem classification:   │       1 ║
      ║                 │ 0: linear                                │         ║
      ║                 │ 1: nonlinear and regular                 │         ║
      ║                 │ 2: nonlinear and "extra sensitive"       │         ║
      ║                 │    (first relax factor is rstart and the │         ║
      ║                 │    nonlinear iteration does not rely     │         ║
      ║                 │    on past convergence)                  │         ║
      ║                 │ 3: fail-early: return immediately upon   │         ║
      ║                 │    (a) two successive non-convergences   │         ║
      ║                 │        or                                │         ║
      ║                 │    (b) after obtaining an error estimate │         ║
      ║                 │        for the first time                │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ RTOL            │ relative *and* absolute accuracy for     │    1e-6 ║
      ║                 │ solution. Must be a vector of length d.  │         ║
      ║                 │ If a scalar is given (like the default   │         ║
      ║                 │ value of 1e-6) then the vector           │         ║
      ║                 │    RTOL*ones(Float64, d)                 │         ║
      ║                 │ is generated.                            │         ║
      ║                 │ Some entries can be NaN! If an entry     │         ║
      ║                 │ is NaN, then no error checking is done   │         ║
      ║                 │ for this component.                      │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ COLLOCATIONPTS  │ number (=k) of collocation points per    │ see left║
      ║                 │ sub-interval.                            │         ║
      ║                 │ Requirement:                             │         ║
      ║                 │   orders[i] ≤ k ≤ 7                      │         ║
      ║                 │ Default:                                 │         ║
      ║                 │   k = max( max(orders)+1, 5-max(orders) )│         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SUBINTERVALS    │ Either a positive integer scalar or a    │       5 ║
      ║                 │ vector of (Float)-times:                 │         ║
      ║                 │                                          │         ║
      ║                 │ (a) scalar: use a "uniform-like" initial │         ║
      ║                 │ grid with the given integer as number    │         ║
      ║                 │ of subintervals.                         │         ║
      ║                 │ Why "uniform-like" and not "uniform"?    │         ║
      ║                 │ Because all values of ζ and all values of│         ║
      ║                 │ OPT_ADDGRIDPOINTS have to be in the grid.│         ║
      ║                 │ If the scalar is too small for all this  │         ║
      ║                 │ values it is increased (internally).     │         ║
      ║                 │                                          │         ║
      ║                 │ (b) vector: all points must be inside    │         ║
      ║                 │ the interval (a,b). Then this points     │         ║
      ║                 │ are used as initial grid. Values of ζ,   │         ║
      ║                 │ OPT_ADDGRIDPOINTS and a and b are added  │         ║
      ║                 │ automatically by this interface.         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ FREEZEINTERVALS │ Only used if OPT_SUBINTERVALS is a       │   false ║
      ║                 │ vector. In this case this flags indicates│         ║
      ║                 │ if colnew is allowed to adaptively       │         ║
      ║                 │ change the grid.                         │         ║
      ║                 │ If true, all grid adaption is turned off │         ║
      ║                 │ and no mesh selection is done.           │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ MAXSUBINTERVALS │ number of maximal subintervals.          │      50 ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ DIAGNOSTICOUTPUT│ diagnostic output for colnew:            │       1 ║
      ║                 │   -1 : full diagnostic printout          │         ║
      ║                 │    0 : selected printout                 │         ║
      ║                 │    1 : no printout                       │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ ADDGRIDPOINTS   │ additional points that are always added  │      [] ║
      ║                 │ to every (time-)grid.                    │         ║
      ║                 │ Every grid contains all values in ζ and  │         ║
      ║                 │ the values in the interval argument.     │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝

  """
function colnew(interval::Vector, orders::Vector, ζ::Vector,
  rhs::Function, Drhs::Function,
  bc::Function, Dbc::Function, guess, opt::AbstractOptionsODE)

  return colnew_impl(interval, orders, ζ, rhs, Drhs, bc, Dbc, guess, opt, 
    ColnewArguments{Int64}(Int64(0)))
end

"""
  colnew with 32bit integers, see colnew.
  """
function colnew_i32(interval::Vector, orders::Vector, ζ::Vector,
  rhs::Function, Drhs::Function,
  bc::Function, Dbc::Function, guess, opt::AbstractOptionsODE)

  return colnew_impl(interval, orders, ζ, rhs, Drhs, bc, Dbc, guess, opt, 
    ColnewArguments{Int32}(Int32(0)))
end

function colnew_impl{FInt<:FortranInt}(
  interval::Vector, orders::Vector, ζ::Vector,
  rhs::Function, Drhs::Function,
  bc::Function, Dbc::Function, guess, opt::AbstractOptionsODE,
  args::ColnewArguments{FInt})

  (lio,l,l_g,l_solver,lprefix) = solver_init("colnew",opt)

  if l_g
    println(lio,lprefix,
      "called with interval=",interval," rhs=", rhs,
      " Drhs=", Drhs, " bc=", bc, " Dbc=", Dbc, " guess=", guess, " and opt")
    show(lio,opt); println(lio)
  end

  (method_colnew, method_appsln) = getAllMethodPtrs(
    (FInt == Int64)? DL_COLNEW : DL_COLNEW_I32 )

  # (time-)vector with interval [a,b]
  t_ab = nothing
  try
    t_ab = getVectorCheckLength(interval, Float64, 2)
  catch e
    throw(ArgumentErrorODE(
      "cannot convert interval to Vector{Float64} with 2 components", 
      :interval, e))
  end
  t_ab[1] ≥ t_ab[2] && throw(ArgumentErrorODE(
    "interval=[a,b] with a<b needed", :interval))
  args.ALEFT = [ t_ab[1] ]
  args.ARIGHT = [ t_ab[2] ]

  # n = NCOMP
  n = length(orders)
  n ≤ 0 && throw(ArgumentErrorODE(
    "orders must be a vector with length ≥ 1", :orders))
  n > 20 && throw(ArgumentErrorODE(
    string("restriction: number of ODEs must be ≤ 20, but ",
      "orders vector had length ",n), :orders))
  args.NCOMP = [ FInt(n) ]

  # orders = M
  try
    orders = getVectorCheckLength(orders, FInt, n) # always copy vector
  catch e
   throw(ArgumentErrorODE(
     "cannot convert orders to Vector{Fortran Ints}", :orders, e))
  end
  any(orders.<1) && throw(ArgumentErrorODE(
    "all components of the orders vector must be positive", :orders))
  any(orders.>5) && throw(ArgumentErrorODE(
    "all components of the orders vector must be ≤ 4", :orders))
  order_max = maximum(orders)
  d = sum(orders)
  d > 40 && throw(ArgumentErrorODE(
    "the dimension=sum(orders) must be ≤ 40", :orders))
  args.M = orders

  # ζ vector
  try
    ζ = getVectorCheckLength(ζ, Float64, d)  # always copy vector
  catch e
    throw(ArgumentErrorODE(
      string("cannot convert ζ to Vector{Float64} of length ", d),
      :ζ, e))
  end
  if !( all(ζ .≥ t_ab[1]) && all(ζ .≤ t_ab[2]) )
    throw(ArgumentErrorODE(
      "all side condition points in ζ must be inside of interval",
      :ζ))
  end
  if !( all( diff(ζ) .≥ 0.0 ) )
    throw(ArgumentErrorODE("ζ must be ascending", ζ))
  end
  args.ZETA = ζ

  max_subintervals = 20
  subintervals = nothing
  # ipar
  args.IPAR = zeros(FInt, 11)
  OPT = nothing
  try
    # RTOL ( => ltol, tol)
    OPT = OPT_RTOL
    rtol = getOption(opt, OPT, 1e-6)
    if isscalar(rtol) 
      rtol = rtol*ones(Float64, d)
    end
    rtol = getVectorCheckLength(rtol, Float64, d, false)
    args.LTOL = zeros(FInt, 0)
    args.TOL = zeros(Float64, 0)
    for j=1:d
      if !isnan(rtol[j])
        push!(args.LTOL, FInt(j))
        push!(args.TOL, rtol[j])
      end
    end
    @assert length(args.LTOL) == length(args.TOL)
    length(args.LTOL)==0 && throw(ArgumentErrorODE(
      "All components in RTOL were NaN!"))
    args.IPAR[4] = FInt( length(args.LTOL) )

    OPT = OPT_BVPCLASS; 
    bvpclass = convert(FInt, getOption(opt, OPT, 1))
    @assert 0 ≤ bvpclass ≤ 3
    args.IPAR[1] = FInt( bvpclass == 0 ? 0 : 1)
    args.IPAR[10] = FInt( bvpclass == 0 ? 0 : bvpclass -1 )

    OPT = OPT_COLLOCATIONPTS
    opt_default = maximum([order_max + FInt(1), FInt(5)-order_max])
    args.IPAR[2] = convert(FInt, getOption(opt, OPT, opt_default))
    @assert order_max ≤ args.IPAR[2] ≤ 7

    OPT = OPT_MAXSUBINTERVALS
    max_subintervals = convert(FInt, getOption(opt, OPT, 50))
    @assert max_subintervals > 5

    OPT = OPT_SUBINTERVALS
    subintervals = getOption(opt, OPT, 5)
    if isscalar(subintervals) 
      args.IPAR[8] = FInt(0) # uniform(-like) initial grid
      args.IPAR[3] = FInt(subintervals)
      @assert 0 < args.IPAR[3] ≤ max_subintervals
    else
      subintervals = getVectorCheckLength(subintervals, Float64, 
        length(subintervals)) # always copy
      append!(subintervals, ζ)
      subintervals = unique(sort!(subintervals))
      subintervals = subintervals[ t_ab[1] .< subintervals .< t_ab[2] ]
      no_sub = length(subintervals)+2-1  # inner_nodes + a-node + b-node
      if no_sub > max_subintervals
        throw(ArgumentErrorODE(
          string("more than OPT_MAXSUBINTERVALS = ", max_subintervals,
            " subintervals in initial grid: ", subintervals)))
      end
      args.IPAR[3] = FInt(no_sub)
      OPT = OPT_FREEZEINTERVALS
      freezeintervals = getOption(opt, OPT, false)
      args.IPAR[8] = FInt( freezeintervals ? 2 : 1 )
      # division in subintervals needs to be saved in fspace (see below)
    end

    if guess == nothing
      args.IPAR[9] = FInt(0)
      guess = dummy_func
    else
      args.IPAR[9] = FInt(1)
    end
    # args.IPAR[9] ∈ {2,3,4} currently not supported by this interface

    OPT = OPT_DIAGNOSTICOUTPUT
    args.IPAR[7] = convert(FInt, getOption(opt, OPT, 1))
    @assert -1 ≤ args.IPAR[7] ≤ 1

    OPT = OPT_ADDGRIDPOINTS
    add_points = getOption(opt, OPT, [])
    add_points = getVectorCheckLength(add_points, Float64, length(add_points))
    @assert all(add_points .≥ t_ab[1]) && all(add_points .≤ t_ab[2])
    sort!(add_points)
    add_points = unique(add_points)
    args.FIXPNT = add_points
    args.IPAR[11] = FInt(length(add_points))
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid", :opt, e))
  end

  # FSPACE and ISPACE
  begin
    kd = FInt(args.IPAR[2]*n)
    kdm = FInt(kd + d)
    nsizei = FInt(3 + kdm)
    nrec = sum( ζ .≥ t_ab[2] )
    nsizef = 4 + 3*d + (5+kd)*kdm + (2*d-nrec)*2*d

    args.IPAR[5] = FInt(max_subintervals*nsizef)
    args.IPAR[6] = FInt(max_subintervals*nsizei) 

    args.FSPACE = zeros(Float64, args.IPAR[5])
    args.ISPACE = zeros(FInt, args.IPAR[6])

    if args.IPAR[8] ≠ 0 && args.IPAR[3] > 1
      # user given initial mesh
      args.FSPACE[2:args.IPAR[3]] = subintervals[:]
    end
  end

  args.IFLAG = zeros(FInt, 1)

  if colnew_global_cbi ≠ nothing
    throw(ArgumentErrorODE(string("The Fortran solver colnew does not ",
      "support 'pass-through' arguments. Hence this julia module does not ",
      "support concurrent/nested colnew calls. Sorry.")))
  end
  cbi = nothing
  try
    global colnew_global_cbi = cbi = ColnewInternalCallInfos(lio, l, n, d, 
      rhs, "unsafe_colnewrhs: ", 0,  Drhs, "unsafe_colnew_Drhs: ",0,
      bc, "unsafe_colnew_bc: ", 0,   Dbc, "unsafe_colnew_Dbc: ", 0,
      guess, "unsafe_colnew_guess")

    if l_solver
      println(lio, lprefix, "call Fortran-colnew $method_colnew with")
      dump(lio, args)
    end

    const fflag = FInt(0)

    ccall( method_colnew, Void,
      (Ptr{FInt}, Ptr{FInt},                      # NCOMP, M
       Ptr{Float64}, Ptr{Float64},                # ALEFT, ARIGHT
       Ptr{Float64}, Ptr{FInt},                   # ZETA, IPAR
       Ptr{FInt}, Ptr{Float64},                   # LTOL, TOL
       Ptr{Float64},                              # FIXPNT
       Ptr{FInt}, Ptr{Float64},                   # ISPACE, FSPACE
       Ptr{FInt},                                 # IFLAG
       Ptr{Void}, Ptr{Void},                      # rhs, Drhs
       Ptr{Void}, Ptr{Void},                      # bc, Dbc
       Ptr{Void},                                 # guess
      ),
      args.NCOMP, args.M,
      args.ALEFT, args.ARIGHT,
      args.ZETA, args.IPAR,
      args.LTOL, args.TOL,
      args.FIXPNT,
      args.ISPACE, args.FSPACE,
      args.IFLAG,
      unsafe_colnew_rhs_c(fflag), unsafe_colnew_Drhs_c(fflag),
      unsafe_colnew_bc_c(fflag), unsafe_colnew_Dbc_c(fflag),
      unsafe_colnew_guess_c(fflag),
    )
    if l_solver
      println(lio, lprefix, "Fortran-colnew $method_colnew returned")
      dump(lio, args)
    end
  finally
    global colnew_global_cbi = nothing
  end
  l_g && println(lio, lprefix, string("done IFALG=", args.IFLAG[1]))
  solobj = ColnewSolution(
    method_appsln, t_ab[1], t_ab[2], n, d, 
    args.ISPACE[1:(7+n)],
    args.FSPACE[1:args.ISPACE[7]])
  stats = Dict{AbstractString,Any}(
    "no_rhs_calls"  => cbi.rhs_count,
    "no_jac_calls"  => cbi.Drhs_count,
    "no_bc_calls"   => cbi.bc_count,
    "no_Dbc_calls"  => cbi.Dbc_count,
  )
  return (solobj, args.IFLAG[1], stats)
end

"""
        function evalSolution{FInt<:FortranInt}(sol::ColnewSolution{FInt},
          t::Real, z::Array{Float64})

  Evaluates an already obtained solution `sol` at time `t`.
  The values of the solution are saved in `z` which must be a 
  vector (of length d). 
  `t` must be in the interval [a,b] where the problem was solved.
  """
function evalSolution{FInt<:FortranInt}(sol::ColnewSolution{FInt},
  t::Real, z::AbstractArray{Float64})

  @assert ndims(z)==1
  @assert length(z)==sol.d
  @assert sol.t_a ≤ t ≤ sol.t_b
  ccall( sol.method_appsln, Void,
    (Ptr{Float64}, Ptr{Float64},                  # t, z
     Ptr{Float64}, Ptr{FInt},                     # FSPACE, ISPACE
    ),
    [Float64(t)], z,
    sol.FSPACE, sol.ISPACE)
  return nothing
end

"""
        function evalSolution{FInt<:FortranInt}(sol::ColnewSolution{FInt}, 
          t::Real)

  Evaluates an already obtained solution `sol` at time `t`.
  A newly allocated vector with the solution values is retured.
  `t` must be in the interval [a,b] where the problem was solved.
  """
function evalSolution{FInt<:FortranInt}(sol::ColnewSolution{FInt}, 
  t::Real)

  z = Vector{Float64}(sol.d)
  evalSolution(sol, t, z)
  return z
end

"""
        function evalSolution{FInt<:FortranInt}(sol::ColnewSolution{FInt}, 
          t::Vector)

  Evaluates an already obtained solution `sol` at time all
  times in the vector `t`.
  A newly allocated matrix of size `(length(t), d)` with the solution 
  values is retured.
  All values of `t` must be in the interval [a,b] where the problem was solved.
  """
function evalSolution{FInt<:FortranInt}(sol::ColnewSolution{FInt}, 
  t::Vector)

  vm = get_view_function()
  tno = length(t)

  Z = zeros(Float64, (tno, sol.d))
  for k=1:tno
    evalSolution(sol, t[k], vm(Z, k, :))
  end
  return Z
end


"""  
  ## Compile COLNEW

  The julia ODEInterface tries to compile and link the solvers
  automatically at the build-time of this module. The following
  calls need only be done, if one uses a different compiler and/or if
  one wants to change/add some compiler options.

  The Fortran source code can be found at:
  
       https://people.sc.fsu.edu/~jburkardt/f77_src/colnew/colnew.html
  
  See `help_colnew_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile BVPSOL with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o colnew.o colnew.f
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o colnew.so colnew.o
       gfortran -shared -fPIC -o colnew.dylib colnew.o
  
  ### Using `gfortran` and 64bit integers (Windows)
  
  Here is an example how to compile BVPSOL with `Float64` reals and
  `Int64` integers with `gfortran`:
  
       gfortran -c -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o colnew.o colnew.f
  
  In order to get create a shared library (from the object file above) use
  
       gfortran -shared -o colnew.dll colnew.o
  """
function help_colnew_compile()
  return Docs.doc(help_colnew_compile)
end

"""
  # Licence

  Colnew is licensed under the GNU LGPL license.

  See [Colnew Hompage](https://people.sc.fsu.edu/~jburkardt/f77_src/colnew/colnew.html).

  ```
                   GNU LESSER GENERAL PUBLIC LICENSE
                         Version 3, 29 June 2007
  
   Copyright (C) 2007 Free Software Foundation, Inc. <http://fsf.org/>
   Everyone is permitted to copy and distribute verbatim copies
   of this license document, but changing it is not allowed.
  
  
    This version of the GNU Lesser General Public License incorporates
  the terms and conditions of version 3 of the GNU General Public
  License, supplemented by the additional permissions listed below.
  
    0. Additional Definitions.
  
    As used herein, "this License" refers to version 3 of the GNU Lesser
  General Public License, and the "GNU GPL" refers to version 3 of the GNU
  General Public License.
  
    "The Library" refers to a covered work governed by this License,
  other than an Application or a Combined Work as defined below.
  
    An "Application" is any work that makes use of an interface provided
  by the Library, but which is not otherwise based on the Library.
  Defining a subclass of a class defined by the Library is deemed a mode
  of using an interface provided by the Library.
  
    A "Combined Work" is a work produced by combining or linking an
  Application with the Library.  The particular version of the Library
  with which the Combined Work was made is also called the "Linked
  Version".
  
    The "Minimal Corresponding Source" for a Combined Work means the
  Corresponding Source for the Combined Work, excluding any source code
  for portions of the Combined Work that, considered in isolation, are
  based on the Application, and not on the Linked Version.
  
    The "Corresponding Application Code" for a Combined Work means the
  object code and/or source code for the Application, including any data
  and utility programs needed for reproducing the Combined Work from the
  Application, but excluding the System Libraries of the Combined Work.
  
    1. Exception to Section 3 of the GNU GPL.
  
    You may convey a covered work under sections 3 and 4 of this License
  without being bound by section 3 of the GNU GPL.
  
    2. Conveying Modified Versions.
  
    If you modify a copy of the Library, and, in your modifications, a
  facility refers to a function or data to be supplied by an Application
  that uses the facility (other than as an argument passed when the
  facility is invoked), then you may convey a copy of the modified
  version:
  
     a) under this License, provided that you make a good faith effort to
     ensure that, in the event an Application does not supply the
     function or data, the facility still operates, and performs
     whatever part of its purpose remains meaningful, or
  
     b) under the GNU GPL, with none of the additional permissions of
     this License applicable to that copy.
  
    3. Object Code Incorporating Material from Library Header Files.
  
    The object code form of an Application may incorporate material from
  a header file that is part of the Library.  You may convey such object
  code under terms of your choice, provided that, if the incorporated
  material is not limited to numerical parameters, data structure
  layouts and accessors, or small macros, inline functions and templates
  (ten or fewer lines in length), you do both of the following:
  
     a) Give prominent notice with each copy of the object code that the
     Library is used in it and that the Library and its use are
     covered by this License.
  
     b) Accompany the object code with a copy of the GNU GPL and this license
     document.
  
    4. Combined Works.
  
    You may convey a Combined Work under terms of your choice that,
  taken together, effectively do not restrict modification of the
  portions of the Library contained in the Combined Work and reverse
  engineering for debugging such modifications, if you also do each of
  the following:
  
     a) Give prominent notice with each copy of the Combined Work that
     the Library is used in it and that the Library and its use are
     covered by this License.
  
     b) Accompany the Combined Work with a copy of the GNU GPL and this license
     document.
  
     c) For a Combined Work that displays copyright notices during
     execution, include the copyright notice for the Library among
     these notices, as well as a reference directing the user to the
     copies of the GNU GPL and this license document.
  
     d) Do one of the following:
  
         0) Convey the Minimal Corresponding Source under the terms of this
         License, and the Corresponding Application Code in a form
         suitable for, and under terms that permit, the user to
         recombine or relink the Application with a modified version of
         the Linked Version to produce a modified Combined Work, in the
         manner specified by section 6 of the GNU GPL for conveying
         Corresponding Source.
  
         1) Use a suitable shared library mechanism for linking with the
         Library.  A suitable mechanism is one that (a) uses at run time
         a copy of the Library already present on the user's computer
         system, and (b) will operate properly with a modified version
         of the Library that is interface-compatible with the Linked
         Version.
  
     e) Provide Installation Information, but only if you would otherwise
     be required to provide such information under section 6 of the
     GNU GPL, and only to the extent that such information is
     necessary to install and execute a modified version of the
     Combined Work produced by recombining or relinking the
     Application with a modified version of the Linked Version. (If
     you use option 4d0, the Installation Information must accompany
     the Minimal Corresponding Source and Corresponding Application
     Code. If you use option 4d1, you must provide the Installation
     Information in the manner specified by section 6 of the GNU GPL
     for conveying Corresponding Source.)
  
    5. Combined Libraries.
  
    You may place library facilities that are a work based on the
  Library side by side in a single library together with other library
  facilities that are not Applications and are not covered by this
  License, and convey such a combined library under terms of your
  choice, if you do both of the following:
  
     a) Accompany the combined library with a copy of the same work based
     on the Library, uncombined with any other library facilities,
     conveyed under the terms of this License.
  
     b) Give prominent notice with the combined library that part of it
     is a work based on the Library, and explaining where to find the
     accompanying uncombined form of the same work.
  
    6. Revised Versions of the GNU Lesser General Public License.
  
    The Free Software Foundation may publish revised and/or new versions
  of the GNU Lesser General Public License from time to time. Such new
  versions will be similar in spirit to the present version, but may
  differ in detail to address new problems or concerns.
  
    Each version is given a distinguishing version number. If the
  Library as you received it specifies that a certain numbered version
  of the GNU Lesser General Public License "or any later version"
  applies to it, you have the option of following the terms and
  conditions either of that published version or of any later version
  published by the Free Software Foundation. If the Library as you
  received it does not specify a version number of the GNU Lesser
  General Public License, you may choose any version of the GNU Lesser
  General Public License ever published by the Free Software Foundation.
  
    If the Library as you received it specifies that a proxy can decide
  whether future versions of the GNU Lesser General Public License shall
  apply, that proxy's public statement of acceptance of any version is
  permanent authorization for you to choose that version for the
  Library.
  ```
  
  """
function help_colnew_license()
  return Docs.doc(help_colnew_license)
end

# Add informations about solvers in global solverInfo-array.
push!(solverInfo,
  SolverInfo("colnew",
    "Multi-Point boundary value solver for mixed order systems with collocation",
    tuple( :OPT_BVPCLASS, :OPT_RTOL, :OPT_COLLOCATIONPTS, 
           :OPT_SUBINTERVALS, :OPT_FREEZEINTERVALS, 
           :OPT_ADDGRIDPOINTS, :OPT_MAXSUBINTERVALS, 
           :OPT_DIAGNOSTICOUTPUT,
          ),
    tuple(
      SolverVariant("colnew_i64",
        "Colnew with 64bit integers",
        DL_COLNEW,
        tuple("colnew", "appsln")),
      SolverVariant("colnew_i32",
        "Colnew with 32bit integers",
        DL_COLNEW_I32,
        tuple("colnew", "appsln")),
    ),
    help_colnew_compile,
    help_colnew_license,
  )
)


# vim:syn=julia:cc=79:fdm=indent:

