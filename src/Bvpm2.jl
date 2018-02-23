# Functions for Boundary-Value Solver: BVP_M-2

import Base: show

"""Name for Loading bvpm2 solver (64bit integers)."""
const DL_BVPM2                = "bvp_m_proxy"

# There is *no* 32-bit variant

"""macro for import bvpm2 solver."""
macro import_bvpm2()
  :(
    using ODEInterface: Bvpm2, bvpm2_copy, bvpm2_terminate, bvpm2_destroy,
                        bvpm2_init, bvpm2_show_details,
                        bvpm2_get_details,
                        bvpm2_get_x, bvpm2_solve, bvpm2_extend,
                        evalSolution,
                        bvpm2_get_params
  )
end

"""macro for import bvpm2 dynamic lib name."""
macro import_DLbvpm2()
  :(
    using ODEInterface: DL_BVPM2
  )
end

"""macro for importing Bvpm2 help."""
macro import_bvpm2_help()
  :(
    using ODEInterface: help_bvpm2_compile, help_bvpm2_license,
                        help_bvpm2_proxy
  )
end

"""
  BVP_M-2 uses global variables during the solution process.
  Hence we can only support one bvpm2 solve-call at a time.
  """
bvpm2_global_cbi = nothing

"""
  structure to save guess-function for use in Fortran callbacks.

  How to make a julia function (which typically cannot be put in `cfunction`,
  because it may be a closure, etc.) callable from Fortran2003 code that
  uses no ISO_C_BINDING and has no pass-through arguments?

  Use a Fortran2003 proxy with ISO_C_BINDING. Here is the calling stack:

  ```
  bvpm2_guess( ...  guess_fcn ... )
  ─────────────────────────────────
  │ [ cbi created with guess_fcn ]
  │ call init_guess3_c( ... guess_fcn_ptr=unsafe_bvpm2_guess_cb_c,
  │ │    ━━━━━━━━━━━━━      guess_pthrough=cbi ... )
  │ │ call bvp_init( ... guess = guess_fcn_proxy ... )
  │ │ │    ════════
  │ │ │ call guess_fcn_proxy(x_point, guess_vector)
  │ │ │ │    ━━━━━━━━━━━━━━━
  │ │ │ │ [nested guess_fcn_proxy in scope of init_guess3_c, so
  │ │ │ │  the guess_pthrough info is available.]
  │ │ │ │ [convert Fortran Arrays to C-Pointer-Array]
  │ │ │ │ call guess_fcn_ptr/unsafe_bvpm2_guess_cb_c(x_point, ...,
  │ │ │ │ │                  ───────────────────────   guess_vector 
  │ │ │ │ │                                            guess_pthrough)
  │ │ │ │ │   [ use guess_pthrough to recover cbi ]
  │ │ │ │ │   call cbi.guess_fcn(x, guess_vector)
  │ │ │ │ │   │    ─────────────
  └ └ └ └ └   └

  Legend:
     ─────────  julia code 
     ━━━━━━━━━  Fortran2003 code in BVP_M_Proxy with ISO_C_BINDING
     ═════════  Fortran2003 code in BVP_M-2 (without ISO_C_BINDING)
  ```

  For more details of this concept/idea call `help_bvpm2_proxy`.
  """
mutable struct Bvpm2_guess_cbi{GUESS_FCN} <: ODEinternalCallInfos
  no_odes       :: Int64            # number of ODEs, used for assertion
  # GUESS:
  guess         :: GUESS_FCN        # Julia-function to call for guess
end

mutable struct Bvpm2_solve_cbi{RHS_FCN, DRHS_FCN,
              BC_FCN, DBC_FCN} <: ODEinternalCallInfos
  logio        :: IO              # where to log
  loglevel     :: UInt64          # log level
  no_odes      :: Int64           # number of ODEs, used for assertion
  no_par       :: Int64           # number of parameters, used for assertion
  no_left_bc   :: Int64           # number of boundary conditions on left side
                                  # used for assertion
  # RHS:
  rhs          :: RHS_FCN         # Julia-function to call for ODEs right-hand
                                  # side
  rhs_lprefix  :: AbstractString  # saved log-prefix for rhs
  rhs_count    :: Int             # count: number of calls to rhs
  # DRHS:
  Drhs         :: DRHS_FCN        # Jacobian(s) of rhs
  Drhs_lprefix :: AbstractString  # saved log-prefix for Drhs
  Drhs_count   :: Int             # count: number of calls to Drhs
  # BC:
  bc_fcn       :: BC_FCN          # Julia-function to call for evaluating
                                  # boundary conditions
  bc_lprefix   :: AbstractString  # saved log-prefix for bc
  bc_count     :: Int             # count: number of calls to bc
  # DBC:
  Dbc          :: DBC_FCN         # Jacobian(s) of side conditions
  Dbc_lprefix  :: AbstractString  # saved log-prefix for Dbc
  Dbc_count    :: Int             # count: number of calls to Dbc
end

"""
  # Bvpm2 object for solving boundary value problems

  This is the Julia part of the BVP_M-2 (Fortran-)solution object. 
  For (nearly) all the operations the corresponding Fortran-Proxy 
  methods are called (call `help_bvpm2_proxy()` to get internal details).

  ## Boundary value problem (BVP)

  BVPs of the following form are considered:

                    1
        y'(x) =  ─────── Sy + f(x, y, p)         for a ≤ x ≤ b   [ODEs]
                  x - a


        ga(y(a), p) = 0,     gb(y(b), p) = 0                     [BCs]

  * y(x) ∈ ℝᵈ and `d` is also called `no_odes` (the number of ordinary
    differential equations). 
  * S ∈ ℝᵈˣᵈ is an optional constant matrix (also 
    called the singularity term) because the whole term S⋅y/(x-a) has a 
    singularity at x=a. If S is not given, then the ODEs are reduced to
    y'(x) = f(x, y, p).
  * p ∈ ℝᵐ (with 0≤m) are unknown parameters of the problem. `m` is also
    called `no_par` (the number of parameters).
  * f(x, y, p) ∈ ℝᵈ is also called the right-hand side (of the ODEs).
  * ga(ya, p) ∈ ℝˡ describes the left boundary conditions. `l` is
    also called `no_left_bc` (the number of the BCs at x=a).
  * ga(yb, p) ∈ ℝⁿ describes the right boundary conditions. It is

          n = d + m - l 
          n = no_odes + no_par - no_left_bc

  ## Initial guess and solutions

  A Bvpm2 object can be used to represent either an initial guess (for a 
  BVP like above) or a solution. It is possible to use a solution of a 
  (different) BVP as initial guess to another BVP.

  Such a Bvpm2 object can be in one of the following states:
  
  * `state==0`: object created (and connected to Fortran-object), but 
    not initialized, i.e. it does neither represent a guess nor an solution.
  * `state==1`: object created, and initialized with an (initial) guess, i.e.
    the object represents a guess.
  * `state==2`: object created and a solution was calculated successfully and
    saved in the object, i.e. the object represents a solution.
  * `state==-1`: object is not connected to a Fortran-Proxy. Either
    `bvpm2_destroy` was called or at creation time, the connection to the 
     Fortran-Proxy couldn't be established, i.e. the object is unusable and
     all associated memory was deallocated.

  The following table shows possible actions and the state-transitions
  initiated by the actions.

  ```
  ╔═══════════════════╤═══════════════════════════╤════════════╤════════════╗
  ║ Action/Function   │ Description               │state before│state after ║
  ╠═══════════════════╪═══════════════════════════╪════════════╪════════════╣
  ║ Bvpm2()           │ create object             │    ---     │     0      ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_init        │ initialize object with    │     0      │     1      ║
  ║                   │ initial guess, etc.       │            │            ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_show_details│ show some details of      │ 0, 1, or 2 │ not changed║
  ║                   │ (Fortran-)BVP_M-2 sol     │            │            ║
  ║                   │ object                    │            │            ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_get_details │ get dict with some details│ -1, 0, 1,  │ not changed║
  ║                   │ of the Bvpm2 object:      │   or 2     │            ║
  ║                   │ e.g. state, number of pts │            │            ║
  ║                   │ in current grid ...       │            │            ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_get_x       │ return current grid of    │ 1, or 2    │ not changed║
  ║                   │ the object.               │            │            ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_get_params  │ return current values of  │ 1, or 2    │ not changed║
  ║                   │ estimated/calculated      │            │            ║
  ║                   │ parameters                │            │            ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_terminate   │ throw away all data and   │ 0, 1, or 2 │      0     ║
  ║                   │ information. Bring to     │            │            ║
  ║                   │ state 0.                  │            │            ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_destroy     │ deallocate all (Fortran-) │ -1, 0, 1,  │     -1     ║
  ║                   │ resources for this object.│   or 2     │            ║
  ║                   │                           │            │            ║
  ╚═══════════════════╧═══════════════════════════╧════════════╧════════════╝ 
  ```

  There are functions that take an Bvpm2-object `obj_in` as input, 
  perhaps change `obj_in` and create an additonal `obj_out`.

  The following table shows possible actions, the change of the state
  of `obj_in` and which `obj_out` object is created:

  ```
  ╔═══════════════════╤═══════════════════════════╤════════════╤════════════╗
  ║ Action/Function   │ Description               │obj_in state│  state of  ║
  ║                   │                           │ from -> to │   obj_out  ║
  ╠═══════════════════╪═══════════════════════════╪════════════╪════════════╣
  ║ bvpm2_solve       │ take obj_in as guess. Do  │ not changed│ 0: no      ║
  ║                   │ not change obj_in.        │            │    success ║
  ║                   │ Produces obj_out object   │            │            ║
  ║                   │ representing the solution.│            │ 2: success ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_copy        │ create deep copy of       │ not changed│ same as    ║
  ║                   │ obj_in                    │            │ obj_in     ║
  ╟───────────────────┼───────────────────────────┼────────────┼────────────╢
  ║ bvpm2_extend      │ extend solution to new    │  2 -> 0    │    1       ║
  ║                   │ interval as new  guess.   │            │            ║
  ║                   │ obj_in will be terminated │            │            ║
  ║                   │ and will be in state 0.   │            │            ║
  ║                   │ Call bvpm2_copy before, if│            │            ║
  ║                   │ you need the solution     │            │            ║
  ║                   │ later on.                 │            │            ║
  ╚═══════════════════╧═══════════════════════════╧════════════╧════════════╝ 
  ```

  """
mutable struct Bvpm2   <: AbstractODESolution{Int64}
  handle                 :: Ptr{Cvoid}
  method_create          :: Ptr{Cvoid}
  method_copy            :: Ptr{Cvoid}
  method_terminate       :: Ptr{Cvoid}
  method_destroy         :: Ptr{Cvoid}
  method_show            :: Ptr{Cvoid}
  method_get_details     :: Ptr{Cvoid}
  method_get_x           :: Ptr{Cvoid}
  method_init_guess1     :: Ptr{Cvoid}
  method_init_guess2     :: Ptr{Cvoid}
  method_init_guess3     :: Ptr{Cvoid}
  method_solve           :: Ptr{Cvoid}
  method_eval_s          :: Ptr{Cvoid}
  method_eval_v          :: Ptr{Cvoid}
  method_get_params      :: Ptr{Cvoid}
  method_extend_s        :: Ptr{Cvoid}
  method_extend_e        :: Ptr{Cvoid}
end

"""
       function Bvpm2(handle::Ptr{Cvoid})

  create bvpm2 object from given handle.
  """
function Bvpm2(handle::Ptr{Cvoid})
  obj = Bvpm2(C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
              C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL,
              C_NULL)
  get_proxy_methods(obj)
  obj.handle = handle
  return obj
end

"""
       function bvpm2_create_handle(obj::Bvpm2)
  
  create Fortran Proxy.
  """
function bvpm2_create_handle(obj::Bvpm2)
  obj.handle = ccall(obj.method_create, Ptr{Cvoid}, () )
  @assert obj.handle ≠ C_NULL
  return nothing
end

"""
       function Bvpm2()

  creates bvpm2 object.
  """
function Bvpm2()
  obj = Bvpm2(C_NULL)
  bvpm2_create_handle(obj)
  return obj
end

"""
  get all (proxy-)methods for bvpm2 object.
  """
function get_proxy_methods(obj::Bvpm2)
  (obj.method_create,
   obj.method_copy,
   obj.method_terminate,
   obj.method_destroy,
   obj.method_show,
   obj.method_get_details,
   obj.method_get_x,
   obj.method_init_guess1,
   obj.method_init_guess2,
   obj.method_init_guess3,
   obj.method_solve, 
   obj.method_eval_s,
   obj.method_eval_v,
   obj.method_get_params,
   obj.method_extend_s,
   obj.method_extend_e,
   ) = getAllMethodPtrs(DL_BVPM2)
  return nothing
end

bvpm2_is_handle_valid(obj::Bvpm2) = obj.handle ≠ C_NULL

"""
  tests if Bvpm2-object is "connected" to a Fortran-proxy handle.
  """
function bvpm2_check_handle(obj::Bvpm2)
  if !bvpm2_is_handle_valid(obj)
    throw(StateErrorODE(string(
      "Bvpm2 object is not connected to a Fortran-proxy. ",
      "Either the connection failed or bvpm2_destroy was called before ",
      "with this Bvpm2 object.")))
  end
  return nothing
end

"""
  tests if the Bvpm2-object is in one of the expected states.

  Throws an Error of the state is not expected.

  Returns details dict.
  """
function bvpm2_check_state(obj::Bvpm2, expected_states)
  details = bvpm2_get_details(obj)
  state = details["state"]

  if state ∉ expected_states
    throw(StateErrorODE(string(
      "Bvpm2 object is in state ", state,", but I need one of ",
      "the states ", expected_states)))
  end

  return details
end

"""
       function bvpm2_terminate(obj::Bvpm2)

  Terminates Bvpm2-object. Put in state as if after creation.
  """
function bvpm2_terminate(obj::Bvpm2)
  bvpm2_check_handle(obj)
  ccall(obj.method_terminate, Cvoid, (Ptr{Cvoid},), obj.handle)
  return nothing
end

"""
       function bvpm2_destroy(obj::Bvpm2)

  destroys Bvpm2-object. Especially free all (Fortran-)allocated
  memory and disconnect from Fortran-Proxy(-handle).
  """
function bvpm2_destroy(obj::Bvpm2)
  if bvpm2_is_handle_valid(obj)
    ccall(obj.method_destroy, Cvoid, (Ptr{Cvoid},), obj.handle)
    obj.handle = C_NULL
  end
  return obj
end

"""
       function bvpm2_copy(obj_in::Bvpm2) -> obj_out

  creates a deep copy `obj_out` of the Bvpm2-object `obj_in`.
  """
function bvpm2_copy(obj::Bvpm2)
  bvpm2_check_state(obj, (0,1,2,))
  new_handle = ccall(obj.method_copy, Ptr{Cvoid}, (Ptr{Cvoid},),
                obj.handle)
  return Bvpm2(new_handle)
end

"""
       function bvpm2_show_details(obj::Bvpm2)

  Debug: call show_details Fortran subroutine.
  """
function bvpm2_show_details(obj::Bvpm2)
  bvpm2_check_handle(obj)
  ccall(obj.method_show, Cvoid, (Ptr{Cvoid},), obj.handle)
  return nothing
end

"""
       function bvpm2_get_details(obj::Bvpm2)

  returns dict with informations about an Bvpm2 object.
  """
function bvpm2_get_details(obj::Bvpm2)
  if !bvpm2_is_handle_valid(obj)
    return Dict("state" => -1 )
  end
  int_slots = zeros(Int64, 16)

  ccall(obj.method_get_details, Cvoid,
    (Ptr{Cvoid},                     # handle
     Int64, Ref{Int64},              # int_slots_len, int_slots
    ),
    obj.handle,
    length(int_slots), int_slots
  )
  return Dict(
    "handle"               => obj.handle,
    "state"                => int_slots[1],
    "no_odes"              => int_slots[2],
    "no_par"               => int_slots[3],
    "no_left_bc"           => int_slots[4],
    "no_pts"               => int_slots[5],
    "info"                 => int_slots[6],
    "max_num_subintervals" => int_slots[7],
    "len_iwork"            => int_slots[8],
    "len_work"             => int_slots[9]
  )
end

function show(io::IO, obj::Bvpm2)
  println(io, typeof(obj))
  details = bvpm2_get_details(obj)
  maxLen = 2+max( 0, map(length, collect(keys(details)))... )
  for key in sort(collect(keys(details)))
    print(io, lpad(key, maxLen), ": ")
    show(io, details[key]); println(io)
  end
end

"""
  returns current vector with x-grid of Bvpm2 object.
  """
function bvpm2_get_x(obj::Bvpm2)
  details = bvpm2_check_state(obj, (1, 2))
  x = Vector{Float64}(uninitialized, details["no_pts"])
  error = ccall(obj.method_get_x, Int64,
    (Ptr{Cvoid}, Int64, Ref{Float64},),   # handle, x_len, x
    obj.handle, length(x), x)
  error ≠ 0 && throw(InternalErrorODE(string(
    "Sorry. Fortran-Proxy returned error ",error)))
  return x
end

"""
       function bvpm2_get_params(obj::Bvpm2)

  returns current vector with parameters of Bvpm2 object.
  """
function bvpm2_get_params(obj::Bvpm2)
  details = bvpm2_check_state(obj, (1, 2))
  p = Vector{Float64}(uninitialized, details["no_par"])
  if length(p) > 0 
    error = ccall(obj.method_get_params, Int64,
      (Ptr{Cvoid}, Int64, Ref{Float64},),   # handle, p_len, p
      obj.handle, length(p), p)
    error ≠ 0 && throw(InternalErrorODE(string(
      "Sorry. Fortran-Proxy returned error ",error)))
  end
  return p
end

"""
  same as bvpm2_get_x.
  """
function getSolutionGrid(sol::Bvpm2)
  return bvpm2_get_x(sol)
end

"""
       function bvpm2_init_tests(obj::Bvpm2, no_odes, no_left_bc, 
           x_grid::Vector, parameters::Vector, max_num_subintervals)

  init tests for common parameters.
  """
function bvpm2_init_tests(obj::Bvpm2, no_odes, no_left_bc, 
    x_grid::Vector, parameters::Vector, max_num_subintervals)

  details = bvpm2_check_state(obj, (0,))

  try
    no_odes = convert(Int64, no_odes)
    @assert no_odes > 0
  catch e
    throw(ArgumentErrorODE("no_odes must be positive integer",
      :no_odes, e))
  end

  try
    no_left_bc = convert(Int64, no_left_bc)
    @assert no_left_bc ≥ 0
  catch e
    throw(ArgumentErrorODE("no_left_bc must be non-negative integer", 
      :no_left_bc, e))
  end

  len_x_grid = length(x_grid)
  if len_x_grid <= 2
    throw(ArgumentErrorODE(string(
      "x_grid must have length > 2. You must specify an initial grid."), 
      :x_grid))
  end
  try
    x_grid = getVectorCheckLength(x_grid, Float64, len_x_grid)
    sort!(x_grid)
  catch e
    throw(ArgumentErrorODE("cannot convert x_grid to Float64 vector",
      :x_grid, e))
  end
  
  try
    parameters = getVectorCheckLength(parameters, Float64,
      length(parameters))
  catch e
    throw(ArgumentErrorODE("cannot convert parameters to a Float64 vector",
      :parameters, e))
  end
  
  try
    max_num_subintervals = convert(Int64, max_num_subintervals)
    @assert max_num_subintervals > 0
    @assert max_num_subintervals ≥ length(x_grid)-1
  catch e
    throw(ArgumentErrorODE("max_num_subintervals must be positive integer", 
      :max_num_subintervals, e))
  end

  return details, no_odes, no_left_bc, x_grid, parameters, max_num_subintervals
end

"""
       function bvpm2_init(obj::Bvpm2,
         no_odes, no_left_bc, x_grid::Vector, constant_guess::Vector, 
         parameters::Vector=[], max_num_subintervals=3000)

  initialize Bvpm2 object with a constant intial guess.
  """
function bvpm2_init(obj::Bvpm2,
  no_odes, no_left_bc, x_grid::Vector, constant_guess::Vector, 
  parameters::Vector=[], max_num_subintervals=3000)
  
  details, no_odes, no_left_bc, x_grid, parameters, max_num_subintervals =
    bvpm2_init_tests(obj, no_odes, no_left_bc, x_grid, parameters, 
                     max_num_subintervals)

  try
    constant_guess = getVectorCheckLength(constant_guess, 
      Float64, no_odes)
  catch e
    throw(ArgumentErrorODE(string("cannot convert constant_guess to a ",
      "Float64 vector of length ",no_odes), :constant_guess, e))
  end
  
  ccall(obj.method_init_guess1, Cvoid,
    (Ptr{Cvoid},                           # handle
     Int64, Int64,                         # no_odes, no_left_bc
     Int64, Ref{Float64},                  # x_len, x
     Int64, Ref{Float64},                  # y_len, y
     Int64, Ref{Float64},                  # p_len, p
     Int64,                                # max_num_subintervals
    ),
    obj.handle,
    no_odes, no_left_bc,
    length(x_grid), x_grid,
    length(constant_guess), constant_guess,
    length(parameters), parameters,
    max_num_subintervals)
  
  return obj
end

"""
       function bvpm2_init(obj::Bvpm2,
         no_odes, no_left_bc, x_grid::Vector, guess::Matrix, 
         parameters::Vector=[], max_num_subintervals=3000)

  initialize Bvpm2 object with a guess for every state at
  every node in x_grid.
  """
function bvpm2_init(obj::Bvpm2,
  no_odes, no_left_bc, x_grid::Vector, guess::Matrix, 
  parameters::Vector=[], max_num_subintervals=3000)
  
  details, no_odes, no_left_bc, x_grid, parameters, max_num_subintervals =
    bvpm2_init_tests(obj, no_odes, no_left_bc, x_grid, parameters, 
                     max_num_subintervals)
  try
    m,n = size(guess)
    if m ≠ no_odes || n ≠ length(x_grid)
      throw(ArgumentErrorODE(string(
        "guess must be a (no_odes, length(x_grid))=(", no_odes,", ",
        length(x_grid)," matrix. But I found a (",m, ", ", n," matrix.")))
    end
    converted_guess = zeros(no_odes, length(x_grid))
    converted_guess[:] = getVectorCheckLength(guess[:],
      Float64, no_odes*length(x_grid), false)
    guess = converted_guess
  catch e
    throw(ArgumentErrorODE(string("cannot convert guess to a ",
      "Float64 (", no_odes, ", ", length(x_grid), ") matrix"), :guess, e))
  end

  ccall(obj.method_init_guess2, Cvoid,
    (Ptr{Cvoid},                           # handle
     Int64, Int64,                         # no_odes, no_left_bc
     Int64, Ref{Float64},                  # x_len, x
     Int64, Int64, Ref{Float64},           # y_dim1, y_dim2, y
     Int64, Ref{Float64},                  # p_len, p
     Int64,                                # max_num_subintervals
    ),
    obj.handle,
    no_odes, no_left_bc,
    length(x_grid), x_grid,
    size(guess, 1), size(guess, 2), guess,
    length(parameters), parameters,
    max_num_subintervals)

  return obj
end

"""
  This function calls `guess_fcn` saved in Bvpm2_guess_cbi.
  """
function bvpm2_guess(x, y, cbi::CI) where CI
  # lprefix = cbi.guess_lprefix
  # (lio,l)=(cbi.logio,cbi.loglevel)
  # l_guess = l & LOG_GUESS > 0

  # l_guess && println(lio, lprefix, "called with x=", x)
  cbi.guess(x, y)
  # l_guess && println(lio, lprefix, "guess result=", y)
  return nothing
end

"""
  This is the guess function given as callback to bvpm2.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_guess_cb(x::Float64, y_len::Int64, y_::Ptr{Float64}, 
     cbi::CI) where CI<:Bvpm2_guess_cbi

  @assert y_len == cbi.no_odes
  y = unsafe_wrap(Array, y_, (y_len,), own=false)
  bvpm2_guess(x, y, cbi)
  return nothing
end

function unsafe_bvpm2_guess_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_guess_cb, Cvoid,
    Tuple{Float64, Int64, Ptr{Float64}, Ref{CI}})
end

"""
       function bvpm2_init(obj, no_odes, no_left_bc, x_grid, 
                           guess<:Function, parameters, 
                           max_num_subintervals=3000)

  The guess function must have the form

       function guess(x,y)

  where inside the function the guess for position x has to be
  stored in y.

  initialize Bvpm2 object where the function `guess` is used
  to get the guesses for the state at different `x` values.
  """
function bvpm2_init(obj::Bvpm2,
  no_odes, no_left_bc, x_grid::Vector, guess::GUESS_FCN, 
  parameters::Vector=[], max_num_subintervals=3000) where GUESS_FCN<:Function

  details, no_odes, no_left_bc, x_grid, parameters, max_num_subintervals =
    bvpm2_init_tests(obj, no_odes, no_left_bc, x_grid, parameters, 
                     max_num_subintervals)

  # Create callback-info
  cbi = Bvpm2_guess_cbi(no_odes, guess)

  ccall(obj.method_init_guess3, Cvoid,
    (Ptr{Cvoid},                           # handle
     Int64, Int64,                         # no_odes, no_left_bc
     Int64, Ref{Float64},                  # x_len, x
     Ptr{Cvoid}, Ref{Bvpm2_guess_cbi},     # guess_fcn, guess_fcn_pthrough
     Int64, Ref{Float64},                  # p_len, p
     Int64,                                # max_num_subintervals
    ),
    obj.handle,
    no_odes, no_left_bc,
    length(x_grid), x_grid,
    unsafe_bvpm2_guess_cb_c(cbi), cbi,
    length(parameters), parameters,
    max_num_subintervals)

  return obj
end

"""
  This function calls `rhs` saved in Bvpm2_solve_cbi.
  """
function bvpm2_rhs(x, y, f, cbi::CI, p=nothing) where CI
  lprefix = cbi.rhs_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_rhs = l & LOG_RHS > 0

  cbi.rhs_count += 1

  l_rhs && println(lio, lprefix, "called with x=", x, " y=", y, " p=", p)
  if p == nothing
    cbi.rhs(x, y, f)
  else
    cbi.rhs(x, y, p, f)
  end
  l_rhs && println(lio, lprefix, "rhs result=", f)
  return nothing
end

"""
  This is the right-hand side given as callback to bvpm2 in the
  case where the problem has no unkown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_rhs_cb(
    x::Float64, y_len::Int64, y_::Ptr{Float64}, 
    f_len::Int64, f_::Ptr{Float64}, cbi::CI) where CI<:Bvpm2_solve_cbi

  @assert y_len == f_len == cbi.no_odes

  y = unsafe_wrap(Array, y_, (y_len,), own=false)
  f = unsafe_wrap(Array, f_, (f_len,), own=false)
  bvpm2_rhs(x, y, f, cbi)
  return nothing
end

function unsafe_bvpm2_rhs_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_rhs_cb, Cvoid,
    Tuple{Float64, Int64, Ptr{Float64}, Int64, Ptr{Float64}, Ref{CI}})
end

"""
  This is the right-hand side given as callback to bvpm2 in the
  case where the problem has unknown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_rhspar_cb(
    x::Float64, y_len::Int64, y_::Ptr{Float64}, 
    p_len::Int64, p_::Ptr{Float64},
    f_len::Int64, f_::Ptr{Float64}, cbi::CI) where CI<:Bvpm2_solve_cbi

  @assert y_len == f_len == cbi.no_odes
  @assert p_len == cbi.no_par

  y = unsafe_wrap(Array, y_, (y_len,), own=false)
  p = unsafe_wrap(Array, p_, (p_len,), own=false)
  f = unsafe_wrap(Array, f_, (f_len,), own=false)
  bvpm2_rhs(x, y, f, cbi, p)
  return nothing
end

function unsafe_bvpm2_rhspar_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_rhspar_cb, Cvoid,
    Tuple{Float64, Int64, Ptr{Float64}, Int64, Ptr{Float64},
     Int64, Ptr{Float64}, Ref{CI}})
end

"""
  This function calls `Drhs` saved in Bvpm2_solve_cbi.
  """
function bvpm2_Drhs(x, y, dfdy, cbi::CI, p=nothing, dfdp=nothing) where CI
  lprefix = cbi.Drhs_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_Drhs = l & LOG_JAC > 0

  cbi.Drhs_count += 1
  l_Drhs && println(lio, lprefix, "called with x=", x, " y=", y, " p=", p)
  if p == nothing
    cbi.Drhs(x, y, dfdy)
  else
    cbi.Drhs(x, y, p, dfdy, dfdp)
  end
  l_Drhs && println(lio, lprefix, "Drhs result: dfdy=", dfdy, " dfdp=", dfdp)
  return nothing
end

"""
  This is the derivative of the right-hand side given as callback to
  bvpm2 in the case where the problem has no unknown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_Drhs_cb(
  x::Float64, y_len::Int64, y_::Ptr{Float64},
  dfdy_dim1::Int64, dfdy_dim2::Int64, dfdy_::Ptr{Float64}, 
  cbi::CI) where CI<:Bvpm2_solve_cbi

  @assert y_len == cbi.no_odes
  @assert dfdy_dim1 == dfdy_dim2 == cbi.no_odes
  
  y = unsafe_wrap(Array, y_, (y_len,), own=false)
  dfdy = unsafe_wrap(Array, dfdy_, (dfdy_dim1, dfdy_dim2,), own=false)
  bvpm2_Drhs(x, y, dfdy, cbi)
  return nothing
end

function unsafe_bvpm2_Drhs_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_Drhs_cb, Cvoid,
    Tuple{Float64, Int64, Ptr{Float64}, Int64, Int64, Ptr{Float64}, Ref{CI}})
end

"""
  This is the derivative of the right-hand side given as callback to
  bvpm2 in the case where the problem has unknown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_Drhspar_cb(
  x::Float64, y_len::Int64, y_::Ptr{Float64},
  p_len::Int64, p_::Ptr{Float64},
  dfdy_dim1::Int64, dfdy_dim2::Int64, dfdy_::Ptr{Float64}, 
  dfdp_dim1::Int64, dfdp_dim2::Int64, dfdp_::Ptr{Float64},
  cbi::CI) where CI<:Bvpm2_solve_cbi

  @assert y_len == cbi.no_odes
  @assert dfdy_dim1 == dfdy_dim2 == cbi.no_odes
  @assert dfdp_dim1 == cbi.no_odes && dfdp_dim2 == cbi.no_par
  
  y = unsafe_wrap(Array, y_, (y_len,), own=false)
  p = unsafe_wrap(Array, p_, (p_len,), own=false)
  dfdy = unsafe_wrap(Array, dfdy_, (dfdy_dim1, dfdy_dim2,), own=false)
  dfdp = unsafe_wrap(Array, dfdp_, (dfdp_dim1, dfdp_dim2,), own=false)
  bvpm2_Drhs(x, y, dfdy, cbi, p, dfdp)
  return nothing
end

function unsafe_bvpm2_Drhspar_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_Drhspar_cb, Cvoid,
    Tuple{Float64, Int64, Ptr{Float64}, Int64, Ptr{Float64},
    Int64, Int64, Ptr{Float64}, Int64, Int64, Ptr{Float64}, Ref{CI}})
end

"""
  This function calls `bc` saved in Bvpm2_solve_cbi.
  """
function bvpm2_bc(ya, yb, bca, bcb, cbi::CI, p=nothing) where CI
  lprefix = cbi.bc_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_bc = l & LOG_BC > 0

  cbi.bc_count += 1
  l_bc && println(lio, lprefix, "called with ya=", ya, " yb=", yb, " p=", p)
  if p == nothing
    cbi.bc_fcn(ya, yb, bca, bcb)
  else
    cbi.bc_fcn(ya, yb, p, bca, bcb)
  end
  l_bc && println(lio, lprefix, "bc result: bca", bca, " bcb=", bcb)
  return nothing
end

"""
  This is the boundary-conditions function given as callback
  to bvpm2 in the case with no unknown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_bc_cb(
     ya_len::Int64, ya_::Ptr{Float64}, 
     yb_len::Int64, yb_::Ptr{Float64},
     bca_len::Int64, bca_::Ptr{Float64},
     bcb_len::Int64, bcb_::Ptr{Float64}, cbi::CI) where CI<:Bvpm2_solve_cbi

  @assert ya_len == yb_len == cbi.no_odes
  @assert bca_len == cbi.no_left_bc
  @assert bcb_len == cbi.no_odes+cbi.no_par-cbi.no_left_bc
  ya = unsafe_wrap(Array, ya_, (ya_len,), own=false)
  yb = unsafe_wrap(Array, yb_, (yb_len,), own=false)
  bca = unsafe_wrap(Array, bca_, (bca_len,), own=false)
  bcb = unsafe_wrap(Array, bcb_, (bcb_len,), own=false)
  bvpm2_bc(ya, yb, bca, bcb, cbi)
  return nothing
end

function unsafe_bvpm2_bc_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_bc_cb, Cvoid,
    Tuple{Int64, Ptr{Float64}, Int64, Ptr{Float64},
     Int64, Ptr{Float64}, Int64, Ptr{Float64}, Ref{CI}})
end

"""
  This is the boundary-conditions function given as callback
  to bvpm2 in the case with unknown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_bcpar_cb(
     ya_len::Int64, ya_::Ptr{Float64}, 
     yb_len::Int64, yb_::Ptr{Float64},
     p_len::Int64, p_::Ptr{Float64},
     bca_len::Int64, bca_::Ptr{Float64},
     bcb_len::Int64, bcb_::Ptr{Float64}, cbi::CI) where CI<:Bvpm2_solve_cbi

   @assert ya_len == yb_len == cbi.no_odes
   @assert bca_len == cbi.no_left_bc
   @assert bcb_len == cbi.no_odes+cbi.no_par-cbi.no_left_bc
   @assert p_len == cbi.no_par
   ya = unsafe_wrap(Array, ya_, (ya_len,), own=false)
   yb = unsafe_wrap(Array, yb_, (yb_len,), own=false)
   p = unsafe_wrap(Array, p_, (p_len,), own=false)
   bca = unsafe_wrap(Array, bca_, (bca_len,), own=false)
   bcb = unsafe_wrap(Array, bcb_, (bcb_len,), own=false)
   bvpm2_bc(ya, yb, bca, bcb, cbi, p)
   return nothing
end

function unsafe_bvpm2_bcpar_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_bcpar_cb, Cvoid,
    Tuple{Int64, Ptr{Float64}, Int64, Ptr{Float64},
     Int64, Ptr{Float64},
     Int64, Ptr{Float64}, Int64, Ptr{Float64}, Ref{CI}})
end

"""
  This function calls `Dbc` saved in Bvpm2_solve_cbi.
  """
function bvpm2_Dbc(ya, yb, dya, dyb, cbi::CI, 
      p=nothing, dpa=nothing, dpb=nothing) where CI
  lprefix = cbi.bc_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_Dbc = l & LOG_JACBC > 0

  cbi.Dbc_count += 1
  l_Dbc && println(lio, lprefix, "called with ya=", ya, " yb=", yb, " p=", p)
  if p == nothing
    cbi.Dbc(ya, yb, dya, dyb)
  else
    cbi.Dbc(ya, yb, dya, dyb, p, dpa, dpb)
  end
  l_Dbc && println(lio, lprefix, "bc result: dya", dya, " dyb=", dyb, 
    " dpa=", dpa, " dpb=", dpb)
  return nothing
end

"""
  This is the derivative of the boundary-conditions given as callback
  to bvpm2 in the case with no unknown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_Dbc_cb(
  ya_len::Int64, ya_::Ptr{Float64}, yb_len::Int64, yb_::Ptr{Float64},
  dya_dim1::Int64, dya_dim2::Int64, dya_::Ptr{Float64},
  dyb_dim1::Int64, dyb_dim2::Int64, dyb_::Ptr{Float64}, 
  cbi::CI) where CI<:Bvpm2_solve_cbi

  @assert ya_len == yb_len == cbi.no_odes
  @assert dya_dim2 == dyb_dim2 == cbi.no_odes
  @assert dya_dim1 == cbi.no_left_bc
  @assert dyb_dim1 == cbi.no_odes+cbi.no_par-cbi.no_left_bc 

  ya = unsafe_wrap(Array, ya_, (ya_len,), own=false)
  yb = unsafe_wrap(Array, yb_, (yb_len,), own=false)
  dya = unsafe_wrap(Array, dya_, (dya_dim1, dya_dim2,), own=false)
  dyb = unsafe_wrap(Array, dyb_, (dyb_dim1, dyb_dim2,), own=false)
  bvpm2_Dbc(ya, yb, dya, dyb, cbi)
  return nothing
end

function unsafe_bvpm2_Dbc_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_Dbc_cb, Cvoid,
    Tuple{Int64, Ptr{Float64}, Int64, Ptr{Float64},
     Int64, Int64, Ptr{Float64}, Int64, Int64, Ptr{Float64}, Ref{CI}})
end


"""
  This is the derivative of the boundary-conditions given as callback
  to bvpm2 in the case with unknown parameters.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_bvpm2_Dbcpar_cb(
  ya_len::Int64, ya_::Ptr{Float64}, yb_len::Int64, yb_::Ptr{Float64},
  p_len::Int64, p_::Ptr{Float64},
  dya_dim1::Int64, dya_dim2::Int64, dya_::Ptr{Float64},
  dyb_dim1::Int64, dyb_dim2::Int64, dyb_::Ptr{Float64}, 
  dpa_dim1::Int64, dpa_dim2::Int64, dpa_::Ptr{Float64},
  dpb_dim1::Int64, dpb_dim2::Int64, dpb_::Ptr{Float64}, 
  cbi::CI) where CI<:Bvpm2_solve_cbi

  @assert ya_len == yb_len == cbi.no_odes
  @assert p_len == cbi.no_par
  @assert dya_dim2 == dyb_dim2 == cbi.no_odes
  @assert dpa_dim2 == dpb_dim2 == cbi.no_par
  @assert dya_dim1 == dpa_dim1 == cbi.no_left_bc
  @assert dyb_dim1 == dpb_dim1 == cbi.no_odes+cbi.no_par-cbi.no_left_bc 

  ya = unsafe_wrap(Array, ya_, (ya_len,), own=false)
  yb = unsafe_wrap(Array, yb_, (yb_len,), own=false)
  p = unsafe_wrap(Array, p_, (p_len,), own=false)
  dya = unsafe_wrap(Array, dya_, (dya_dim1, dya_dim2,), own=false)
  dyb = unsafe_wrap(Array, dyb_, (dyb_dim1, dyb_dim2,), own=false)
  dpa = unsafe_wrap(Array, dpa_, (dpa_dim1, dpa_dim2,), own=false)
  dpb = unsafe_wrap(Array, dpb_, (dpb_dim1, dpb_dim2,), own=false)
  bvpm2_Dbc(ya, yb, dya, dyb, cbi, p, dpa, dpb)
  return nothing
end

function unsafe_bvpm2_Dbcpar_cb_c(cbi::CI) where CI
  return cfunction(unsafe_bvpm2_Dbcpar_cb, Cvoid,
    Tuple{Int64, Ptr{Float64}, Int64, Ptr{Float64},
     Int64, Ptr{Float64},
     Int64, Int64, Ptr{Float64}, Int64, Int64, Ptr{Float64}, 
     Int64, Int64, Ptr{Float64}, Int64, Int64, Ptr{Float64}, Ref{CI}})
end


"""
       function bvpm2_solve(guess_obj::Bvpm2, rhs, bc, 
         opt::AbstractOptionsODE) -> (obj_out, retcode, stats)

  ## Right-hand side for the ODEs: `rhs`

  The function `rhs` must have the form:

       function rhs(x, y, f)              [no_par == 0]
       function rhs(x, y, p, f)           [no_par != 0]

  where

        x::Float, y::Vector{Float64}(no_odes), p::Vector{Float64}(no_par)
        f::Vector{Float64}(no_odes)

  Inside the function, the values of the right-hand side must be saved
  in `f`.

  ## Derivatives of right-hand side: `Drhs`

  The function `Drhs` is optional. If not given finite differences are
  used to approximate the derivatives. If `Drhs` is given it must have
  the form:

       function Drhs(x, y, dfdy)          [no_par == 0]
       function Drhs(x, y, p, dfdy, dfdp) [no_par != 0]

  where

        x::Float, y::Vector{Float64}(no_odes), p::Vector{Float64}(no_par)
        dfdy::Matrix{Float64}(no_odes, no_odes)
        dfdp::Matrix{Float64}(no_odes, no_par)

  Inside the function, the values of the derivatives must be saved
  in `dfdy` and `dfdp`.

  ## Boundary conditions: `bc`

  The function `bc` must have the form:

       function bc(ya, yb, bca, bcb)      [no_par == 0]
       function bc(ya, yb, p, bca, bcb)   [no_par != 0]

  where

        ya::Vector{Float64}(no_odes), yb::Vector{Float64}(no_odes), 
        p::Vector{Float64}(no_par),
        bca::Vector{Float64}(no_left_bc),
        bcb::Vector{Float64}(no_odes - no_left_bc)

  Inside the function, the values of the boundary conditions must be
  saved in `bca` and `bcb`.

  ## Derivatives of the boundary conditons: `Dbc`

  The function `Dbc` is optional. If not given finite differences are
  used to approximate the derivatives. If `Dbc` is given it must have
  the form:

       function Dbc(ya, yb, dya, dyb)                 [no_par == 0]
       function Dbc(ya, yb, dya, dyb, p, dpa, dpb)    [no_par != 0]

  where

        ya::Vector{Float64}(no_odes), yb::Vector{Float64}(no_odes), 
        p::Vector{Float64}(no_par),
        dya::Matrix{Float64}(no_left_bc, no_odes)
        dyb::Matrix{Float64}(no_odes+no_par-no_left_bc, no_odes)
        dpa::Matrix{Float64}(no_left_bc, no_par)
        dpb::Matrix{Float64}(no_odes+no_par-no_left_bc, no_par)

  Inside the function, the values of the derivatives of the boundary 
  conditions must be saved in `dya`, `dyb`, `dpa` and `dpb`.

  ## Options `opt`

  In `opt` the following options are used:
  
      ╔═════════════════╤══════════════════════════════════════════╤═════════╗
      ║  Option OPT_…   │ Description                              │ Default ║
      ╠═════════════════╪══════════════════════════════════════════╪═════════╣
      ║ RTOL            │ relative accuracy for solution.          │    1e-6 ║
      ║                 │ solution. Must be a scalar.              │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ METHODCHOICE    │ Choice for IVP-solvers:                  │       4 ║
      ║                 │ 2: Runge-Kutta method of order 2         │         ║
      ║                 │ 4: Runge-Kutta method of order 4         │         ║
      ║                 │ 6: Runge-Kutta method of order 6         │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ DIAGNOSTICOUTPUT│ diagnostic output for bvpm2:             │      -1 ║
      ║                 │   -1 : no output                         │         ║
      ║                 │    0 : only output if computation fails  │         ║
      ║                 │    1 : intermediate output               │         ║
      ║                 │    2 : full output                       │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ ERRORCONTROL    │ determines the error-estimation for      │       1 ║
      ║                 │ which RTOL is used:                      │         ║
      ║                 │    1 : defect                            │         ║
      ║                 │    2 : global error                      │         ║
      ║                 │    3 : defect and then global error      │         ║
      ║                 │    4 : linear combination of defect      │         ║
      ║                 │        and global error                  │         ║
      ╟─────────────────┼──────────────────────────────────────────┼─────────╢
      ║ SINGULARTERM    │ either nothing if the ODEs have no       │ nothing ║
      ║                 │ singular terms at the left boundary or   │         ║
      ║                 │ a constant (d,d) matrix for the          │         ║
      ║                 │ singular term.                           │         ║
      ╚═════════════════╧══════════════════════════════════════════╧═════════╝

  ## Return-Code `retcode`

  `retcode` can have to following values:

        <0: failure
        ≥0: computation successful
  """
function bvpm2_solve(guess_obj::Bvpm2, rhs, bc, 
  opt::AbstractOptionsODE; Drhs=nothing, Dbc=nothing)

  (lio,l,l_g,l_solver,lprefix) = solver_init("bvpm2", opt)

  if l_g
    println(lio, lprefix,
      "called with rhs=", rhs, "bc=", bc, "guess_obj=", guess_obj)
  end

  # singular-term, method, tol, trace error_control in opt

  details = bvpm2_check_state(guess_obj, (1,2))
  no_odes = details["no_odes"]
  no_par = details["no_par"]
  no_left_bc = details["no_left_bc"]
  tol = NaN; method = 4; trace =  -1; error_control = 1
  si_dim = 0; si_matrix = Matrix{Float64}(uninitialized, 0,0);

  OPT = nothing
  try
    OPT = OPT_RTOL
    tol = convert(Float64, getOption(opt, OPT, 1e-6))
    @assert isscalar(tol) && tol > 0.0

    OPT = OPT_METHODCHOICE
    method = convert(Int64, getOption(opt, OPT, 4))
    @assert method ∈ (2, 4, 6)

    OPT = OPT_DIAGNOSTICOUTPUT
    trace = convert(Int64, getOption(opt, OPT, -1))
    @assert trace ∈ (-1, 0, 1, 2)

    OPT = OPT_ERRORCONTROL
    error_control = convert(Int64, getOption(opt, OPT, 1))
    @assert error_control ∈ (1, 2, 3, 4)

    OPT = OPT_SINGULARTERM
    singular_term = getOption(opt, OPT, nothing)
    if singular_term ≠ nothing
      singular_term = convert(Matrix{Float64}, singular_term)
      @assert no_odes == size(singular_term, 1) == size(singular_term, 2)
      si_dim = no_odes
      si_matrix = copy(singular_term)
    end
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid", :opt, e))
  end
  
  error_ret = Vector{Float64}(uninitialized, 5)
  handle_out = [ C_NULL ]

  if bvpm2_global_cbi ≠ nothing
    throw(ArgumentErrorODE(string("The Fortran solver BVP_M-2 uses ",
      "global variables during the solution process. Hence this julia ",
      "module does not support concurrent/nested bvpm2_solve calls. Sorry."),
      :opt))
  end
  cbi = nothing; error = 0
  try
    global bvpm2_global_cbi = cbi = Bvpm2_solve_cbi(lio, l, 
      no_odes, no_par, no_left_bc,
      rhs, "unsafe_bvpm2_rhs: ", 0, 
      (Drhs ≠ nothing) ? Drhs : dummy_func, "unsafe_bvpm2_Drhs: ", 0,
      bc, "unsafe_bvpm2_b: ", 0,
      (Dbc ≠ nothing) ? Dbc : dummy_func, "unsafe_bvpm2_Dbc: ", 0)

    if l_solver
      println(lio, lprefix, "call Fortran bvpm2_solver ", 
        guess_obj.method_solve)
    end

    rhs_fcn_ptr = (no_par > 0) ? unsafe_bvpm2_rhspar_cb_c(cbi) :
                                 unsafe_bvpm2_rhs_cb_c(cbi)
    Drhs_fcn_ptr = (Drhs ≠ nothing) ? (
                     (no_par > 0) ? unsafe_bvpm2_Drhspar_cb_c(cbi) :
                                    unsafe_bvpm2_Drhs_cb_c(cbi)
                         ) : C_NULL
    bc_fcn_ptr =  (no_par > 0) ? unsafe_bvpm2_bcpar_cb_c(cbi) :
                                 unsafe_bvpm2_bc_cb_c(cbi)
    Dbc_fcn_ptr = (Dbc ≠ nothing) ? (
                     (no_par > 0) ? unsafe_bvpm2_Dbcpar_cb_c(cbi) :
                                    unsafe_bvpm2_Dbc_cb_c(cbi)
                         ) : C_NULL

    error = ccall(guess_obj.method_solve, Int64,
      (Ptr{Cvoid}, Ref{Ptr{Cvoid}},          # handle_guess, handle_out
       Ptr{Cvoid}, Ptr{Cvoid},               # rhs_fcn_ptr, bc_fcn_ptr
       Int64, Ptr{Float64},                  # si_dim, si_matrix
       Int64, Float64,                       # method, tol
       Ptr{Cvoid}, Ptr{Cvoid},               # dfdy_fcn_ptr, dbcdy_fcn_ptr
       Int64, Int64,                         # trace, error_control
       Int64, Ptr{Float64},                  # error_ret_len, error_ret
       Ref{Bvpm2_solve_cbi},                 # calls_pthrough
      ), 
      guess_obj.handle, handle_out,
      rhs_fcn_ptr, bc_fcn_ptr,
      si_dim, si_matrix,
      method, tol,
      Drhs_fcn_ptr, Dbc_fcn_ptr,
      trace, error_control,
      length(error_ret), error_ret,
      cbi)

    if l_solver
      println(lio, lprefix, "call Fortran bvpm2_solver ", 
        guess_obj.method_solve, " returned ",error)
    end
  finally
    global bvpm2_global_cbi = nothing
  end

  error ≠ 0 && throw(InternalErrorODE(string(
    "Internal error during solve call: error = ",error)))
  
  @assert handle_out[1] ≠ C_NULL
  obj_out = Bvpm2(handle_out[1])
  out_details = bvpm2_get_details(obj_out)
  stats = Dict(
    "no_rhs_calls"  => cbi.rhs_count,
    "no_jac_calls"  => cbi.Drhs_count,
    "no_bc_calls"   => cbi.bc_count,
    "no_Dbc_calls"  => cbi.Dbc_count,
    "cond"          => error_ret[1],
    "cerror"        => error_ret[2],
    "reerror"       => error_ret[3],
    "hoerror"       => error_ret[4],
    "dcerror"       => error_ret[5]
  )
  return obj_out, out_details["info"], stats
end

"""
       function bvpm2_extend(sol_obj::Bvpm2, anew, bnew, 
                yanew::Vector, ybnew::Vector; 
                p_new=[], max_num_subintervals=0)

  extends a solution (`state == 2`) to a new interval. Take the
  two given states `yanew` and `ybnew` as new states for the new
  interval. (If `anew ≥ a` then `yanew` is ignored. If `bnew ≤ b` then
  ybnew is ignored.)

  You can change the parameter guess also, by using `p_new` and you
  can change the maximal number of subintervals, too.

  sol_obj will be "terminated". After the call `sol_obj` will be
  in state 0. Use `bvpm2_copy` before, if you need the solution afterwards.

  A new Bvpm2-object `guess_obj` will be created an returned.
  """
function bvpm2_extend(sol_obj::Bvpm2, anew, bnew, 
         yanew::Vector, ybnew::Vector; 
         p_new::Vector=[], max_num_subintervals=0)
  
  details = bvpm2_check_state(sol_obj, (2,))
  try
    anew = convert(Float64, anew)
  catch e
    throw(ArgumentErrorODE("Cannot convert anew to Float64", :anew, e))
  end
  try
    bnew = convert(Float64, bnew)
  catch e
    throw(ArgumentErrorODE("Cannot convert bnew to Float64", :bnew, e))
  end
  no_odes = details["no_odes"]
  no_par = details["no_par"]
  try
    yanew = getVectorCheckLength(yanew, Float64, no_odes)
  catch e
    throw(ArgumentErrorODE("Cannot convert yanew", :yanew, e))
  end
  try
    ybnew = getVectorCheckLength(ybnew, Float64, no_odes)
  catch e
    throw(ArgumentErrorODE("Cannot convert ybnew", :ybnew, e))
  end
  if length(p_new) > 0
    try
      p_new = getVectorCheckLength(p_new, Float64, no_par)
    catch e
      throw(ArgumentErrorODE("Cannot convert p_new", :p_new, e))
    end
  end
  try
    max_num_subintervals = convert(Int64, max_num_subintervals)
  catch e
    throw(ArgumentErrorODE("Cannot convert max_num_subintervals to Float64",
      :max_num_subintervals, e))
  end
  result = [ C_NULL ]
  error = ccall(sol_obj.method_extend_s, Int64, 
    (Ptr{Cvoid}, Ref{Ptr{Cvoid}},     # handle_in, handle_out
     Float64, Float64,                # anew, bnew
     Int64, Ptr{Float64},             # yanew_len, yanew
     Int64, Ptr{Float64},             # ybnew_len, ybnew
     Int64, Ptr{Float64},             # p_len, p
     Int64,),                         # max_num_subintervals
    sol_obj.handle, result, 
    anew, bnew,
    length(yanew), yanew,
    length(ybnew), ybnew,
    length(p_new), p_new,
    max_num_subintervals)

  error ≠ 0 && throw(InternalErrorODE(string(
    "Internal error during extend call: error = ",error)))

  return Bvpm2(result[1])
end

"""
       function bvpm2_extend(sol_obj::Bvpm2, anew, bnew, order; 
                p_new::Vector=[], max_num_subintervals=0)

  extends a solution (`state == 2`) to a new interval using 
  constant (`order==0`) or linear (`order==1`) extrapolation.

  You can change the parameter guess also, by using `p_new` and you
  can change the maximal number of subintervals, too.

  sol_obj will be "terminated". After the call `sol_obj` will be
  in state 0. Use `bvpm2_copy` before, if you need the solution afterwards.

  A new Bvpm2-object `guess_obj` will be created an returned.
  """
function bvpm2_extend(sol_obj::Bvpm2, anew, bnew, order; 
         p_new::Vector=[], max_num_subintervals=0)
  details = bvpm2_check_state(sol_obj, (2,))
  try
    anew = convert(Float64, anew)
  catch e
    throw(ArgumentErrorODE("Cannot convert anew to Float64", :anew, e))
  end
  try
    bnew = convert(Float64, bnew)
  catch e
    throw(ArgumentErrorODE("Cannot convert bnew to Float64", :bnew, e))
  end
  try
    order = convert(Int64, order)
    @assert order ∈ (0, 1,)
  catch e
    throw(ArgumentErrorODE("Cannot convert order to Int64: 0 or 1", 
      :order, e))
  end
  no_odes = details["no_odes"]
  no_par = details["no_par"]
  if length(p_new) > 0
    try
      p_new = getVectorCheckLength(p_new, Float64, no_par)
    catch e
      throw(ArgumentErrorODE("Cannot convert p_new", :p_new, e))
    end
  end
  result = [ C_NULL ]
  error = ccall(sol_obj.method_extend_e, Int64,
    (Ptr{Cvoid}, Ref{Ptr{Cvoid}},     # handle_in, handle_out
     Float64, Float64,                # anew, bnew
     Int64,                           # order
     Int64, Ptr{Float64},             # p_len, p
     Int64,),                         # max_num_subintervals
    sol_obj.handle, result,
    anew, bnew,
    order,
    length(p_new), p_new,
    max_num_subintervals)

  error ≠ 0 && throw(InternalErrorODE(string(
    "Internal error during extend call: error = ",error)))

  return Bvpm2(result[1])
end

"""
       function evalSolution(sol::Bvpm2, x::Real, z::Vector{Float64}, 
         dz::Vector{Float64})

  Evaluates an already obtained solution `sol` at scalar point `x`.
  The values of the solution are saved in `z` which must be a 
  vector (of length d).
  For the vector `dz` there are two cases: If `dz` is a empty vector (of
  length 0) then the derivates of z are not calculated, otherwise
  `dz` has to be a vector (of length d) where the derivates are
  stored.
  """
function evalSolution(sol::Bvpm2, x::Real, z::Vector{Float64}, 
  dz::Vector{Float64}=Vector{Float64}(uninitialized, 0))
  
  details = bvpm2_check_state(sol, (2,))
  no_odes = details["no_odes"]
  length(z) ≠ no_odes && throw(ArgumentErrorODE(string(
    "Expected for z a vector of length ", no_odes, ", but found a vector ",
    "with length ",length(z)), :z))
  if length(dz) ≠ 0 && length(dz) ≠ no_odes
    throw(ArgumentErrorODE(string(
       "Expected for dz a vector of length 0 or ", no_odes, ", but found ",
       "a vector with length ",length(dz)), :dz))
  end
  error = ccall(sol.method_eval_s, Int64,
    (Ptr{Cvoid}, Float64,                  # handles, x
     Int64, Ptr{Float64},                  # z_len, z
     Int64, Ptr{Float64},                  # dz_len, dz
    ),
    sol.handle, convert(Float64, x),
    length(z), z,
    length(dz), dz)

  error ≠ 0 && throw(InternalErrorODE(string(
    "Internal error during solve call: error = ",error)))
  return nothing
end

"""
       function evalSolution(sol::Bvpm2, x::Real)  -> z
  
  Allocates vector `z` and calls `evalSolution(sol, x, z)`.
  """
function evalSolution(sol::Bvpm2, x::Real)
  details = bvpm2_check_state(sol, (2,))
  no_odes = details["no_odes"]
  z = Vector{Float64}(no_odes)
  evalSolution(sol, x, z)
  return z
end

"""
       function evalSolution(sol::Bvpm2, x::Vector{Float64}, 
         z::Matrix{Float64}, dz::Matrix{Float64})

  Evaluates an already obtained solution `sol` at scalar point `x`.
  The values of the solution are saved in `z` which must be a 
  vector (of length d).
  For the vector `dz` there are two cases: If `dz` is a empty vector (of
  length 0) then the derivates of z are not calculated, otherwise
  `dz` has to be a vector (of length d) where the derivates are
  stored.
  """
function evalSolution(sol::Bvpm2, x::Vector{Float64}, z::Matrix{Float64}, 
  dz::Matrix{Float64}=Matrix{Float64}(uninitialized, 0,0))
  
  details = bvpm2_check_state(sol, (2,))
  no_odes = details["no_odes"]
  x_len = length(x)
  size(z) ≠ (no_odes, x_len) && throw(ArgumentErrorODE(string(
    "Expected for z a (", no_odes, ", ", x_len, ") matrix, but found a (",
    size(z,1), ", ", size(z,2), ") matrix "), :z))
  if size(dz) ≠ (0,0) && size(dz) ≠ (no_odes, x_len) 
    throw(ArgumentErrorODE(string(
    "Expected for dz an empty matrix or a (", no_odes, ", ", x_len, ") ",
    "matrix, but found a (", size(dz,1), ", ", size(dz,2), ") matrix "), :dz))
  end
  error = ccall(sol.method_eval_v, Int64,
    (Ptr{Cvoid},                            # handle 
     Int64, Ptr{Float64},                  # x_len, x
     Int64, Int64, Ptr{Float64},           # z_dim1, z_dim2
     Int64, Int64, Ptr{Float64},           # dz_dim1, dz_dim2
    ),
    sol.handle,
    length(x), x,
    size(z, 1), size(z, 2), z, 
    size(dz, 1), size(dz, 2), dz
    )

  error ≠ 0 && throw(InternalErrorODE(string(
    "Internal error during solve call: error = ",error)))
  return nothing
end

"""
       function evalSolution(sol::Bvpm2, x::Vector{Float64}) -> z

  Allocates matrix `z` and calls `evalSolution(sol, x, z)`.
  """
function evalSolution(sol::Bvpm2, x::Vector{Float64})
  details = bvpm2_check_state(sol, (2,))
  no_odes = details["no_odes"]
  x_len = length(x)
  z = Matrix{Float64}(uninitialized, no_odes, x_len)
  evalSolution(sol, x, z)
  return z
end

"""  
  ## Compile BVP_M-2

  The julia ODEInterface tries to compile and link the solvers
  automatically at the build-time of this module. The following
  calls need only be done, if one uses a different compiler and/or if
  one wants to change/add some compiler options.

  The Fortran source code can be found at:
  
       http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml
  
  See `help_bvpm3_license` for the licsense information.
  
  ### Using `gfortran` and 64bit integers (Linux and Mac)
  
  Here is an example how to compile BVP_M-2 with `Float64` reals (and
  `Int64` integers with `gfortran`):

       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o bvp_la-2.o bvp_la-2.f
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o bvp_m-2.o bvp_m-2.f90
       gfortran -c -fPIC -fdefault-integer-8 
                -fdefault-real-8 -fdefault-double-8 
                -o bvp_m_proxy.o bvp_m_proxy.f90

  The last file `bvp_m_proxy.f90` is a Julia/C-Proxy and is part of this
  `ODEInterface` package.
  
  In order to get create a shared library (from the object file above) use
  one of the forms below (1st for Linux, 2nd for Mac):

       gfortran -shared -fPIC -o bvp_m_proxy.so 
                bvp_m_proxy.o bvp_m-2.o bvp_la-2.o
  """
function help_bvpm2_compile()
  return Docs.doc(help_bvpm2_compile)
end

"""
  Please call this `help_bvpm2_proxy` function to get an detailed description
  of the Fortran-Proxy written for BVP_M-2.
  """
function help_bvpm2_proxy()
  src_file = joinpath(guess_path_of_module(), "bvp_m_proxy.f90")
  markdown = ""
  fh = open(src_file, "r")
    while !eof(fh)
      line = readline(fh)
      if length(line) > 0 && line[1] == '!'
        markdown *= line[ 
          ((length(line) > 1 && line[2] == ' ') ? 3 : 2):end ] * "\n"
      else
        break
      end
    end
  close(fh)
  return Base.Markdown.parse(markdown)
end

"""
  # Licence

  BVP_SOLVER, Release 2, with global error estimation and control options.
  Copyright (c) 2012, Jason Boisvert, Paul Muir, Ray Spiteri.
  Jason Boisvert, Ray Spiteri, Department of Computer Science, University of Saskatchewan.
  Paul Muir, Mathematics and Computing Science, Saint Mary's University.
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
      * Neither Saint Mary's University nor Southern Methodist University nor 
        the names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.
 
  THIS SOFTWARE IS PROVIDED BY Jason Boisvert, Paul Muir, and Ray Spiteri ''AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL Paul Muir and Larry Shampine BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
  See documentation below for FUNCTION BVP_SOLVER for a description of the changes
  to the argument list for BVP_SOLVER. 
 
  Copyright (c) 2006, Paul Muir and Larry Shampine.
  Paul Muir, Mathematics and Computing Science, Saint Mary's University.
  Larry Shampine, Mathematics, Southern Methodist University.
  All rights reserved.
 
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:
      * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
      * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
      * Neither Saint Mary's University nor Southern Methodist University nor 
        the names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.
 
  THIS SOFTWARE IS PROVIDED BY Paul Muir and Larry Shampine ''AS IS'' AND ANY
  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL Paul Muir and Larry Shampine BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
  L.F. Shampine, P.H. Muir, H. Xu, A user-friendly Fortran BVP solver, 
  J. Numer. Anal. Indust. Appl. Math., 1, 2006, 201--217.
  """
function help_bvpm2_license()
  return Docs.doc(help_bvpm2_license)
end


# Add informations about solvers in global solverInfo-array.
push!(solverInfo,
  SolverInfo("bvpm2",
    "Boundary Value Problem Solver with defect and global error control",
    tuple( :OPT_RTOL,  :OPT_METHODCHOICE, :OPT_DIAGNOSTICOUTPUT,
           :OPT_ERRORCONTROL, :OPT_SINGULARTERM, ),
    tuple(
      SolverVariant("bvpm2_i64",
        "Bvpm2 with 64bit integers",
        DL_BVPM2,
        tuple("create_sol_wrapper_c", 
              "copy_sol_wrapper_c",
              "terminate_sol_wrapper_c",
              "destroy_sol_wrapper_c", 
              "show_sol_wrapper_c", 
              "get_sol_wrapper_details_c", 
              "get_sol_wrapper_x_c",
              "init_guess1_c", 
              "init_guess2_c", 
              "init_guess3_c",
              "solve_c",
              "eval_s_sol_c",
              "eval_v_sol_c",
              "get_sol_wrapper_params_c",
              "extend_sol_s_c",
              "extend_sol_e_c",
              ))
    ),
    help_bvpsol_compile,
    help_bvpsol_license,
  )
)

# vim:syn=julia:cc=79:fdm=indent:
