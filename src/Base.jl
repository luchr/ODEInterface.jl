# Basic functions for ODEInterface

"""
  supported (signed) Integer types for Fortran codes.
  """
FortranInt = Union{Int32,Int64}

"""
  dummy function returning nothing.
  """
function dummy_func()
  return nothing
end

"""
  tests if `cand` is a number.
  """
function isscalar(cand)
  return isa(cand,Base.Number)
end

"""Bitmask: log nothing."""
const LOG_NOTHING    = UInt64(0)

"""Bitmask: log some general info (esp. for main call)."""
const LOG_GENERAL    = UInt64(1)<<0

"""Bitmask: log calls to right-hand side."""
const LOG_RHS        = UInt64(1)<<1

"""Bitmask: log arguments passed to C/Fortran solvers."""
const LOG_SOLVERARGS = UInt64(1)<<2

"""Bitmask: log calls to solout function."""
const LOG_SOLOUT     = UInt64(1)<<3

"""Bitmask: log calls to julia output function."""
const LOG_OUTPUTFCN  = UInt64(1)<<4

"""Bitmask: log calls to eval_sol_fcn function."""
const LOG_EVALSOL    = UInt64(1)<<5

"""Bitmask: log calls to mass function."""
const LOG_MASS       = UInt64(1)<<6

"""Bitmask: log calls to jacobian function."""
const LOG_JAC        = UInt64(1)<<7

"""Bitmask: log calls to boundary condition function."""
const LOG_BC         = UInt64(1)<<8

"""Bitmask: during boundary value problems: 
  log calls to initial value solver."""
const LOG_BVPIVPSOL  = UInt64(1)<<9

"""Bitmask: log calls to right-hand side derivative function."""
const LOG_RHSDT      = UInt64(1)<<10

"""Bitmask: log everything."""
const LOG_ALL        = UInt64(0xFFFFFFFFFFFFFFFF)

"""
  macro for importing all the LOG bitmasks.
  """
macro import_LOG()
  :(
    using ODEInterface: LOG_NOTHING, LOG_GENERAL, LOG_RHS,
                        LOG_SOLVERARGS, LOG_SOLOUT, LOG_OUTPUTFCN,
                        LOG_EVALSOL, LOG_MASS, LOG_JAC, LOG_BC, 
                        LOG_BVPIVPSOL, LOG_RHSDT,
                        LOG_ALL
  )
end

"""
       function getVectorCheckLength(vec,T::DataType,d::Integer,copy=true)
                         -> Vector{T}
  
  try to convert to `Vector{T}` and checks given length.
  If the `docopy` argument is `true` then the return value will
  always be a different object than `vec`: If `convert` didn't need to
  create a copy then this is done by this function.
  
  throws ArgumentErrorODE this is not possible.
  """
function getVectorCheckLength(vec,T::DataType,d::Integer,docopy=true)
  result = nothing
  try
    result = convert(Vector{T},vec)
  catch e
    throw(ArgumentErrorODE(
      string("Cannot convert vec with type ",typeof(vec),
             " to Vector{",T,"}"),:vec,e))
  end
  length(result)!=d && throw(ArgumentErrorODE(
    string("vec has wrong length: expected $d found ",length(result)),:vec))
  if docopy && result ≡ vec
    result = copy(vec)
  end
  return result
end

"""
       function getMatrixCheckSize(mat,T::DataType,
                    m::Integer,n::Integer,docopy=true) -> Matrix{T}
  
  try to convert to `Matrix{T}' and checks the size.
  if the `docopy` argument is `true` then the return value will
  always be a different object than `mat`: If `convert` didn't need
  to create a copy then this is done by this function.
  
  throws ArgumentErrorODE this is not possible.
  """
function getMatrixCheckSize(mat,T::DataType,m::Integer,n::Integer,docopy=true)
  result = nothing
  try
    result = convert(Matrix{T},mat)
  catch e
    throw(ArgumentErrorODE(
      string("Cannot convert mat with type ",typeof(mat),
             "to Matrix{",T,"}"),:mat,e))
  end
  size(mat,1) ≠ m && throw(ArgumentErrorODE(
    string("mat has wrong number of rows: expedted $m found ",size(mat,1)),
    :mat))
  size(mat,2) ≠ n && throw(ArgumentErrorODE(
    string("mat has wrong number of columns: expedted $n found ",size(mat,2)),
    :mat))
  if docopy && result ≡ mat
    result = copy(mat)
  end
  return result
end

if !isdefined(Base, :unsafe_wrap)
  """Define unsafe_wrap for old julia versions.
    """
  function unsafe_wrap(target, p, dims, own)
    return pointer_to_array(p, dims, own)
  end
end

"""
  Type encapsulating all required data for ODE-Solver-Callbacks.
  
  For further explanation see, `GlobalCallInfoDict`
  """
@ABSTRACT(ODEinternalCallInfos,Any)

# vim:syn=julia:cc=79:fdm=indent:
