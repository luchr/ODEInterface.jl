# Basic functions for ODEInterface

import Base: dump

function dump(io::IO, x::Tuple, n::Int, indent)
  println(io,typeof(x)," len ",length(x))
  if n>0
    i = 1
    for elem in x
      print(io,indent,"  ", i, ": ")
      dump(io,elem,n-1,string(indent,"  "))
      if i > 10
        println(io,indent,"  ...")
        break
      end
      i += 1
    end
  end
end

"""
  dummy function returning nothing.
  """
function dummy_func()
  return nothing
end

"""
      function uniqueToken() -> Array{UInt64}
  
  Generate an unique token `T` with unique `UInt64` token-value `T[1]`.
  
  As long as the return value of this function will exist (i.e. is not garbage
  collected) this function will never return the same token value.
  
  Note:
  For the implementation no token-pool is needed. The idea is: just
  return an newly created `UInt64` Array (of length 1) and save
  the `object_id` in this array. As long as this object/array is not
  garbage collected, every newly created token/array will have a different 
  `object_id`.
  """
function uniqueToken()
  token = Array{UInt64}(1)
  token[1] = object_id(token)
  return token
end

"""
       function packUInt64ToVector!{T}(vec::Vector{T}, value::UInt64)

  stores a UInt64 in the given integer Vector. sizeof(T)∊{4,8}
  
  see unpackUInt64FromVector or unpackUInt64FromPtr for the opposite action.
  """
function packUInt64ToVector!{T}(vec::Vector{T}, value::UInt64)
  if 4==sizeof(T)
    vec[1] = reinterpret(T,UInt32(value>>32))
    vec[2] = reinterpret(T,UInt32(value&0xffffffff))
  elseif 8==sizeof(T)
    vec[1] = reinterpret(T,value)
  else
    throw(InternalErrorODE("Can only handle sizeof(T)∊{4,8}"))
  end
  return nothing
end

"""
       function unpackUInt64FromVector{T}(vec::Vector{T}) -> UInt64
  
  retrieves a UInt64 from the given integer Vector. sizeof(T)∊{4,8}
  
  see packUInt64ToVector for the opposite action.
  """
function unpackUInt64FromVector{T}(vec::Vector{T})
  if 4==sizeof(T)
    return UInt64(reinterpret(UInt32,vec[1]))<<32  |
           UInt64(reinterpret(UInt32,vec[2]))
  elseif 8==sizeof(T)
    return UInt64(reinterpret(UInt64,vec[1]))
  else
    throw(InternalErrorODE("Can only handle sizeof(T)∊{4,8}"))
  end
end

"""
       function unpackUInt64FromPtr{T}(p::Ptr{T}) -> UInt64
  
  retrieves a UInt64 from the given (integer) pointer. sizeof(T)∊{4,8} 
  
  see packUInt64ToVector for the opposite action.
  """
function unpackUInt64FromPtr{T}(p::Ptr{T})
  if 4==sizeof(T)
    return UInt64(reinterpret(UInt32,unsafe_load(p,1)))<<32 |
           UInt64(reinterpret(UInt32,unsafe_load(p,2)))
  elseif 8==sizeof(T)
    return UInt64(reinterpret(UInt64,unsafe_load(p,1)))
  else
    throw(InternalErrorODE("Can only handle sizeof(T)∊{4,8}"))
  end
end

"""
         function int2logstr(no::Unsigned) -> UTF8String
    
    Generate Log-String for integer (token-)number.
  """
int2logstr(no::Unsigned) = "⟬$(hex(no))⟭"

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

"""Bitmask: log everything."""
const LOG_ALL        = UInt64(0xFFFFFFFFFFFFFFFF)

"""
  macro for importing all the LOG bitmasks.
  """
macro import_LOG()
  :(
    using ODEInterface: LOG_NOTHING, LOG_GENERAL, LOG_RHS,
                        LOG_SOLVERARGS, LOG_SOLOUT, LOG_OUTPUTFCN,
                        LOG_EVALSOL, LOG_ALL
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


# vim:syn=julia:cc=79:fdm=indent:
