# Errors and Exceptions for ODEInterface

import Base: showerror

"""macro for importing the ODE Exceptions."""
macro import_exceptions()
  :(
    using ODEInterface: WrappedODEException, ArgumentErrorODE,
                        OutputErrorODE, SolverODEnotLoaded, 
                        FunctionCallNotSupported, FeatureNotSupported,
                        InternalErrorODE
  )
end

"""
  The ancestor for all wrapped exceptions in ODEInterface.
  
  Required fields: msg, error
  """
abstract WrappedODEException <: Base.WrappedException

function showerror(io::IO,e::WrappedODEException)
  println(io,e.msg)
  if e.error!=nothing
    println(io,"Wrapped exception:")
    showerror(io,e.error)
  end
end

"""
  This error indicates that one input argument is invalid.
  
  This is a WrappedException: If the invalidity of the argument
  was detected by some error/exception then, this initial
  error/exception can be found in the `error` field.
  """
type ArgumentErrorODE <: WrappedODEException
  msg      :: AbstractString
  argname  :: Symbol
  error
end

function ArgumentErrorODE(msg,argname)
  ArgumentErrorODE(msg,argname,nothing)
end

function showerror(io::IO,e::ArgumentErrorODE)
  println(io,e.msg)
  e.argname!=nothing && println(io,string("in argument ",e.argname)) 
  if e.error!=nothing
    println(io,"Wrapped exception:")
    showerror(io,e.error)
  end
end

"""
  This error indicates that a function returned invalid output.
  """
type OutputErrorODE <: WrappedODEException
  msg      :: AbstractString
  func     :: Function
  error
end

function OutputErrorODE(msg,func)
  OutputErrorODE(msg,func,nothing)
end

function showerror(io::IO,e::OutputErrorODE)
  println(io,e.msg)
  e.func!=nothing && println(io,string("function ",e.func)) 
  if e.error!=nothing
    println(io,"Wrapped exception:")
    showerror(io,e.error)
  end
end

"""
  This error indicates that a Fortran/C-solver is not loaded.
  """
type SolverODEnotLoaded <: WrappedODEException
  msg      :: AbstractString
  error
end

function SolverODEnotLoaded(msg)
  SolverODEnotLoaded(msg,nothing)
end

"""
  This error indicates that a function was called at a time, where
  this is not possible.
  """
type FunctionCallNotSupported <: WrappedODEException
  msg      :: AbstractString
  error
end

function FunctionCallNotSupported(msg)
  FunctionCallNotSupported(msg,nothing)
end

"""
  This error indicates that a requested feature is not supported or 
  is not possible.
  """
type FeatureNotSupported <: WrappedODEException
  msg      :: AbstractString
  error
end

function FeatureNotSupported(msg)
  FeatureNotSupported(msg,nothing)
end

"""
  This error indicates an internal error.
  """
type InternalErrorODE <: WrappedODEException
  msg      :: AbstractString
  error
end

function InternalErrorODE(msg)
  InternalErrorODE(msg,nothing)
end


# vim:syn=julia:cc=79:fdm=indent:

