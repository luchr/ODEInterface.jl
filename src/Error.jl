# Errors and Exceptions for ODEInterface

import Base: showerror

"""macro for importing the ODE Exceptions."""
macro import_exceptions()
    return :(
        using ODEInterface: WrappedODEException, ArgumentErrorODE,
            OutputErrorODE, SolverODEnotLoaded,
            FunctionCallNotSupported, FeatureNotSupported,
            StateErrorODE,
            InternalErrorODE
    )
end

"""
The ancestor for all wrapped exceptions in ODEInterface.

Required fields: msg, error
"""
abstract type WrappedODEException <: Base.WrappedException end

function showerror(io::IO, e::WrappedODEException)
    println(io, e.msg)
    return if e.error !== nothing
        println(io, "Wrapped exception:")
        showerror(io, e.error)
    end
end

"""
This error indicates that one input argument is invalid.

This is a WrappedException: If the invalidity of the argument
was detected by some error/exception then, this initial
error/exception can be found in the `error` field.
"""
mutable struct ArgumentErrorODE <: WrappedODEException
    msg::AbstractString
    argname::Symbol
    error
end

function ArgumentErrorODE(msg, argname)
    return ArgumentErrorODE(msg, argname, nothing)
end

function showerror(io::IO, e::ArgumentErrorODE)
    println(io, e.msg)
    e.argname !== nothing && println(io, string("in argument ", e.argname))
    return if e.error !== nothing
        println(io, "Wrapped exception:")
        showerror(io, e.error)
    end
end

"""
This error indicates that a function returned invalid output.
"""
mutable struct OutputErrorODE <: WrappedODEException
    msg::AbstractString
    func
    error
end

function OutputErrorODE(msg, func)
    return OutputErrorODE(msg, func, nothing)
end

function showerror(io::IO, e::OutputErrorODE)
    println(io, e.msg)
    e.func !== nothing && println(io, string("function ", e.func))
    return if e.error !== nothing
        println(io, "Wrapped exception:")
        showerror(io, e.error)
    end
end

"""
This error indicates that a Fortran/C-solver is not loaded.
"""
mutable struct SolverODEnotLoaded <: WrappedODEException
    msg::AbstractString
    error
end

function SolverODEnotLoaded(msg)
    return SolverODEnotLoaded(msg, nothing)
end

"""
This error indicates that a function was called at a time, where
this is not possible.
"""
mutable struct FunctionCallNotSupported <: WrappedODEException
    msg::AbstractString
    error
end

function FunctionCallNotSupported(msg)
    return FunctionCallNotSupported(msg, nothing)
end

"""
This error indicates that a requested feature is not supported or
is not possible.
"""
mutable struct FeatureNotSupported <: WrappedODEException
    msg::AbstractString
    error
end

function FeatureNotSupported(msg)
    return FeatureNotSupported(msg, nothing)
end

"""
This error indicates that an object is in the wrong state, e.g.
is not initialized.
"""
mutable struct StateErrorODE <: WrappedODEException
    msg::AbstractString
    error
end

function StateErrorODE(msg)
    return StateErrorODE(msg, nothing)
end

"""
This error indicates an internal error.
"""
mutable struct InternalErrorODE <: WrappedODEException
    msg::AbstractString
    error
end

function InternalErrorODE(msg)
    return InternalErrorODE(msg, nothing)
end


# vim:syn=julia:cc=79:fdm=indent:
