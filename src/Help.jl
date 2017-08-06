# Some help texts for the ODEInterface

"""macro for importing all the help-functions."""
macro import_help()
  quote
    using ODEInterface: help_overview, help_callsolvers, help_options, 
          help_outputfcn, help_solversupport, help_install, 
          help_loadsolvers, help_specialstructure,
          help_internals
    @ODEInterface.import_dopri5_help
    @ODEInterface.import_odex_help
    @ODEInterface.import_radau5_help
    @ODEInterface.import_radau_help
    @ODEInterface.import_seulex_help
    @ODEInterface.import_rodas_help
    @ODEInterface.import_ddeabm_help
    @ODEInterface.import_ddebdf_help
    @ODEInterface.import_bvpsol_help
    @ODEInterface.import_colnew_help
  end
end


# takebuf_string deprecation: 358c4419
if VERSION < v"0.6.0-dev+1254"
  buf2str(buf) = takebuf_string(buf)
else
  buf2str(buf) = String(take!(buf))
end

using Base.Markdown

"""
  # Overview

  ## Importing all help topics
  
  You can use
  
       using ODEInterface
       @ODEInterface.import_help
  
  to have all `help_...` commands imported; so you can use
  
       ?help_overview
       help_overview()
  
  ## Help topics
  
       help_install          requirements, installation, compiling the solvers
       help_solversupport    supported ODE solvers
       help_loadsolvers      loading the ODE solvers
       help_callsolvers      how to call the ODE solvers
       help_options          how to use parameters/options for solvers
       help_outputfcn        how to use "output functions", dense output
       help_specialstructure support for problems with "special structure"
       
       help_internals        some internal/developer information
  
  ## Help for each solver
  
  Each solver has its own help page. Just look at the documentation of
  `dopri5`, `dop853`, `odex`, `radau5`, `radau`, `rodas`, `seulex`,
  `bvpsol`, `ddeabm` and `Bvpm2`.
  """
function help_overview()
  return Docs.doc(help_overview)
end

"""
       function help_solversupport()
  
  This function (when called) produces a (markdown) object with
  informations about the supported solvers.
  """
function help_solversupport()
  io = IOBuffer()
  write(io,"## Supported solvers\n\n")
  for solver in solverInfo
    write(io,"### ",solver.name,"\n\n",solver.description,"\n\n");
    for variant in solver.variants
      write(io,"* ",variant.name,": ",variant.description,"\n");
    end
    write(io,"\nHelp for compiling: ",string(solver.help_compile));
    write(io,"\n\nLicense: ",string(solver.help_license),"\n\n");
  end
  write(io,"## Load status\n\n",
           "     ╔═════════════╤═══════════════╤══════════",
           "╤═══════════════════════════════╗\n",
           "     ║    Name     │   libname     │  loaded  ",
           "│       method                  ║\n",
           "     ╠═════════════╪═══════════════╪══════════",
           "╪═══════════════════════════════╣\n")
  first = Vector{Bool}(4)
  first_solver = true
  for solver in solverInfo
    if !first_solver
      write(io,"     ╟─────────────┼───────────────┼──────────",
               "┼───────────────────────────────╢\n")
    end
    first_solver = false
    first[1] = true
    first_variant = true
    for variant in solver.variants
      if !first_variant
        write(io,"     ║             ├───────────────┼──────────",
                 "┼───────────────────────────────╢\n")
      end
      loaded = true
      try
        v = dlSolversInfo[variant.libname]
        @assert v.error == nothing
        @assert v.libhandle ≠ C_NULL
      catch 
        loaded = false
      end
      
      first[2] = first[3] = true
      for method in variant.methods
        found = true
        try
          index = findfirst( x -> x.generic_name == method, 
                             dlSolversInfo[variant.libname].methods )
          m = dlSolversInfo[variant.libname].methods[index]
          @assert m.error == nothing
          @assert m.method_ptr ≠ C_NULL
        catch e
          found = false
        end
        method_name = method
        if length(method_name) > 22
          method_name = string(method[1:19],"...")
        end
        write(io,"     ║ ",rpad(first[1]?solver.name:"",12),"│ ",
                 rpad(first[2]?variant.libname:"",14)      ,"│ ",
                 rpad(first[3]?(loaded?"yes":" no"):"",9)  ,"│ ",
                 rpad(method_name,23),found?"  OK  ":"FAILED",
                 " ║\n")
        first[1] = first[2] = first[3] = false
      end
      first[1] = false
      first_variant = false
    end
  end
  write(io,"     ╚═════════════╧═══════════════╧══════════",
           "╧═══════════════════════════════╝\n\n",
           "For more informations: see result of `loadODESolvers`.\n")
  return Markdown.parse(buf2str(io))
end

"""
  ## Installation
  
  Because this module is an *interface* for C-/Fortran-ODE-solvers, you
  need to compile/link the solvers before you can use them.
  
  All solvers are dynamically loaded (at julia runtime),
  see `loadODESolvers`. Therefor a shared library is needed for every
  solver, i.e. you have to compile each solver and create a shared library.
  
  See `help_solversupport` for a list with supported solvers and for
  further informations how to compile/link the solvers.

  This module has its own build script `deps/build.jl` which tries
  to compile and link the shared libraries automatically.
  """
function help_install()
  return Docs.doc(help_install)
end

"""
  ## Loading the solvers
  
  All ODE solvers are dynamically loaded. See `help_install` for 
  informations how to create such shared libraries for the solvers.
  
  Before using a solver for the 1st time, it has to be loaded by
  a call of `loadODESolvers`.
  """
function help_loadsolvers()
  return Docs.doc(help_loadsolvers)
end

"""
  ## Calling the Solvers
  
  There are two ways of calling the solvers.
  
  1. A calling convention close to the original Fortran-call,
     trying to provide/expose all the features the Fortran-codes have.
  1. A simplified version, closer to odecalls like in MATLAB.
  
  ### The full-featured calling-method
  
  All ODE-solvers have the same calling convention:
  
      (t,x,retcode,stats) = 
          odesolver(rhs, t0::Real, T::Real,
                    x0::Vector, opt::AbstractOptionsODErhs)
  
      function rhs(t::Float64,x::Vector{Float64}) -> Vector{Float64}
          if OPT_RHS_CALLMODE == RHS_CALL_RETURNS_ARRAY
  
      function rhs(t::Float64,x::Vector{Float64},dx::Vector{Float64}) 
          if OPT_RHS_CALLMODE == RHS_CALL_INSITU
  
  The input arguments are:
  
  1. a julia function `rhs` for evaluating the right-hand side of the ODE,
     see below.
     It's OK to return something, that `convert` can transform
     to a `Vector{Float64}`.
  1. the initial time `t0`. `(t0,x0)` is the initial value of the 
     initial value problem.
  1. the final time `T`.
  1. the initial state `x0`. `(t0,x0)` is the initial value of the 
     initial value problem.
  1. further parameters/options in `opt` for the solver and for the interface. 
     There is a separate section for the explanation of the options, see
     help_options.
  
  The output arguments are:
  
  1. `t` the *last* time for which the solution has been computed 
     (if the whole computation was successfull, then `t==T`)
  1. `x` the numerical solation at time `t`
  1. `retcode` the return code of the solver (interpretation is solver dependent)
  
  There are two possible ways to provide the Julia right-hand side:
  
      function rhs(t::Float64,x::Vector{Float64}) -> Vector{Float64}
  
  This is used, if `OPT_RHS_CALLMODE == RHS_CALL_RETURNS_ARRAY`. So you can
  use anonymous functions like `(t,x) -> x` as right-hand sides. But this
  form has a price: every time the right-hand side is called, a temporary
  Array is created (the result). The form
  
      function rhs(t::Float64,x::Vector{Float64},dx::Vector{Float64}) 
                   -> nothing
  
  used if `OPT_RHS_CALLMODE == RHS_CALL_INSITU` does not have this problem.
  But the right-hand side must be a function filling in the values of `x'`
  in `dx`.
  
  ### The simplified version
  
  There is a simplified calling convention (using the methods above) to
  provide a method like odecalls in MATLAB, 
  see `odecall`.
  """
function help_callsolvers()
  return Docs.doc(help_callsolvers)
end

"""
  ## The `opt` Argument: the Options

  All options are handled by `OptionsODE`. See `OptionsODE` how to
  query, set and change options.
  
  There are the following classes of options.
  
  1. Options for this ODEInterface (common for all solvers)
  1. Options for the ODE solvers
  
  ### Options of this ODEInterface
  
  * `OPT_RHS_CALLMODE`:
    There are two possible ways to call the Julia right-hand side: 
    `RHS_CALL_RETURNS_ARRAY` and `RHS_CALL_INSITU`, see
    `help_callsolvers` for an explanation.
    difference.
  * `OPT_LOGIO`:
    This option sets the `IO` that is used for logging messages
  * `OPT_LOGLEVEL`:
    This is a bitmask for activating different logging messages. 
    The following bitmasks are available.
  
           LOG_NOTHING     log nothing
           LOG_GENERAL     log some general information, 
                           especially the main julia call of the solver
           LOG_RHS         log all calls of the right-hand side
           LOG_SOLVERARGS  log the arguments for the C-/Fortran-calls
                           before and after the call
           LOG_OUTPUTFCN   log calls of the julia output function
           LOG_SOLOUT      log calls of the solution output routine
           LOG_EVALSOL     log calls of the eval_sol_fcn
           LOG_MASS        log call(s) of the mass function
           LOG_JAC         log calls of the jacobian function of RHS
           LOG_BC          log calls of the boundary condition function
           LOG_BVPIVPSOL   log (during boundary value problems) calls to
                           initial value solvers
           LOG_RHSDT       log calls of the right-hand side time-derivative
           LOG_JACBC       log calls of the jacobian of the boundary condition
           LOG_GUESS       log calls to the guess function
           LOG_ALL         all of the above
  
  ### Options for the solvers
  
  Different solvers support different options. All the options a solver
  supports are listed in the help-command of the specific solver, e.g.
  `help_dopri5`.
  
  To get an overview, what options are supported by what solvers,
  call `ODEInterface.help_options()`.
  """
function help_options()
  # sort solvers alphabetically and give each solver a number
  solvers = sort(solverInfo,by=x -> x.name)
  # get for every option all the solvers that use this option
  opt_usage = Dict{Symbol,Set{SolverInfo}}()
  for solver in solvers
    for opt in solver.supported_opts
      if !haskey(opt_usage,opt)
        opt_usage[opt] = Set{SolverInfo}()
      end
      push!(opt_usage[opt],solver)
    end
  end
  # produce output: alphab. order for options
  io = IOBuffer()
  write(io,"## Solvers\n\n")
  for solver in solvers;       write(io,"1. ",solver.name,"\n"); end
  write(io,"\n## Options used by Solvers\n\n")

  opts = collect(keys(opt_usage))
  sort!(opts)

  write(io,"      ╔══════════════════════")
  for k in 1:length(solvers);  write(io,"╤══");                  end
  write(io,"╗\n")

  write(io,"      ║        Option        ")
  for k in 1:length(solvers);  write(io,"│",lpad(string(k),2));  end
  write(io,"║\n")

  write(io,"      ╠══════════════════════")
  for k in 1:length(solvers);  write(io,"╪══");                  end
  write(io,"╣\n")

  for opt in opts
    write(io,"      ║ ",rpad(opt,21))
    for solver in solvers; 
      write(io,"│ ",(solver ∈ opt_usage[opt])?"✔":" ")
    end
    write(io,"║\n")
  end

  write(io,"      ╚══════════════════════")
  for k in 1:length(solvers);  write(io,"╧══");                  end
  write(io,"╝\n")

  return Markdown.parse(buf2str(io))
end

"""
  ## `OPT_OUTPUTMODE` 
  
  This option determines if the `OPT_OUTPUTFCN` is
  called, and if dense output (the `eval_sol_fcn`) is prepared/supported.
  
  * `OUTPUTFCN_NEVER`: don't call the output function
  * `OUTPUTFCN_WODENSE`: call the output function, but `eval_sol_fcn`
    is not used
  * `OUTPUTFCN_DENSE`: call the output function and prepare `eval_sol_fcn`

  ## `OPT_OUTPUTFCN` 
  
       function outputfcn(reason::OUTPUTFCN_CALL_REASON,
        told::Float64,t::Float64, x::Vector{Float64},eval_sol_fcn::Function,
        extra_data::Dict)  -> OUTPUTFCN_RETURN_VALUE
  
  A (julia) function that is called 
  
  1. at beginning of the solution process with
     `reason == OUTPUTFCN_CALL_INIT`, `told=t0`, `t`=`T`, `x=x0`,
     `eval_sol_fcn` a dummy function throwing an error if called,
     `extra_data` a `Dict` persistent until the last call of the output 
     function. The return value is *ignored*.
  1. after every successfull integration step with
     `reason == OUTPUTFCN_CALL_STEP`, `[told,t]` the time interval of the
     last step, `x` the numerical solution at time `t`,
     `eval_sol_fcn` a function to evaluate the solution in `t1 ∊ [told,t]`, if
     requested by `OPT_OUTPUTMODE`, otherwise a call to this function will
     result in an error.
  1. at the end (after the last step) with
     `reason == OUTPUTFCN_CALL_DONE`.
     The return value is *ignored*.

  With `eval_sol_fcn`
  
          function eval_sol_fcn(t1::Float64) -> Vector{Float64}
  
  the numerical solution can be evaluted in the
  time interval `[told,t]` (if `OPT_OUTPUTMODE == OUTPUTFCN_DENSE`).
  
  If supported by the solver, the numerical solution may be changed
  in the `outputfcn` (if `reason == OUTPUTFCN_CALL_STEP`) and the solver
  continues the process with the changed solution.
  The return value `OUTPUTFCN_RET_CONTINUE_XCHANGED` indicates 
  this. `OUTPUTFCN_RET_CONTINUE` tells the solver to continue (without
  changes in `x`) and `OUTPUTFCN_RET_STOP` tells the solver to stop the
  solver.
  """
function help_outputfcn()
  return Docs.doc(help_outputfcn)
end

"""
  ## Special Structure
  
  Some solvers (e.g. radau5 and radau) supports ODEs with a
  "special structure". In this context an ODE has special structure 
  if there exists M1>0 and M2>0 and M1 = M⋅M2 for some integer M
  and
  
                   x'(k) = x(k+M2)   for all k = 1,…,M1            (*)
  
  There are the options `OPT_M1` and `OPT_M2` to tell the solvers, if the
  ODE has this special structure. In this case only the non-trivial parts
  of the right-hand side, of the mass matrix and of the Jacobian matrix 
  had to be supplied.
  
  If an ODE has the special structure, then the right-hand side 
  for `OPT_RHS_CALLMODE == RHS_CALL_RETURNS_ARRAY` has to return
  a vector of length d-M1 because the right-hand side for the first M1 entries
  are known from (*). For `OPT_RHS_CALLMODE == RHS_CALL_INSITU` the
  right-hand side gets (a reference) to the array of length d, but only
  the last d-M1 components need to be filled in.
  
  The mass matrix has the form
  
           ⎛ 1    │    ⎞ ⎫
           ⎜  ⋱   │    ⎟ ⎪
           ⎜   ⋱  │    ⎟ ⎬ M1
           ⎜    ⋱ │    ⎟ ⎪
      M =  ⎜     1│    ⎟ ⎭
           ⎜──────┼────⎟
           ⎜      │    ⎟ ⎫
           ⎜      │ M̃  ⎟ ⎬ d-M1
           ⎝      │    ⎠ ⎭
  
  Then there has to be only the (d-M1)×(d-M1) matrix M̃ in `OPT_MASSMATRIX`.
  Of course, M̃ can be banded. Then `OPT_MASSMATRIX` is a `BandedMatrix` with
  lower bandwidth < d-M1.
  
  If an ODE has the special structure, then the Jacobian matrix has the form

           ⎛0   1      ⎞ ⎫
           ⎜ ⋱   ⋱     ⎟ ⎪
           ⎜  ⋱   ⋱    ⎟ ⎬ M1
           ⎜   ⋱   ⋱   ⎟ ⎪
      J =  ⎜    0   1  ⎟ ⎭
           ⎜───────────⎟
           ⎜           ⎟ ⎫
           ⎜      J̃    ⎟ ⎬ d-M1
           ⎝           ⎠ ⎭
  
  Then there has to be only the (d-M1)×d matrix J̃ in `OPT_JACOBIMATRIX`.
  In this case banded Jacobian matrices are only supported for the
  case M1 + M2 == d. Then in this case J̃ is divided into d/M2 = 1+(M1/M2)
  blocks of the size M2×M2. All this blocks can be banded with a (common)
  lower bandwidth < M2. 
  
  The option `OPT_JACOBIBANDSTRUCT` is used to describe the banded structure
  of the Jacobian. It is eigher `nothing` if the Jacobian is full or a
  tuple `(l,u)` with the lower and upper bandwidth.
  
  The function for providing the Jacobian for ∂f/∂x can have the following
  forms:
  
       function (t,x,J)       -> nothing       (A)
       function (t,x,J1,…,JK) -> nothing       (B)
  
  The following table shows when `OPT_JACOBIMATRIX` has the form (A) 
  and when it has the form (B):
  
       ╔══════╤═════════════════════════════╤════════════════════════════════╗
       ║      │ JACOBIBANDSTRUCT == nothing │ JACOBIBANDSTRUCT == (l,u)      ║
       ╟──────┼─────────────────────────────┼────────────────────────────────╢
       ║ M1==0│ (A), J full, size: d×d      │ (A) J (l,u)-banded, size d×d   ║
       ╟──────┼─────────────────────────────┼────────────────────────────────╢
       ║ M1>0 │ (A), J full, size: (d-M1)×d │ (B) J1,…,JK (l,u)-banded       ║
       ║      │                             │     each with size M2×M2 and   ║
       ║      │                             │     K = 1 + M1/M2              ║
       ║      │                             │     M1 + M2 == d               ║
       ╚══════╧═════════════════════════════╧════════════════════════════════╝
  """
function help_specialstructure()
  return Docs.doc(help_specialstructure)
end

"""
  ## Internals

  1. What is the typical "call stack" for all this callbacks? see
     documentation of `DopriInternalCallInfos`,  `OdexInternalCallInfos` and
     `Radau5InternalCallInfos`.
  1. What closures (and how many) are generated to support the eval_sol_fcn?
     see `create_radau_eval_sol_fcn_closure`
  """
function help_internals()
  return Docs.doc(help_internals)
end

# vim:syn=julia:cc=79:fdm=indent:
