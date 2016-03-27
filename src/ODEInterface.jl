__precompile__(false)

"""
  # ODEInterface
  
  This julia module provides an interface to solvers for 
  ordinary differential equations (ODEs) written in Fortran
  for solving initial value problems of the form
  
      x' = rhs(t,x),      x(t₀) = x₀
  
  or (for solvers supporting a "mass matrix" M)
  
      M⋅x' = rhs(t,x),    x(t₀) = x₀.
  
  ## What does "Interface" mean?
  
  This julia module does *not* contain code for solving initial value
  problems, but this module does contain code for interacting with
  compiled Fortran-solvers. That's the reason, why this module is not called
  ODESuite.
  
  ## What solvers are currently supported?
  
  Currently the following Fortran-solvers, written by
  Prof. E. Hairer and Prof. G. Wanner, are supported:
  
  * dopri5: explicit Runge-Kutta method of order 5(4) due to Dormand & Prince
  * dop853: explicit Runge-Kutta method of order 8(5,3) due to Dormand & Prince
  * odex: GBS extrapolation-algorithm based on the explicit midpoint rule
  * radau5: implicit Runge-Kutta method (Radau IIA) of order 5
  * radau: implicit Runge-Kutta method (Radau IIA) of variable order 
  between 5 and 13
  
  see [Software page of Prof. Hairer](http://www.unige.ch/~hairer/software.html).
  
  The following features of this solvers are supported by this ODEInterface:
  
  * providing an output function (e.g. for dense output or for event location)
  to the solvers
  * providing mass- and jacobi-matrices for the solvers (with support for
  banded matrices)
  * all the solvers' parameters for fine-tuning them
  * support for problems with "special structure", see `help_specialstructure`
  
  ## What are the requirements for this module
  
  In order to use this module, you have to *compile* the supported
  Fortran solvers and provide a shared library for each solver. Just call
  `ODEInterface.help_solversupport` for further informations (help topics)
  on how to compile the solvers and how to create shared libraries.
  
  ## Further help
  
  see `ODEInterface.help_overview` for an overview of some help topics. 
  
  ## Contacting the author of this module
  
  The author of this julia module is 
  
       Dr. Christian Ludwig
       email: ludwig@ma.tum.de
         (Faculty of Mathematics, Technische Universität München)
  
  """
module ODEInterface

const VERSION = "0.0.8-dev"

using Base

include("./Error.jl")
include("./Options.jl")
include("./DLSolvers.jl")
include("./Base.jl")
include("./SoloutEval.jl")
include("./Banded.jl")

"""macro for importing the *huge* set of symbols."""
macro import_huge()
  quote
    @ODEInterface.import_RHScallmode
    @ODEInterface.import_bandedmatrixfuncs
    @ODEInterface.import_LOG
    @ODEInterface.import_dynamicload
    @ODEInterface.import_dopri5
    @ODEInterface.import_DLdopri5
    @ODEInterface.import_dop853
    @ODEInterface.import_DLdop853
    @ODEInterface.import_odex
    @ODEInterface.import_exceptions
    @ODEInterface.import_options
    @ODEInterface.import_outputfcn
    @ODEInterface.import_OPTcommon
    @ODEInterface.import_odecall
    @ODEInterface.import_radau5
    @ODEInterface.import_DLradau5
    @ODEInterface.import_radau
    @ODEInterface.import_DLradau
    @ODEInterface.import_seulex
    @ODEInterface.import_DLseulex
    @ODEInterface.import_bvpsol
    @ODEInterface.import_DLbvpsol
  end
end

"""macro for importing the *normal* set of symbols."""
macro import_normal()
  quote
    @ODEInterface.import_RHScallmode
    @ODEInterface.import_bandedmatrix
    @ODEInterface.import_dynamicload
    @ODEInterface.import_dopri5
    @ODEInterface.import_dop853
    @ODEInterface.import_odex
    @ODEInterface.import_radau5
    @ODEInterface.import_radau
    @ODEInterface.import_seulex
    @ODEInterface.import_bvpsol
    @ODEInterface.import_options
    @ODEInterface.import_OPTcommon
  end
end

"""
  Type describing a "variant" of a solver.
  
  What is a variant of a solver? Some solvers support more than one
  forms of dynamic libraries, e.g. with 32bit integers and with 64bit
  integers. The purpose of this type is to have enough fields for 
  describing such an variant.
  """
immutable SolverVariant
  name          :: AbstractString   # name of the variant
  description   :: AbstractString   # a short description of the variant
  libname       :: AbstractString   # the name of the dynamic library
  methods       :: Tuple            # tuples with strings of methods to load
end

"""
  Type describing a solver.
  """
immutable SolverInfo
  name          :: AbstractString   # name of the Solver
  description   :: AbstractString   # a short description of the solver
  supported_opts:: Tuple            # Supported OPT_... Symbols
  variants      :: Tuple            # Tuple with SolverVariant(s)
  help_compile  :: Function         # Function with compile info
  help_license  :: Function         # Function with license info
end

"""
  In this array every (supported) solver appends its SolverInfo-record.
  """
const solverInfo = Vector{SolverInfo}() 

"""
  Type encapsulating all required data for ODE-Solver-Callbacks.
  
  For further explanation see, `GlobalCallInfoDict`
  """
abstract ODEinternalCallInfos

"""
  Global `Dict` for saving ODE-Server-Callback-Infos stored with `UInt64` keys.
  
  All the ODE-solvers of this ODEInterface use callback functions:
  In order to evalute the right-hand-side (RHS) the Fortran-Codes want to have
  a code-pointer (i.e. address) for a compiled C or Fortran-Function to
  call. There are other callback mechanisms, too: e.g. calling a
  output function after every successfull integration step.
  
  This interface wants to enable the user to supply julia functions as
  callback functions, e.g. a julia anonymous function `(t,x)->x` for the RHS.
  But a julia function is much more than a simple code-pointer. One can try to
  'convert' a julia function to a code-pointer/address with `cfunction`, but
  this is (by principle) not always possible (e.g. if the julia function has
  some closure information attached to it), see
  "Passing Callback Functions to C": 
  http://julialang.org/blog/2013/05/callback/ 
  for a more detailed explanation.
  
  In order to support even anonymous julia functions (and closures) this
  interface uses the following idea:
  All necessary informations (e.g. a julia function reference to the RHS) is
  stored globally and a global function (e.g. `unsafe_dopri5RHSCallback`)
  is used (as `unsafe_dopri5RHSCallback_c`) as code-pointer for the
  Fotran routines.
  
  How are the necessary informations stored globally?
  Using global variables can have the negative side-effect of being not
  re-entrant safe. Imagine a julia-RHS that in turn solves a ODE for
  calculating the RHS (i.e. imagine nested usages of a solver).
  In order to make this re-entrant safe the following idea is used:
  Every call of an ODE-solver gets an unique `UInt64`-id (the call-id). In
  this `GlobalCallInfoDict` all necessary informations for this call
  are stored at this call-id-key. Nested calls gat different call-ids and
  are stored at different keys in this `Dict`.
  This call-id is passed (as pass-through paramter) to the Fortran solvers.
  The solvers pass this call-id to every call of the RHS. There (e.g. in
  `unsafe_dopri5RHSCallback`) all in informations needed can be found in the
  `GlobalCallInfoDict`.
  At the end, when the Fortran-solver has finished, the key for the
  call-id is deleted.
  """
const GlobalCallInfoDict = Dict{UInt64,ODEinternalCallInfos}()

"""
  Ancestor for all types storing arguments for ODE-(C-/Fortran-)solvers.
  """
abstract AbstractArgumentsODESolver{FInt} <: Any

# Common options
"""macro for importing common OPT options."""
macro import_OPTcommon()
  :(
    using ODEInterface:   OPT_LOGIO, OPT_LOGLEVEL, OPT_RHS_CALLMODE,
                          OPT_RTOL, OPT_ATOL,
                          OPT_MAXSTEPS, OPT_EPS, OPT_OUTPUTFCN,
                          OPT_OUTPUTMODE, OPT_STEST, OPT_RHO, OPT_SSMINSEL,
                          OPT_SSMAXSEL, OPT_SSBETA, OPT_MAXSS, OPT_INITIALSS,
                          OPT_MAXEXCOLUMN, OPT_MAXSTABCHECKS, 
                          OPT_MAXSTABCHECKLINE, OPT_INTERPOLDEGREE,
                          OPT_ORDERDECFRAC, OPT_ORDERINCFRAC, 
                          OPT_STEPSIZESEQUENCE,
                          OPT_SSREDUCTION, OPT_SSSELECTPAR1, OPT_SSSELECTPAR2,
                          OPT_RHO2, OPT_DENSEOUTPUTWOEE, OPT_TRANSJTOH,
                          OPT_MAXNEWTONITER, OPT_NEWTONSTARTZERO,
                          OPT_DIMOFIND1VAR, OPT_DIMOFIND2VAR, OPT_DIMOFIND3VAR,
                          OPT_STEPSIZESTRATEGY, OPT_M1, OPT_M2,
                          OPT_JACRECOMPFACTOR, OPT_NEWTONSTOPCRIT, 
                          OPT_FREEZESSLEFT, OPT_FREEZESSRIGHT, OPT_MASSMATRIX,
                          OPT_JACOBIMATRIX, OPT_JACOBIBANDSTRUCT,
                          OPT_MAXSTAGES, OPT_MINSTAGES, OPT_INITSTAGES,
                          OPT_ORDERDECFACTOR, OPT_ORDERINCFRAC,
                          OPT_ORDERDECSTEPFAC1, OPT_ORDERDECSTEPFAC2,
                          OPT_RHSAUTONOMOUS, OPT_LAMBDADENSE,
                          OPT_WORKFORRHS, OPT_WORKFORJAC, OPT_WORKFORDEC,
                          OPT_WORKFORSOL, OPT_BVPCLASS, OPT_SOLMETHOD,
                          OPT_IVPOPT
  )
end

"""macro for importing RHS_CALL_MODE."""
macro import_RHScallmode()
  :(
    using ODEInterface: RHS_CALL_MODE, RHS_CALL_RETURNS_ARRAY,
                        RHS_CALL_INSITU
  )
end

const OPT_LOGIO            = "logio"
const OPT_LOGLEVEL         = "loglevel"
const OPT_RHS_CALLMODE     = "RightHandSideCallMode"

const OPT_RTOL             = "RelTol"
const OPT_ATOL             = "AbsTol"
const OPT_MAXSTEPS         = "MaxNumberOfSteps"
const OPT_EPS              = "eps"

const OPT_OUTPUTFCN        = "OutputFcn"
const OPT_OUTPUTMODE       = "OutputFcnMode"

const OPT_STEST            = "StiffTestAfterStep"
const OPT_RHO              = "rho"
const OPT_SSMINSEL         = "StepSizeMinSelection"
const OPT_SSMAXSEL         = "StepSizeMaxSelection"
const OPT_SSBETA           = "StepSizeBeta"
const OPT_MAXSS            = "MaxStep"
const OPT_INITIALSS        = "InitialStep"


const OPT_MAXEXCOLUMN      = "MaxExtrapolationColumn"
const OPT_MAXSTABCHECKS    = "MaxNumberOfStabilityChecks"
const OPT_MAXSTABCHECKLINE = "MaxLineForStabilityCheck"
const OPT_INTERPOLDEGREE   = "DegreeOfInterpolation"
const OPT_ORDERDECFRAC     = "OrderDecreaseFraction"
const OPT_ORDERINCFRAC     = "OrderIncreaseFraction"
const OPT_STEPSIZESEQUENCE = "StepSizeSequence"
const OPT_SSREDUCTION      = "StepSizeReduction"
const OPT_SSSELECTPAR1     = "StepSizeSelectionParam1"
const OPT_SSSELECTPAR2     = "StepSizeSelectionParam2"
const OPT_RHO2             = "rho2"
const OPT_DENSEOUTPUTWOEE  = "DeactivateErrorEstInDenseOutput"

const OPT_TRANSJTOH        = "TransfromJACtoHess"
const OPT_MAXNEWTONITER    = "MaxNewtonIterations"
const OPT_NEWTONSTARTZERO  = "StartNewtonWithZeros"
const OPT_DIMOFIND1VAR     = "DimensionOfIndex1Vars"
const OPT_DIMOFIND2VAR     = "DimensionOfIndex2Vars"
const OPT_DIMOFIND3VAR     = "DimensionOfIndex3Vars"
const OPT_STEPSIZESTRATEGY = "StepSizeStrategy"
const OPT_M1               = "M1"
const OPT_M2               = "M2"
const OPT_JACRECOMPFACTOR  = "RecomputeJACFactor"
const OPT_NEWTONSTOPCRIT   = "NewtonStopCriterion"
const OPT_FREEZESSLEFT     = "FreezeStepSizeLeftBound"
const OPT_FREEZESSRIGHT    = "FreezeStepSizeRightBound"
const OPT_MASSMATRIX       = "MassMatrix"
const OPT_JACOBIMATRIX     = "JacobiMatrix"
const OPT_JACOBIBANDSTRUCT = "JacobiBandStructure"

const OPT_MAXSTAGES        = "MaximalNumberOfStages"
const OPT_MINSTAGES        = "MinimalNumberOfStages"
const OPT_INITSTAGES       = "InitialNumberOfStages"
const OPT_ORDERINCFACTOR   = "OrderIncreaseFactor"
const OPT_ORDERDECFACTOR   = "OrderDecreaseFactor"
const OPT_ORDERDECSTEPFAC1 = "OrderDecreaseStepFactor1"
const OPT_ORDERDECSTEPFAC2 = "OrderDecreaseStepFactor2"

const OPT_RHSAUTONOMOUS    = "AutonomousRHS"
const OPT_LAMBDADENSE      = "LambdaForDenseOutput"
const OPT_WORKFORRHS       = "WorkForRightHandSide"
const OPT_WORKFORJAC       = "WorkForJacobimatrix"
const OPT_WORKFORDEC       = "WorkForLuDecomposition"
const OPT_WORKFORSOL       = "WorkForSubstitution"

const OPT_BVPCLASS         = "BoundaryValueProblemClass"
const OPT_SOLMETHOD        = "SolutionMethod"
const OPT_IVPOPT           = "OptionsForIVPsolver"


@enum(RHS_CALL_MODE,
      RHS_CALL_RETURNS_ARRAY, RHS_CALL_INSITU)

@doc """
  The right-hand side has to return `x'` as Array.
  
  The right-hand side must be a function of the form
  
       funtion (t,x) -> dx
  
  `dx` is a `Vector{Float64}` with the same length as `x`.
  """ -> RHS_CALL_RETURNS_ARRAY

@doc """
  The right-hand side has to return `nothing`. It gets an additional Array
  where it has to save the the values of `x'`.

  The right-hand side must be a function of the form
  
       funtion (t,x,dx) -> nothing
  
  `dx` is a `Vector{Float64}` with the same length as `x`. In `dx` the
  function has to fill in the values of `x'`.
  """ -> RHS_CALL_INSITU

"""
       function extractTOLs(d::Integer, opt::AbstractOptionsODE) 
              -> (scalarFlag,rtol,atol)
  
  extract `OPT_RTOL` and `OPT_ATOL` and convert to `Vector{Float64}`.
  
  Supports scalar `OPT_RTOL` and `OPT_ATOL` and converts them to a
  `Vector{Float64}` of length 1.
  """
function extractTOLs(d::Integer,opt::AbstractOptionsODE)
  rtol = getOption(opt,OPT_RTOL,1e-3)
  atol = getOption(opt,OPT_ATOL,1e-6)
  scalarFlag=false

  if isscalar(rtol) && isscalar(atol)
    scalarFlag=true
    try
      rtol = [ convert(Float64,rtol) ]
      atol = [ convert(Float64,atol) ]
    catch e
      throw(ArgumentErrorODE(
        string("OPT_RTOL and OPT_ATOL are scalars. ",
               "But I've a problem converting rtol=$rtol,",
               " and/or atol=$atol to Float64"),:opt,e))
    end
  else
    try
      rtol = getVectorCheckLength(rtol,Float64,d)
      atol = getVectorCheckLength(atol,Float64,d)
    catch e
      throw(ArgumentErrorODE(
        string("OPT_RTOL and/or OPT_ATOL is not a real scalar. ",
               "But I've also a problem converting them to Float64 ",
               "vectors of length $d"),:opt,e))
    end
  end

  return (scalarFlag,rtol,atol)
end

"""
       function extractLogOptions(opt::AbstractOptionsODE) -> (lio, lev)
  
  Extract options for logging.
  
  throws ArgumentErrorODE if logio is not an IO or
  if loglevel is not convertable to UInt64.
  """
function extractLogOptions(opt::AbstractOptionsODE)
  lio=getOption(opt,OPT_LOGIO,STDERR)
  if !isa(lio,IO)
    throw(ArgumentErrorODE("option '$OPT_LOGIO' was not an Base.IO",:opt))
  end
  lev=nothing
  try
    lev=convert(UInt64,getOption(opt,OPT_LOGLEVEL,0))
  catch e
    throw(ArgumentErrorODE(
      "option '$OPT_LOGLEVEL' cannot be converted to UInt64",:opt,e))
  end
  return (lio,lev)
end

"""
       function extractOutputFcn(opt::AbstractOptionsODE) 
              -> (output_mode, output_fcn)
  """
function extractOutputFcn(opt::AbstractOptionsODE)
  OPT = nothing
  output_mode = OUTPUTFCN_NEVER; output_fcn = output_fcn_donothing
  try 
    OPT = OPT_OUTPUTMODE; output_mode = getOption(opt,OPT,OUTPUTFCN_NEVER)
    @assert isa(output_mode,OUTPUTFCN_MODE)
    
    if output_mode ≠ OUTPUTFCN_NEVER
      OPT = OPT_OUTPUTFCN; 
      output_fcn = getOption(opt,OPT,nothing)
      @assert isa(output_fcn,Function)
    end
  catch e
    throw(ArgumentErrorODE("Option '$OPT' invalid",:opt,e))
  end
  
  return (output_mode,output_fcn)
end

"""
       function solver_init(solver_name::AbstractString, 
                            opt::AbstractOptionsODE, cid=nothing)
          ->  (lio,l,l_g,l_solver,lprefix,cid,cid_str)
  """
function solver_init(solver_name::AbstractString, opt::AbstractOptionsODE,
                     cid=nothing)
  (lio,l) = extractLogOptions(opt);
  (l_g,l_solver)=( l & LOG_GENERAL>0, l & LOG_SOLVERARGS>0 )
  if cid==nothing
    cid = uniqueToken()
  end
  cid_str = int2logstr(cid[1])
  lprefix = string(cid_str,solver_name,": ")
  return (lio,l,l_g,l_solver,lprefix,cid,cid_str)
end

"""
       function solver_start(solver_name::AbstractString, rhs::Function, 
                   t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE)
          ->  (lio,l,l_g,l_solver,lprefix,cid,cid_str)
  
  initialization for a (typical) solver call/start.
  """
function solver_start(solver_name::AbstractString, rhs::Function, 
            t0::Real, T::Real, x0::Vector, opt::AbstractOptionsODE)
  
  (lio,l,l_g,l_solver,lprefix,cid,cid_str) = solver_init(solver_name,opt)

  if l_g
    println(lio,lprefix,
      "called with rhs=",rhs," t0=",t0," T=",T," x0=",x0," and opt")
    show(lio,opt); println(lio)
  end
  
  return (lio,l,l_g,l_solver,lprefix,cid,cid_str)
end

"""
       function solver_extract_rhsMode(opt::AbstractOptionsODE)
                   -> rhs_mode
  
  extracts `rhs_mode` from `opt`.
  """
function solver_extract_rhsMode(opt::AbstractOptionsODE)
  try
    rhs_mode = getOption(opt,OPT_RHS_CALLMODE,RHS_CALL_RETURNS_ARRAY)
    @assert rhs_mode ∈ (RHS_CALL_RETURNS_ARRAY,RHS_CALL_INSITU)
    return rhs_mode
  catch e
    throw(ArgumentErrorODE("Option 'OPT_RHS_CALLMODE' is invalid",:opt,e))
  end
end

"""
       function solver_extract_commonOpt{FInt}(t0::Real, T::Real, x0::Vector, 
             opt::AbstractOptionsODE, args::AbstractArgumentsODESolver{FInt})
         -> (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn)

  get d, fill args.N, args.x, args.t, args.tEnd, args.RTOL, args.ATOL
  """
function solver_extract_commonOpt{FInt}(t0::Real, T::Real, x0::Vector, 
             opt::AbstractOptionsODE, args::AbstractArgumentsODESolver{FInt})
  
  d = FInt(0)
  try
    d=convert(FInt,length(x0))
  catch e
    throw(ArgumentErrorODE("Cannot convert length(x0) to $FInt",:x0,e))
  end
  0==d && throw(ArgumentErrorODE("x0 was empty vector",:x0))
  args.N = [d]
  try
    args.x = getVectorCheckLength(x0,Float64,d)
  catch e
    throw(ArgumentErrorODE("Cannot convert x0 to Vector{Float64}",:x0,e))
  end

  try
    args.t = [ convert(Float64,t0) ]
  catch e
    throw(ArgumentErrorODE("Cannot convert t0=$t0 to Float64",:t0,e))
  end
  try
    args.tEnd = [ convert(Float64,T) ]
  catch e
    throw(ArgumentErrorODE("Cannot convert T=$T to Float64",:T,e))
  end

  (scalarFlag,args.RTOL,args.ATOL) = extractTOLs(d,opt)

  rhs_mode = solver_extract_rhsMode(opt)

  (output_mode,output_fcn) = extractOutputFcn(opt)
  if output_mode == OUTPUTFCN_NEVER
    nrdense = FInt(0);
  else
    nrdense = d; 
  end

  return (d,nrdense,scalarFlag,rhs_mode,output_mode,output_fcn)
end

include("./HWcommon.jl")

# Solvers:
include("./Dopri.jl")
include("./Dopri5.jl")
include("./Dop853.jl")
include("./Odex.jl")
include("./Radau.jl")
include("./Seulex.jl")
include("./Bvpsol.jl")

include("./Call.jl")
include("./Help.jl")

end


# vim:syn=julia:cc=79:fdm=indent:
