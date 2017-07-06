# Functions used for SLATEC Solvers like DDEABM

"""
  Some solvers have the conecpt of continuation calls (CC). This can be used
  to get some "kind" of dense output support. How does this work?

  Typically the solver "overshoots" beyond T for calculating the solution
  at T (then with the help of interpolation):

                                                  ↓ intermediate steps
       ──╫───┼────┼────┼─╫───┼─────             ──┼──
         t₀              T

  [Use `OPT_TSTOP` if the right-hand side is not defined for some t>T or
  if there are singularities.]

  If you do a CC and you want the solution at new Tₙ, then in the situtation

       ──╫───┼────┼────┼─╫──╫┼─────
         t₀              T  Tₙ

  where Tₙ is in the already computed interval, the solution at Tₙ can be
  computed simply by interpolation.
  If Tₙ is beyond the last computed interval then the solver continues with
  the process of numerical integration.

  So CC can be cheap and the solver is optimized for CC. CC are much more
  efficient than a "restart" at T.
  """
const SLATEC_continuation_call = nothing

"""
  ## License

  The solver ddeabm is part of the SLATEC Common Mathematical Library which
  is in the public domain. More informations can be found at

       http://www.netlib.org/slatec/guide

  """
const SLATEC_license = nothing

"""
  checks if output mode is supported and gives hint to
  `OPT_OUTPUTATTIMES` if the mode `OUTPUTFCN_DENSE` was requested.
  """
function check_slatec_output_mode(output_mode::OUTPUTFCN_MODE)
  if output_mode == OUTPUTFCN_DENSE
    throw(ArgumentErrorODE(
      string("Option 'OPT_OUTPUTMODE': Dense output not supported by this ",
        "solver. See description of 'OPT_OUTPUTMODE' and ",
        "'OPT_OUTPUTATTIMES' for the possibility to get the solution at ",
        "predefined time values"),
      :opt))
  end
  return nothing
end

"""
       function extractSlatecJacobiOpt{FInt<:FortranInt}(d::FInt,
               opt::AbstractOptionsODE)

  extracts jacobi options for some SLATEC solvers, like ddebdf.
  """
function extractSlatecJacobiOpt{FInt<:FortranInt}(d::FInt,
        opt::AbstractOptionsODE)
  OPT, jacobimatrix, jacobibandstruct = nothing, nothing, nothing
  try
    OPT = OPT_JACOBIMATRIX
    jacobimatrix = getOption(opt, OPT, nothing)
    @assert (jacobimatrix == nothing) || isa(jacobimatrix,Function)

    OPT = OPT_JACOBIBANDSTRUCT
    bs = getOption(opt, OPT, nothing)

    if bs ≠ nothing
      jacobibandstruct = ( convert(FInt,bs[1]), convert(FInt,bs[2]) )
      @assert 0 ≤ jacobibandstruct[1] ≤ d
      @assert 0 ≤ jacobibandstruct[2] ≤ d
    end
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid", :opt, e))
  end
  return (jacobimatrix, jacobibandstruct)
end

"""
  handles OPT_OUTPUTATTIMES for some SLATEC solvers, like ddeabm, ddebdf.

  returns always an Float64-Vector. Even if nothing was given by the user.
  If there was a vector given then it is copied and checked for
  monotonicity.
  """
function extractSlatecOutputAtTimes(t0, T, opt::AbstractOptionsODE)
  t_values = getOption(opt, OPT_OUTPUTATTIMES, nothing)
  if t_values == nothing 
        t_values = zeros(Float64, 0)
  end
  output_attimes_given = length(t_values)>0

  # convert *and* copy given OPT_OUTPUTATTIMES
  # "copy", because t_values will (typically) be changed by the caller
  try
    t_values = getVectorCheckLength(t_values, Float64, length(t_values), true)
  catch e
    throw(ArgumentErrorODE(
      string("I've a problem converting OPT_OUTPUTATTIMES to a ",
        "Float64-Vector"), :opt, e))
  end

  if output_attimes_given
    s = sign(T-t0)
    if any( x -> sign(x) != s, diff(t_values) )
      throw(ArgumentErrorODE(string("Because sign(T-t0)=sign($T-$t0)=$s ",
        "the vector in OPT_OUTPUTATTIMES has to be ",
        (T≥t0?"ascending":"descending"), ". But this is not the case."),:opt))
    end
    if s*(t0-t_values[1])≥0 || s*(t_values[end]-T)≥0
      throw(ArgumentErrorODE(string("All values of OPT_OUTPUTATTIMES ",
        "must be between t0=$t0 and T=$T"), :opt))
    end
  end

  return t_values :: Vector{Float64}
end

"""
       function unsafe_SLATEC1RHSCallback{CI<:ODEinternalCallInfos}(
                t_::Ptr{Float64}, x_::Ptr{Float64},
                f_::Ptr{Float64}, rpar_::Ptr{Float64}, cbi::CI)
  
  This is the right-hand side given as callback to SLATEC solvers,
  e.g. ddeabm.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  """
function unsafe_SLATEC1RHSCallback{CI<:ODEinternalCallInfos}(
        t_::Ptr{Float64}, x_::Ptr{Float64},
        f_::Ptr{Float64}, rpar_::Ptr{Float64}, cbi::CI)
  n = cbi.N
  t = unsafe_load(t_)
  x = unsafe_wrap(Array, x_, (n,), false)
  f = unsafe_wrap(Array, f_, (n,), false)

  cbi.rhs_count += 1
  hw1rhs(n, t, x, f, cbi)
  return nothing
end

"""
        function unsafe_SLATEC1RHSCallback_c{FInt,CI}(cbi::CI, fint_flag::FInt)
          -> C-callable function pointer

  This method generates a Pointer to C-callable instructions.
  The two method type parameters `FInt` and `CI` are important:
  `FInt` is the used Fortran integer type and `CI` is the used 
  `ODEinternalCallInfos` *SubType*.
  Because `unsafe_SLATEC1RHSCallback` is a parameterized method,
  special variants are compiled, if `FInt` or `CI` changes.
  If `CI` itself is a parameterized type (depending on all the
  user-given Julia-functions like right-hand side, etc.) then
  calls to such Julia-functions can be resolved at compile-time
  (instead of dynamic calls during run-time).
  """
function unsafe_SLATEC1RHSCallback_c{FInt,CI}(cbi::CI, fint_flag::FInt) 
  return cfunction(unsafe_SLATEC1RHSCallback, Void, (Ptr{Float64},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{CI}))
end


function unsafe_SLATEC1JacCallback{FInt<:FortranInt,
        CI<:ODEinternalCallInfos}(t_::Ptr{Float64}, x_::Ptr{Float64}, 
        dfx_::Ptr{Float64}, dfx_rows_::Ptr{FInt}, rpar_::Ptr{Float64},
        cbi::CI)
  n = cbi.N
  t = unsafe_load(t_)
  x = unsafe_wrap(Array, x_, (n,), false)
  dfx_rows = unsafe_load(dfx_rows_)
  M = unsafe_wrap(Array, dfx_, (dfx_rows,n,), false)

  cbi.jac_count += 1

  lprefix = cbi.jac_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_jac = l & LOG_JAC>0

  l_jac && println(lio,lprefix,"called with dfx_rows=",dfx_rows," t=",t)
  jac = cbi.jacobimatrix
  jb = cbi.jacobibandstruct
  if jb == nothing
    @assert dfx_rows == n
    J = M
    jac(t, x, J)
  else
    @assert dfx_rows ≥ 1+2*jb[1]+jb[2]
    banded_memory = cbi.submatrix(M, 1+jb[1]:1+2*jb[1]+jb[2], :)
    J = BandedMatrix{Float64}(n, n, jb[1], jb[2], banded_memory)
    jac(t, x, J)
  end
  l_jac && println(lio,lprefix,"dfx=",M)
  return nothing
end

"""
       function unsafe_SLATEC1JacCallback_c{FInt, CI}
                (cbi::CI, fint_flag::FInt)
          -> C-callable function pointer
  """
function unsafe_SLATEC1JacCallback_c{FInt, CI}(cbi::CI, fint_flag::FInt)
  return cfunction(unsafe_SLATEC1JacCallback, Void, (Ptr{Float64},
    Ptr{Float64}, Ptr{Float64}, Ptr{FInt}, Ptr{Float64}, Ref{CI}))
end

# vim:syn=julia:cc=79:fdm=indent:
