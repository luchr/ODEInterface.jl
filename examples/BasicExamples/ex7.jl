using ODEInterface
@ODEInterface.import_huge

# We look at (again) the boundary value problem (see, ex6.jl)
#
#                       sin⁴(x)
#    ε⋅y'(x) = sin²(x)-λ──────  ;   y(-π/2)=y(π/2)=1
#                         y
# 
# Here ε is given (e.g. ε=0.1) and λ is an unknown parameter.
#
# We use bvpm2 and a "homotopy": starting with ε=0.1, using this
# solution as start-guess for smaller ε.

a, b = -pi/2, pi/2

found_pyplot = true
try
  using PyPlot
catch
  found_pyplot = false
end

global ε = nothing 
global ε_old = nothing
global sol_old = nothing

function rhs(x, y, p, f)
  f[1] = ( sin(x)^2 - p[1]*sin(x)^4/y[1] ) / ε
end

function Drhs(x, y, p, dfdy, dfdp)
  dfdy[1,1] = ( p[1]*sin(x)^4/y[1]^2 ) / ε
  dfdp[1,1] = ( -sin(x)^4/y[1] ) / ε
end

function bc(ya, yb, p, bca, bcb)
  bca[1] = ya[1] - 1.0
  bcb[1] = yb[1] - 1.0
end

function Dbc(ya, yb, dya, dyb, p, dpa, dpb)
  dya[1,1] = 1.0
  dyb[1,1] = 1.0
  dpa[1,1] = 0.0
  dpb[1,1] = 0.0
end

opt = OptionsODE("example 7",
          OPT_RTOL => 1e-6,
          OPT_METHODCHOICE => 4,)
xx = collect(linspace(a,b, 400));

initial_guess = Bvpm2()
bvpm2_init(initial_guess, 1, 1, collect(linspace(a, b, 20)), [0.5,], [1.0,])

sol = nothing
for ε = [0.1, 0.05, 0.02, 0.01]
  guess = (sol_old≠nothing) ? sol_old : initial_guess
  sol, retcode, stats = bvpm2_solve(guess, rhs, bc, opt, Drhs=Drhs, Dbc=Dbc)
  @printf("ε=%g, retcode=%i\n", ε, retcode)
  @assert retcode ≥ 0
  zz_new = evalSolution(sol, xx)
  if sol_old == nothing
      zz_old = copy(zz_new); zz_old[1,:] = 0.5;
  else
      zz_old = evalSolution(sol_old, xx)
      bvpm2_destroy(sol_old)
  end    
  if found_pyplot
    plot(xx, zz_old[1,:], xx, zz_new[1,:]);grid(true)
  end
  sol_old = sol; ε_old = ε
end

xx = collect(linspace(a, b, 9))
zz = evalSolution(sol, xx)
println("Solution at some points:")

display( [ xx' ; zz[1,:]' ] )
println()

bvpm2_destroy(initial_guess)
bvpm2_destroy(sol)

# vim:syn=julia:cc=79:fdm=indent:
