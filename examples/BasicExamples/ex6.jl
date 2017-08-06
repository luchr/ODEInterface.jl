using ODEInterface
@ODEInterface.import_huge

# We look at the boundary value problem
#
#                       sin⁴(x)
#    ε⋅y'(x) = sin²(x)-λ──────  ;   y(-π/2)=y(π/2)=1
#                         y
# 
# Here ε is given (e.g. ε=0.1) and λ is an unknown parameter: λ'=0
#
# We use colsys/colnew and a "homotopy": starting with ε=1.0, using this
# solution as start-guess for ε=0.5, then using this as start guess
# for ε=0.2 and then ε=0.1

a, b = -pi/2, pi/2
orders = [1, 1,]
ζ = [a, b]

found_pyplot = true
try
  using PyPlot
catch
  found_pyplot = false
end

global ε = nothing 
global ε_old = nothing
global sol_old = nothing

function rhs(x, z, f)
    s² = sin(x)^2
    f[1] = (s²-z[2]*s²*s²/z[1])/ε
    f[2] = 0.0
end

function Drhs(x, z, df)
    df[:]=0.0
    s⁴ = sin(x)^4
    df[1,1] = z[2]*s⁴/(z[1]^2)
    df[1,2] = -s⁴/z[1]
end

function bc(i, z, bc)
    bc[1] = z[1]-1.0
end

function Dbc(i, z, dbc)
    dbc[1] = 1.0
    dbc[2] = 0.0
end

function initial_guess(x, z, dmz)
    z[1] = 0.5
    z[2] = 1.0
    rhs(x, z, dmz)
end

opt = OptionsODE("example 6",
      OPT_BVPCLASS => 2, OPT_COLLOCATIONPTS => 7,
      OPT_RTOL => [1e-4, 1e-4], OPT_MAXSUBINTERVALS => 200)
xx = collect(linspace(a,b, 400));

sol = nothing
for ε = [1.0, 0.5, 0.2, 0.1]
  guess = (sol_old≠nothing) ? sol_old : initial_guess    
  sol, retcode, stats = colnew([a,b], orders, ζ, rhs, Drhs, bc, Dbc, guess ,opt);
  @printf("ε=%g, retcode=%i\n", ε, retcode)
  @assert retcode>0
  zz_new = evalSolution(sol, xx)
  if sol_old==nothing
      zz_old = copy(zz_new); zz_old[:,1] = 0.5; zz_old[:,2] = NaN
  else
      zz_old = evalSolution(sol_old, xx)
  end    
  if found_pyplot
    plot(xx, zz_old[:,1], xx, zz_new[:,1]);grid(true)
  end
  sol_old = sol; ε_old = ε
end

xx = collect(linspace(a, b, 9))
zz = evalSolution(sol, xx)
println("Solution at some points:")

display( [ xx' ; zz[:,1]' ] )
println()

# vim:syn=julia:cc=79:fdm=indent:
