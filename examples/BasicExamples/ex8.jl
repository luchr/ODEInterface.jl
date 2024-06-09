using ODEInterface
@ODEInterface.import_huge

# We look at the boundary value differential algebraic problem
#
#    x₁' = (ϵ+x₂-p₂(t))y+p₁'(t)
#    x₂' = p₂'(t)
#    x₃' = y
#    0   = (x₁-p₁(t))(y-eᵗ)
#
# with boundary conditions
#
#    x₁(0)=p₁(0), x₂(1)=p₂(1), x₃(0)=1
# 
# Here ε is given (e.g. ε=1)
#
# We use coldae to solve this BVDAE

a, b = 0.0, 1.0
orders = [1, 1, 1]
ζ = [0.0, 0.0, 1.0]

function rhs(x, z, y, f)
    e = 2.7
    f[1] = (1+z[2]-sin(x))*y[1] + cos(x)
    f[2] = cos(x)
    f[3] = y[1]
    f[4] = (z[1]-sin(x))*(y[1]-e^x)
end

function Drhs(x, z, y, df)
    df[:]=0.0
    df[1,1] = 0.0
    df[1,2] = y[1]
    df[1,3] = 0.0
    df[1,4] = 1+z[2]-sin(x)

    df[2,1] = 0.0
    df[2,2] = 0.0
    df[2,3] = 0.0
    df[2,4] = 0.0
    
    df[3,1] = 0.0
    df[3,2] = 0.0
    df[3,3] = 0.0
    df[3,4] = 1.0
    
    df[4,1] = y[1]-e^x
    df[4,2] = 0.0
    df[4,3] = 0.0
    df[4,4] = z[1]-sin(x)
end

function bc(i, z, bc)
    if i == 1
        bc[1] = z[1]
    elseif i == 2
        bc[1] = z[3] - 1.0
    elseif i == 3
        bc[1] = z[2] - sin(1.0)
    end
end

function Dbc(i, z, dbc)
    if i == 1
        dbc[1] = 1.0
        dbc[2] = 0.0
        dbc[3] = 0.0
    elseif i == 2
        dbc[1] = 0.0
        dbc[2] = 0.0
        dbc[3] = 1.0
    elseif i == 3
        dbc[1] = 0.0
        dbc[2] = 1.0
        dbc[3] = 0.0
    end
end

ny = 1
index = 1

opt = OptionsODE("example 8",
      OPT_BVPCLASS => 2, OPT_COLLOCATIONPTS => 7,
      OPT_RTOL => [1e-4, 1e-4, 1e-4], OPT_MAXSUBINTERVALS => 200)
xx = collect(linspace(a,b, 400));

guess = nothing  
sol, retcode, stats = coldae([a,b], orders, ny, index, ζ, rhs, Drhs, bc, Dbc, guess, opt);
@printf("ε=%g, retcode=%i\n", ε, retcode)
@assert retcode>0

xx = collect(linspace(a, b, 9))
zz = evalSolution(sol, xx)
