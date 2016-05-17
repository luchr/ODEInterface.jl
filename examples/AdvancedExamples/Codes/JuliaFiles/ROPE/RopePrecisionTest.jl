
# Check if all the packages are installed or not
cond = "Gadfly" in keys(Pkg.installed()) &&
"Colors" in keys(Pkg.installed()) &&
"ODEInterface" in keys(Pkg.installed());
@assert cond "Please check if the following package(s) are installed:\
    Gadfly\
    Colors\
    ODEInterface"

# Load all the required packages
using Gadfly
using Colors
using ODEInterface
@ODEInterface.import_huge
loadODESolvers();

######################## Function for saving plots #################################
# Input:
# fileName = Name of the file where the plot is to be stored
#            (with or without extension)
# f_e = Array containing function evaluations as columns for each solver
# err = Array containing erros as columns for each solver
# solverNames = Array containing the names of solvers used in respective order
# plotSize = size of the plot to be created
# 
# Values have been tuned for a graph similar to the one in 
# Solving Ordinary Differential Equations I by
# Hairer, Ernst, Nørsett, Syvert P., Wanner, Gerhard
# page: 252
###################################################################################
function savePlotPNG(fileName,f_e,err,solverNames,
    plotSize=[30cm,30cm])
    
    numOfLayers = length(solverNames);
    
    if !contains(fileName,".")
        fileName = string(fileName,".png");
    end
    
    plotColorsHex = ["#4D4D4D","#5DA5DA","#FAA43A","#60BD68",
        "#F17CB0","#B2912F","#B276B2", "#DECF3F","#F15854"];
    plotColors = [parse(Colorant,c) for c in plotColorsHex];
    
    majorFontSize = 24pt;
    minorFontSize = 20pt;
    pointSize = 5pt;
    
    myplot = plot(Scale.x_log10,Scale.y_log10,
        Coord.cartesian(xflip=true),
        Guide.manual_color_key("Legend",solverNames,plotColorsHex[1:numOfLayers]),
        Guide.xlabel("error"),Guide.ylabel("#Function Evaluations"),
        Guide.xticks(ticks=[0:-3:-13;]),Guide.yticks(ticks=[3:1:5.1;]),
        Theme(major_label_font_size=majorFontSize,panel_stroke=colorant"black",
        minor_label_font_size=minorFontSize,key_title_font_size=majorFontSize,
        key_label_font_size=minorFontSize,key_position=:top,key_max_columns=1));
    
    for i = 1:numOfLayers
        push!(myplot,layer(x=err[:,i],y=f_e[:,i],Geom.point,Geom.path,
        Theme(default_color=plotColors[i],default_point_size=pointSize)));
    end
    
    draw(PNG(fileName,plotSize[1],plotSize[2]),myplot)
    return nothing
end

# Number of subdivisions of the rope
global n = 40;

# Define the system of ODEs
function rope(t,x,dx)
    n2 = n*n; # n^2
    n3by4 = convert(Int64,3*n/4); # 3*n/4
    
    # Force in x-direction
    Fx = 0.4;
    # Force in y-direction
    Fy = cosh(4*t-2.5)^(-4);
    
    # Compute required matrices
    c = -cos(x[1:n-1]-x[2:n]);
    cDiag = [1;2*ones(n-2);3];
    C = spdiagm((c,cDiag,c),(-1,0,1));
    
    d = -sin(x[1:n-1]-x[2:n]);
    D = spdiagm((-d,d),(-1,1));
    
    # Compute the inhomogeneous term
    v = -(n2+n/2-n*[1:n;]).*sin(x[1:n])-n2*sin(x[1:n])*Fx;
    v[1:n3by4] = v[1:n3by4] + n2*cos(x[1:n3by4])*Fy; 
    
    w = D*v+x[n+1:2*n].^2;
    u = C\w;
    
    # Write down the system
    dx[1:n] = x[n+1:2*n];
    dx[n+1:2*n] = C*v + D*u;
    
    return nothing
end

# Initial Conditions
t0 = 0.0; T = 3.723; x0=zeros(2*n);

# Compute the "reference solution"
opt = OptionsODE(OPT_EPS => 1.11e-16,OPT_RHS_CALLMODE => RHS_CALL_INSITU,
OPT_RTOL => 1e-16,OPT_ATOL=>1e-16);
(t,x_ref,retcode,stats) = dop853(rope,t0, T, x0, opt);

if retcode != 1
    println("Reference solution failed")
else
    # Initialization for the loop
    # f_e = function evaluations
    f_e = zeros(Int32,89,3);
    # err = error for last step using infinity norm
    err = zeros(Float64,89,3);

    # solverNames = names of the solvers used for the plot
    solverNames = ["DOPRI5","DOP853","ODEX"];

    # Compute all the solutions
    for i=0:88
        
        # Set up the tolerance
        Tol = 10^(-3-i/8);
        
        # Set up solver options
        opt = OptionsODE(OPT_EPS => 1.11e-16,OPT_RHS_CALLMODE => RHS_CALL_INSITU,
        OPT_RTOL => Tol,OPT_ATOL => Tol);

        # Solve using DOPRI5
        (t,x,retcode,stats) = dopri5(rope,t0, T, x0, opt);
        # Check if solver was successful
        if retcode != 1
            printFlag = false;
            break;
        end
        f_e[i+1,1] = stats.vals[13];
        err[i+1,1] = norm(x_accurate[1:n] - x[1:n],Inf);

        # Solve using DOP853
        (t,x,retcode,stats) = dop853(rope,t0, T, x0, opt);
        # Check if solver was successful
        if retcode != 1
            printFlag = false;
            break;
        end
        f_e[i+1,2] = stats.vals[13];
        err[i+1,2] = norm(x_accurate[1:n] - x[1:n],Inf);

        # Solve using ODEX
        (t,x,retcode,stats) = odex(rope,t0, T, x0, opt);
        # Check if solver was successful
        if retcode != 1
            printFlag = false;
            break;
        end
        f_e[i+1,3] = stats.vals[13];
        err[i+1,3] = norm(x_accurate[1:n] - x[1:n],Inf);
    end

    # Save the plot in PNG format
    if printFlag
        savePlotPNG("RopeConvTest",f_e,err,solverNames);
    else
        println("Cannot generate plot due to solver failure")
    end
end
