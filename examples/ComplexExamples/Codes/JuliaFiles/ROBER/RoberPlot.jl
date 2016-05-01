
# Load all the required packages
using ODEInterface
using ForwardDiff
using Gadfly
using Colors
@ODEInterface.import_huge
loadODESolvers();

# Define the right-hand function for Automatic Differentiation
function roberAD(x)
    return [-0.04*x[1]+1e4*x[2]*x[3],
        0.04*x[1]-1e4*x[2]*x[3]-3e7*(x[2])^2,
        3*10^7*(x[2])^2]
end

# Define the system for the solver
function rober(t,x,dx)
    dx[1] = -0.04*x[1]+1e4*x[2]*x[3];
    dx[2] = 0.04*x[1]-1e4*x[2]*x[3]-3e7*(x[2])^2;
    dx[3] = 3e7*(x[2])^2;
    return nothing
end

# Automatic Differentiation for a more general problem
function getJacobian(t,x,J)
    J[:,:] = ForwardDiff.jacobian(roberAD,x);
    return nothing
end

# Flag to check whether plot is to be generated and saved or not
# Also checks if all solvers are successful
printFlag = true;

# Initial conditions
t0 = 0.0; T = 10.^[0.0:11.0;]; x0=[1.0,0.0,0.0];

# Get "reference solution"
# TolMin < 1e-14 gives error "TOLERANCES ARE TOO SMALL" 
Tol = 1e-14;
opt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
    OPT_RTOL => Tol, OPT_ATOL=>Tol*1e-6,OPT_EPS => 1.11e-16,
    OPT_JACOBIMATRIX => getJacobian);

(t,x,retcode,stats) = odecall(seulex,rober,[t0, T[end]], x0, opt);
t = [t0;t];
x = [x0';x];

plotColorsHex = ["#4D4D4D","#5DA5DA","#FAA43A"];
plotColors = [parse(Colorant,c) for c in plotColorsHex];

p = plot(layer(x=t,y=x[:,1],Geom.path,Theme(line_width=2pt,default_color=plotColors[1])),
layer(x=t,y=1e4*x[:,2],Geom.path,Theme(line_width=2pt,default_color=plotColors[2])),
layer(x=t,y=x[:,3],Geom.path,Theme(line_width=2pt,default_color=plotColors[3])),Scale.x_log10,
Guide.xlabel("time (s)"), Guide.ylabel("Species concentration"),
Guide.manual_color_key("Legend",["Specie 1", "(Specie 2)*1e4", "Specie 3" ],plotColorsHex[1:3]),
Theme(panel_stroke=colorant"black",highlight_width=0.1pt,key_position=:top,
key_max_columns = 1,major_label_font_size=14pt,minor_label_font_size=10pt,
key_title_font_size=12pt,key_label_font_size=10pt),
Guide.xticks(ticks = [11:-1:-8;]),Guide.yticks(ticks=[1.0:-0.1:0.0;]));

draw(PNG("../../ImagesAndPDFs/Plots/RoberPlot.png",20cm,15cm),p)


