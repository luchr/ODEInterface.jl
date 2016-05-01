
# Load all the required packages
using ODEInterface
using ForwardDiff
using Gadfly
using Colors
@ODEInterface.import_huge
loadODESolvers();

# Define the right-hand function for automatic differentiation
function vdpolAD(x)
    return [x[2],((1-x[1]^2)*x[2]-x[1])*1e6]
end

# Define the system for the solver
function vdpol(t,x,dx)
    dx[:] = vdpolAD(x);
    return nothing
end

# Define the Jacobian function using AD
function getJacobian(t,x,J)
    J[:,:] = ForwardDiff.jacobian(vdpolAD,x);
    return nothing
end

# Flag to check whether plot is to be generated and saved or not
# Also checks if all solvers are successful
printFlag = true;

# Initial conditions
t0 = 0.0; T = [1.0:11.0;]; x0 = [2.0,0.0];

# Get "reference solution"
Tol = 1e-14;
# for Tol < 1e-14 we get the error "TOLERANCES ARE TOO SMALL"
opt = OptionsODE(OPT_EPS=>1.11e-16,OPT_RTOL=>Tol, OPT_ATOL=>Tol,
OPT_RHS_CALLMODE => RHS_CALL_INSITU,
OPT_JACOBIMATRIX => getJacobian);

(t,x,retcode,stats) = odecall(seulex,vdpol,[t0, T[end]], x0, opt);
t = [t0;t];
x = [x0';x];

p = plot(x=t,y=x[:,1],Geom.path,
Theme(line_width=2pt,default_color=colorant"black",
panel_stroke=colorant"black",key_position=:top,
key_max_columns = 1,major_label_font_size=24pt,minor_label_font_size=22pt,
key_title_font_size=22pt,key_label_font_size=20pt),
Coord.cartesian(ymin=-2.5,ymax=2.5,xmin=0,xmax=11),
Guide.xlabel("time (s)"),Guide.ylabel("Position (x)"),
Guide.xticks(ticks=[0:11;]),Guide.yticks(ticks=[-2:0.5:2;]));

draw(PNG("../../ImagesAndPDFs/Plots/vdpolPlot.png",30cm,20cm),p);


