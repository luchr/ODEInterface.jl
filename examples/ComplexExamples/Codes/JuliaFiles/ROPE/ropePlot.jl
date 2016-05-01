
# Load required packages
using Gadfly
using Colors
using ODEInterface
using HDF5
@ODEInterface.import_huge
loadODESolvers();

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
(t,x_ref,retcode,stats) = odecall(dop853,rope,[t0, T], x0, opt);

x = zeros(length(t),n+1);
y = zeros(length(t),n+1);
l=1/n;

for j = 1:length(t)
    for i=1:n
        x[j,i+1] = x[j,i]+l*sin(x_ref[j,i]);
        y[j,i+1] = y[j,i]+l*cos(x_ref[j,i]);
    end
end

h5write("ropePlotx.h5","x",x)
h5write("ropePloty.h5","y",y)
h5write("ropePlott.h5","t",t)


