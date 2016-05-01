# graphics routine for the output of the state variables and Lagrange
# multipliers for the reentry problem
# following Stoer/Bulirsch 1973

# also suitable for representing the solution of the auxiliary problem with the
# associated initial guesses of the Lagrange multipliers

# output J is an approximation of the value of the functional to be minimized
# along the given trajectory

# authors: Folkmar Bornemann, Vishal Sontakke, 2016/04/23

using Gadfly
include("utility.jl");

function showSolution(t,x,T,aux,plotFlag,color,fName)
    K = 1000;
    M = length(t);
    J = 0.0;

    for j = 1:M-1

      # calculating the piecewise trajectories

      tSpan = Array{Float64}(K);
      tSpan[:] = linspace(t[j],t[j+1],K);

      if aux
      	# Solver options
        # IVP Solver options
  			opt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
  			OPT_RTOL => tol, OPT_ATOL => tol);

        (tGrid,xGrid,_,_) =
            odecall(odesolver,aux_f,tSpan,x[:,j],opt);

        xGrid = xGrid';

        if dir == "backward"
            tGrid = 1-tGrid[end:-1:1];
            xGrid = xGrid[:,end:-1:1,:];
        end
        u = Array{Float64}(length(tGrid));
        for k=1:length(tGrid)
            (xGrid[4:6,k],u[k]) = multiplier_start(tGrid[k],xGrid[:,k]);
        end
      else
        # solver options
        opt = OptionsODE(OPT_RHS_CALLMODE => RHS_CALL_INSITU,
  			OPT_RTOL => tol, OPT_ATOL => tol);
        
        (tGrid,xGrid,_,_) =
                odecall(odesolver,reentry_f,tSpan,x[:,j],opt);
        xGrid = xGrid';
      end

      # evaluation of Phi, Hamilton function and control
      H = zeros(tGrid);
      Phi = zeros(tGrid);
      u = zeros(tGrid);
      for k = 1:length(H)
          (Phi[k],H[k],u[k]) = utility(xGrid[:,k]);
      end

      # trapezoidal rule for the approximation of functions in J in unit time
      for k = 1:(length(tGrid)-1)
          J = J + 0.5*T*(tGrid[k+1]-tGrid[k])*(Phi[k]+Phi[k+1]);
      end

      # convert units
      xGrid[2,:] = 180/pi*xGrid[2,:]; # Flight-path angle in Degrees
      xGrid[3,:] = R*xGrid[3,:];      # Height in 1e5 ft

      # plot the 3 state variables and the 3 Lagrange multipliers
      lineWidth = 3pt;

      for k=1:6
          push!(plotVar[k+(k>3)],layer(x=tGrid,y=xGrid[k,:],
          	Geom.line,Theme(line_width = lineWidth,default_color=color)))
      end

      # plot the control
      push!(plotVar[4],layer(x=tGrid,y=180/pi*u,
      	Geom.line,Theme(line_width = lineWidth,default_color=color)));

      # plot the Hamilton-function
      push!(plotVar[8],layer(x=tGrid,y=H,
        Geom.line,Theme(line_width = lineWidth,default_color=color)));
    end

    plotVar[1].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=0.26,ymax=0.38);
    plotVar[2].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=-10,ymax=6);
    plotVar[3].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=1.5,ymax=4.5);
    plotVar[4].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=-90,ymax=90);
    plotVar[5].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=-2,ymax=-0.4);
    plotVar[6].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=-1,ymax=2.5);
    plotVar[7].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=-15,ymax=15);
    plotVar[8].coord = Gadfly.Coord.Cartesian(xmin=0,xmax=1,ymin=-2.5e-10,ymax=2.5e-10);

    if plotFlag
      minFontSize = 16pt;
      majFontSize = 18pt;
    	for i=1:8
			plotVar[i].theme = Gadfly.Theme(panel_stroke = colorant"black",
			    minor_label_font_size=minFontSize,
		        major_label_font_size=majFontSize);
		end
      draw(PDF(fName,1855pt,1001pt),vstack(hstack(plotVar[1],plotVar[2],plotVar[3],plotVar[4]),hstack(plotVar[5],plotVar[6],plotVar[7],plotVar[8])))
    end

    return J
end
