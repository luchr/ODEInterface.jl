
using ODEInterface
using Gadfly
using Colors
@ODEInterface.import_huge
loadODESolvers();

function threebody(t,x,dx)
    with_bigfloat_precision(113) do
        mu = parse(BigFloat,"0.012277471"); ms = 1 - mu;
        r1 = vecnorm(x[1:2]-[-mu,0]);
        r2 = vecnorm(x[1:2]-[ ms,0]);
        dx[1] = x[3];
        dx[2] = x[4];
        dx[3] = x[1] + 2*x[4] - ms*(x[1]+mu)/r1^3 - mu*(x[1]-ms)/r2^3;
        dx[4] = x[2] - 2*x[3] - ms*     x[2]/r1^3 - mu*     x[2]/r2^3;
    
        return nothing
    end
end

set_bigfloat_precision(113);
opt=OptionsODE(OPT_EPS => 1.11e-16,OPT_RTOL => 1e-16,OPT_ATOL=>1e-16,OPT_RHS_CALLMODE => RHS_CALL_INSITU);

t0 = 0.0; T = parse(BigFloat,"17.0652165601579625588917206249");
x0=[0.994, 0.0, 0.0, parse(BigFloat,"-2.00158510637908252240537862224")];

(t,x,retcode,stats) = odecall(dop853, threebody, [t0, T], x0, opt);

plotColors = ["black", "#4d5991", "#4bdf8c"];

arenstorfPlot = plot(layer(x=x[:,1],y=x[:,2],Geom.point,
Theme(default_point_size=2pt,default_color=parse(Colorant,plotColors[1])),order=1),
layer(x=[0],y=[0],Geom.point,
Theme(default_point_size=10pt,default_color=parse(Colorant,plotColors[2])),order=2),
layer(x=[0.994],y=[0],Geom.point,
Theme(default_point_size=8pt,default_color=parse(Colorant,plotColors[3])),order=2),
Theme(major_label_font_size=30pt,panel_stroke=colorant"black",highlight_width=0.1pt,
minor_label_font_size=28pt,key_title_font_size=30pt,key_label_font_size=24pt,key_position=:top,
key_max_columns = 1),
Coord.Cartesian(xmin=-1.4,xmax=1.1),Guide.xlabel("x",orientation=:horizontal),
Guide.manual_color_key("Legend",["Spaceship", "Earth", "Moon"],["black", "#4d5991", "#4bdf8c"]))

draw(SVGJS(30cm, 30cm), arenstorfPlot)
