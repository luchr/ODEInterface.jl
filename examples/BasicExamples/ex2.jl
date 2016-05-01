using ODEInterface
@ODEInterface.import_huge

loadODESolvers()

opt=OptionsODE(
  OPT_RTOL => 1e-8,
  OPT_ATOL => 1e-3)

(t,x,retcode,stats) = odecall( dopri5_i32, (t,x) -> x, 
  collect(linspace(0,1,15)),[1,2], opt)

println("retcode=$retcode")
println("t=$t");
println("x=",x);
println("stats:");display(stats);println();


# vim:syn=julia:cc=79:fdm=indent:
