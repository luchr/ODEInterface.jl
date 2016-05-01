using ODEInterface
@ODEInterface.import_huge

loadODESolvers()

function myoutputfcn(reason,told,t,x,eval_sol_func,extra_data)
  println("myoutputfcn called with reason=$reason, told=$told, t=$t, x=$x")
  if reason == OUTPUTFCN_CALL_STEP
    tmid = told + 0.5*(t-told)
    xmid = eval_sol_func(tmid)
    println("solution at $tmid: $xmid")
  end
  return OUTPUTFCN_RET_CONTINUE
end

opt=OptionsODE(
  OPT_RTOL       => 1e-8,
  OPT_ATOL       => 1e-3,
  OPT_OUTPUTFCN  => myoutputfcn,
  OPT_OUTPUTMODE => OUTPUTFCN_DENSE,
  )

(t,x,retcode,stats) = odecall( dop853, (t,x) -> x,[0,2],[1,2], opt)

println("retcode=$retcode")
println("t=$t");
println("x=",x);
println("stats:");display(stats);println();


# vim:syn=julia:cc=79:fdm=indent:
