using ODEInterface
@ODEInterface.import_huge

loadODESolvers()


function testrhs(t,x,dx)
  dx[1] = x[2]
  dx[2] = t*x[1]
  return nothing
end

function myjac(t,x,J)
  J[1,1]=0; J[1,2]=1;
  J[2,1]=t; J[2,2]=0;
  return nothing
end

opt = OptionsODE("airytest",
      OPT_RTOL     => 1e-9,
      OPT_ATOL     => 1e-9,
      OPT_MAXSTEPS => 1e6,
      OPT_RHS_CALLMODE => RHS_CALL_INSITU,
      OPT_JACOBIMATRIX => myjac,
      OPT_JACRECOMPFACTOR => -1,
      )

t0 = 0
T = -1000
x0 = [ airyai(t0), airyaiprime(t0) ]

@time  (t,x,retcode,stats) = radau5(testrhs,t0,T,x0,opt)
@time  (t,x,retcode,stats) = radau5(testrhs,t0,T,x0,opt)

@assert retcode==1
@assert isapprox(x[1],airyai(T),rtol=1e-6,atol=1e-6)
dump(stats)

# Dop853:
# 1st: 3.120123 seconds (15.58 M allocations: 621.371 MB, 3.29% gc time)
# to
# 1st: 1.951694 seconds (5.09 M allocations: 171.603 MB, 1.72% gc time)
# 2nd: 0.321590 seconds (3.99 M allocations: 121.809 MB, 4.94% gc time)

# Radau5:
# 1st: 4.720030 seconds (16.58 M allocations: 583.899 MB, 2.18% gc time)
# 2nd: 2.374964 seconds (14.99 M allocations: 514.093 MB, 3.23% gc time)
# to
# 1st: 4.390065 seconds (13.73 M allocations: 444.348 MB, 1.56% gc time)
# 2nd: 2.081976 seconds (12.12 M allocations: 374.365 MB, 1.77% gc time)

# vim:syn=julia:cc=79:fdm=indent:
