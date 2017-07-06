"""
  Module for testing ODEInterface
  """
module ODEInterfaceTest

using Base.Test

using ODEInterface
@ODEInterface.import_huge

if VERSION >= v"0.6.0-dev"
  # @testloop was merged with @testset 93502e0f7
  macro testloop(ex::Expr)
    quote
      @testset( $(esc(ex))  )
    end
  end
else
  # v0.4 compatibility
  if !isdefined(Base.Test, Symbol("@testset"))
    macro testset(name, start_tag)
      return Expr(:block, :( println("Testing ",$name) ), start_tag)
    end
  end
  
  if !isdefined(Base.Test, Symbol("@testloop"))
    macro testloop(ex::Expr)
      quote
        $(esc(ex))
      end
    end
  end
end


const dl_solvers = (DL_DOPRI5, DL_DOPRI5_I32, 
                    DL_DOP853, DL_DOP853_I32,
                    DL_ODEX, DL_ODEX_I32,
                    DL_RADAU5, DL_RADAU5_I32, DL_RADAU, DL_RADAU_I32,
                    DL_SEULEX, DL_SEULEX_I32, 
                    DL_RODAS, DL_RODAS_I32,
                    DL_BVPSOL, DL_BVPSOL_I32,
                    DL_DDEABM, DL_DDEABM_I32,
                    DL_DDEBDF, DL_DDEBDF_I32,
                    ) 
const solvers = (dopri5, dopri5_i32, 
                 dop853, dop853_i32,
                 odex, odex_i32, 
                 radau5, radau5_i32, radau, radau_i32,
                 seulex, seulex_i32, 
                 rodas, rodas_i32,
                 ddeabm, ddeabm_i32,
                 ddebdf, ddebdf_i32,
                )

const solvers_without_dense_output = (ddeabm, ddeabm_i32, ddebdf, ddebdf_i32)

const solvers_without_special_struct_support = (ddebdf, ddebdf_i32)

const solvers_mas = ( radau5, radau5_i32, radau, radau_i32,
                      seulex, seulex_i32, 
                      rodas, rodas_i32,
                    )

const solvers_jac = ( radau5, radau5_i32, radau, radau_i32,
                      seulex, seulex_i32, 
                      rodas, rodas_i32,
                      ddebdf, ddebdf_i32,
                    )

const solvers_rhsdt = ( rodas, rodas_i32
                      )

const solvers_bv  = ( bvpsol, bvpsol_i32 
                    )

function test_ode1(solver::Function)
  opt = OptionsODE("ode1",
        OPT_RTOL => 1e-10,
        OPT_ATOL => 1e-10)
  t0 = 0; T = 1; x0=[1,2]; rhs = (t,x) -> x

  (t,x,retcode,stats) = solver(rhs, t0, T, x0, opt)
  @assert 1 == retcode
  @assert t == T
  @assert x0 == [1,2]
  @assert isapprox(x[1],exp(1),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[2],2*exp(1),rtol=1e-7,atol=1e-7)
  if haskey(stats,"step_predict")
    @assert isa(stats["step_predict"],Number)
  end
  return true
end

function test_ode2(solver::Function)
  dense_flag = !(solver in solvers_without_dense_output)
  x3 = NaN
  called_init = false
  called_done = false

  function outputfcn(reason,told,t,x,eval_sol_fcn,extra_data)
    if reason == OUTPUTFCN_CALL_INIT 
      called_init = true
      extra_data["test_ode2_data"] = 56
    end
    if reason == OUTPUTFCN_CALL_DONE 
      called_done = true
      @assert extra_data["test_ode2_data"] == 56
    end
    if reason == OUTPUTFCN_CALL_STEP
      if told ≤ 3.0 ≤ t
        if dense_flag
          x3 = eval_sol_fcn(3.0)[1]
        end
        return OUTPUTFCN_RET_STOP
      end
    end
    return OUTPUTFCN_RET_CONTINUE
  end
  opt = OptionsODE("ode2",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8,
        OPT_OUTPUTFCN => outputfcn,
        OPT_OUTPUTMODE => dense_flag ? OUTPUTFCN_DENSE : OUTPUTFCN_WODENSE,
        )
  t0 = 0; T = 5000; x0=[1,2]; rhs = (t,x) -> x
  (t,x,retcode,stats) = solver(rhs, t0, T, x0, opt)
  @assert called_init && called_done
  @assert 2 == retcode
  @assert t < T
  if dense_flag
    @assert isapprox(x3,exp(3),rtol=1e-7,atol=1e-7)
  end
  return true
end

function test_ode3(solver::Function)
  opt = OptionsODE("ode3",
        OPT_RTOL => 1e-10,
        OPT_ATOL => 1e-10,
        OPT_RHS_CALLMODE => RHS_CALL_INSITU,
        )
  function rhs(t,x,dx)
    dx[1]=x[1]; dx[2]=x[2];
    return nothing
  end
  t0 = 0; T = 1; x0=[1,2]; 

  (t,x,retcode,stats) = solver(rhs, t0, T, x0, opt)
  @assert 1 == retcode
  @assert t == T
  @assert x0 == [1,2]
  @assert isapprox(x[1],exp(1),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[2],2*exp(1),rtol=1e-7,atol=1e-7)
  return true
end

function test_massode1(solver::Function)
  mas = [ 2.0 1.0; 1.0 2.0]
  x1_exact = t -> 1.5*exp(t/3)-0.5*exp(t)
  x2_exact = t -> 1.5*exp(t/3)+0.5*exp(t)
  opt = OptionsODE("massode1",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8,
        OPT_MASSMATRIX => mas,
        )
  t0 = 0; T = 1; x0=[1,2]; rhs = (t,x) -> x
  (t,x,retcode,stats) = solver(rhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[1],x1_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[2],x2_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_massode2(solver::Function)
  mas = BandedMatrix(5,5, 1,1, 0.0)
  setdiagonal!(mas,0,2); setdiagonal!(mas,1,1); setdiagonal!(mas,-1,1)

  x1_exact = t -> 0.5*(-3*exp(t/3)+2*exp(t/2)-exp(t)+
                        (2+sqrt(3))*exp((2-sqrt(3))*t)-
                        (sqrt(3)-2)*exp((2+sqrt(3))*t) )
  x3_exact = t -> -exp(t/2)+(2+sqrt(3))*exp((2-sqrt(3))*t) - 
                  (sqrt(3)-2)*exp((2+sqrt(3))*t)

  opt = OptionsODE("massode2",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8,
        OPT_MASSMATRIX => mas,
        )
  t0 = 0; T = 1; x0=[1,2,3,4,5]; rhs = (t,x) -> x

  (t,x,retcode,stats) = solver(rhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[1],x1_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[3],x3_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_massode3(solver::Function)
  mas = [2 1 ; 1 2 ]

  x5_exact = t -> -0.5*exp(t/3)*(exp(2*t/3)-11)
  x6_exact = t ->  0.5*exp(t/3)*(exp(2*t/3)+11)

  opt = OptionsODE("massode3",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8,
        OPT_M1       => 4,
        OPT_M2       => 2,
        OPT_MASSMATRIX => mas,
        )
  t0 = 0; T = 1; x0=[1,2,3,4,5,6]; 
  rhs = (t,x) -> [ x[5],x[6] ]

  (t,x,retcode,stats) = solver(rhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[5],x5_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[6],x6_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_massode4(solver::Function)
  mas = BandedMatrix(3,3, 1,1, 0.0)
  setdiagonal!(mas,0,2); setdiagonal!(mas,1,1); setdiagonal!(mas,-1,1)

  x4_exact = t -> 4*exp(t)*(cosh(t/sqrt(2))-sqrt(2)*sinh(t/sqrt(2)))

  opt = OptionsODE("massode4",
        OPT_RTOL => 1e-10,
        OPT_ATOL => 1e-10,
        OPT_M1       => 2,
        OPT_M2       => 2,
        OPT_MASSMATRIX => mas,
        )
  t0 = 0; T = 1; x0=[1,2,3,4,5]; 
  rhs = (t,x) -> [ x[3],x[4],x[5] ]

  (t,x,retcode,stats) = solver(rhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[4],x4_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_jacode1(solver::Function)
  x1_exact = t -> exp(2*t)/(2-exp(2*t))
  x2_exact = t -> 2*sqrt(1.0/(2-exp(2*t)))
  
  function myjac(t,x,J)
    @assert isa(J,Array{Float64})
    J[1,1] = x[2]^2
    J[1,2] = 2*x[1]*x[2]
    J[2,1] = x[2]
    J[2,2] = x[1]
  end

  opt = OptionsODE("odejac1",
        OPT_RTOL => 1e-10,
        OPT_ATOL => 1e-10,
        OPT_JACOBIMATRIX => myjac,
        )
  t0 = 0; T = 0.2; x0=[1,2]; 
  (t,x,retcode,stats) = solver( (t,x)-> [x[1]*(x[2])^2,x[1]*x[2]], 
                              t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[1],x1_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[2],x2_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_jacode2(solver::Function)
  function myrhs(t,x)
    return [ x[1], x[1]+x[2], x[2]+x[3], x[3]+x[4] ]
  end
  
  function myjac(t,x,J)
    @assert isa(J,BandedMatrix{Float64})
    setdiagonal!(J,0,1.0)
    setdiagonal!(J,-1,1.0)
  end

  x1_exact = t -> exp(t)
  x2_exact = t -> (2+t)*exp(t)
  x4_exact = t -> exp(t)/6*(24+18*t+6*t^2+t^3)

  opt = OptionsODE("odejac2",
        OPT_RTOL => 1e-10,
        OPT_ATOL => 1e-10,
        OPT_JACOBIMATRIX => myjac,
        OPT_JACOBIBANDSTRUCT => (1,0),
        OPT_JACRECOMPFACTOR => -1,
        )
  t0 = 0; T = 1; x0=[1,2,3,4]; 

  (t,x,retcode,stats) = solver(myrhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[1],x1_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[2],x2_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[4],x4_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_jacode3(solver::Function)
  if solver ∈ solvers_without_special_struct_support
    return true
  end
  function myrhs(t,x)
    return [ # x[3], x[4],
             x[1]+x[2]+x[4],
             x[1]+x[2]+x[3] ]
  end
  
  function myjac(t,x,J)
    @assert isa(J,Array{Float64})
    @assert (2,4)==size(J)
    J[1,1] = 1; J[1,2] = 1; J[1,3]=0; J[1,4]=1;
    J[2,1] = 1; J[2,2] = 1; J[2,3]=1; J[2,4]=0;
  end
  
  x1_exact = t -> exp(-t)/3*(1-3*exp(t)+5*exp(3*t))
  x4_exact = t -> 2*exp(-t)/3*(1+5*exp(3*t))

  opt = OptionsODE("odejac3",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8,
        OPT_M1 => 2,
        OPT_M2 => 2,
        OPT_JACOBIMATRIX => myjac,
        OPT_JACOBIBANDSTRUCT => nothing,
        OPT_JACRECOMPFACTOR => -1,
        )
  t0 = 0; T = 1; x0=[1,2,3,4]; 

  (t,x,retcode,stats) = solver(myrhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[1],x1_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[4],x4_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_jacode4(solver::Function)
  if solver ∈ solvers_without_special_struct_support
    return true
  end
  function myrhs(t,x)
    return [ # x[3], x[4],
             2*x[1] + 2*x[3],
               x[2] +   x[4] ]
  end

  function myjac(t,x,J1,J2)
    @assert isa(J1,BandedMatrix{Float64}) && isa(J2,BandedMatrix{Float64})
    @assert (2,2)==size(J1) && (2,2)==size(J2)
    setdiagonal!(J1,0,[2,1])
    setdiagonal!(J2,0,[2,2])
  end

  hp = 1+sqrt(3); hm = 1-sqrt(3)

  x3_exact = t -> (-5*exp(hm*t)+3*sqrt(3)*exp(hm*t)+ 
                    5*exp(hp*t)+3*sqrt(3)*exp(hp*t) )/(2*sqrt(3))

  opt = OptionsODE("odejac3",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8,
        OPT_M1 => 2,
        OPT_M2 => 2,
        OPT_JACOBIMATRIX => myjac,
        OPT_JACOBIBANDSTRUCT => (0,0,),
        OPT_JACRECOMPFACTOR => -1,
        )
  t0 = 0; T = 1; x0=[1,2,3,4]; 

  (t,x,retcode,stats) = solver(myrhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[3],x3_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_rhstimederiv1(solver::Function)
  myrhs = (t,x) -> [ t*x[2], 4*t*x[1] ]

  function myjac(t,x,J)
    @assert isa(J,Array{Float64})
    @assert (2,2)==size(J)
    J[1,1] = 0; J[1,2] = t;
    J[2,1] = 4*t; J[2,2] = 0;
    return nothing
  end

  function myrhstimederiv(t,x,drhsdt)
    @assert isa(drhsdt,Array{Float64})
    @assert (2,)==size(drhsdt)
    drhsdt[1] = x[2]
    drhsdt[2] = 4*x[1]
    return nothing
  end

  x1_exact = t -> cosh(t*t) + 0.5*sinh(t*t)
  x2_exact = t -> cosh(t*t) + 2.0*sinh(t*t)

  opt = OptionsODE("odetimederiv1",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8,
        OPT_JACOBIMATRIX => myjac,
        OPT_RHSTIMEDERIV => myrhstimederiv,
        )
  t0 = 0; T = 1; x0=[1.0,1]

  (t,x,retcode,stats) = solver(myrhs, t0, T, x0, opt)
  @assert 1==retcode
  @assert isapprox(x[1],x1_exact(T),rtol=1e-7,atol=1e-7)
  @assert isapprox(x[2],x2_exact(T),rtol=1e-7,atol=1e-7)
  return true
end

function test_odecall1(solver::Function)
  opt = OptionsODE("odecall1",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8)
  t = [linspace(0,1,10)...]
  x0=[1,2]; rhs = (t,x) -> x
  (tVec,xVec,retcode,stats) = odecall(solver,rhs,t,x0,opt)
  @assert 1 == retcode
  @assert t == tVec
  return true
end

function test_odecall2(solver::Function)
  opt = OptionsODE("odecall2",
        OPT_RTOL => 1e-8,
        OPT_ATOL => 1e-8)
  t = [0,1]
  x0=[1,2]; rhs = (t,x) -> x
  (tVec,xVec,retcode,stats) = odecall(solver,rhs,t,x0,opt)
  @assert 1 == retcode
  @assert length(tVec)>2
  return true
end

function test_bvp1(solver::Function)
  ivpopt = OptionsODE("ivpoptions",
                   OPT_RHS_CALLMODE => RHS_CALL_INSITU)
                 
  opt = OptionsODE("bvp1",
                   OPT_RHS_CALLMODE => RHS_CALL_INSITU,
                   OPT_MAXSTEPS     => 10,
                   OPT_RTOL         => 1e-6,
                   OPT_BVPCLASS     => 0,
                   OPT_SOLMETHOD    => 0,
                   OPT_IVPOPT       => ivpopt)

  tNodes = [0,5]
  xInit = [ 5.0  0.45938665265299; 1  1];
  odesolver = nothing

  function f(t,x,dx)
    dx[1] =  x[2]
    dx[2] = -x[1]
    return nothing
  end

  function bc(xa,xb,r)
    r[1] = xa[1] - 5
    r[2] = xb[1] - 0.45938665265299
    return nothing
  end

  (t,x,code,stats) = solver(f,bc,tNodes,xInit,odesolver,opt)
  @assert tNodes == [0,5]
  @assert code>0
  @assert t == tNodes
  @assert isapprox(x[1,1],5.0,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[2,1],1.0,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[1,2],0.459387,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[2,2],5.07828,rtol=1e-4,atol=1e-4)
  
  return true
end

function test_bvp2(solver::Function)
  ivpopt = OptionsODE("ivpoptions",
                   OPT_RHS_CALLMODE => RHS_CALL_INSITU)
                 
  opt = OptionsODE("bvp2",
                   OPT_RHS_CALLMODE => RHS_CALL_INSITU,
                   OPT_MAXSTEPS     => 10,
                   OPT_RTOL         => 1e-6,
                   OPT_BVPCLASS     => 0,
                   OPT_SOLMETHOD    => 0,
                   OPT_IVPOPT       => ivpopt)

  tNodes = [0,5]
  xInit = [ 5.0  0.45938665265299; 1  1];
  odesolver = dop853

  function f(t,x,dx)
    dx[1] =  x[2]
    dx[2] = -x[1]
    return nothing
  end

  function bc(xa,xb,r)
    r[1] = xa[1] - 5
    r[2] = xb[1] - 0.45938665265299
    return nothing
  end

  (t,x,code,stats) = solver(f,bc,tNodes,xInit,odesolver,opt)
  @assert tNodes == [0,5]
  @assert code>0
  @assert t == tNodes
  @assert isapprox(x[1,1],5.0,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[2,1],1.0,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[1,2],0.459387,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[2,2],5.07828,rtol=1e-4,atol=1e-4)
  
  return true
end

function test_bvp3(solver::Function)
  ivpopt = OptionsODE("ivpoptions",
                   OPT_RHS_CALLMODE => RHS_CALL_INSITU)
                 
  opt = OptionsODE("bvp3",
                 OPT_RHS_CALLMODE => RHS_CALL_INSITU,
                 OPT_MAXSTEPS     => 100,
                 OPT_RTOL         => 1e-6,
                 OPT_BVPCLASS     => 2,
                 OPT_SOLMETHOD    => 1,
                 OPT_IVPOPT       => ivpopt)

  tNodes = collect(linspace(0,5,11))
  xInit = [ones(1,length(tNodes)-1) 0 ; ones(1,length(tNodes)) ]
  odesolver = dop853

  function f(t,x,dx)
    dx[1] = t*x[2]
    dx[2] = 4*max(0,x[1])^1.5
    return nothing
  end

  function bc(xa,xb,r)
    r[1] = xa[1] - 1
    r[2] = xb[1] - 0
    return nothing
  end

  (t,x,code,stats) = solver(f,bc,tNodes,xInit,odesolver,opt)
  @assert code>0
  @assert isapprox(x[1,1],1,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[2,1],-3.17614,rtol=1e-4,atol=1e-4)
  @assert isapprox(x[1,2],0.755201,rtol=1e-4,atol=1e-4)
  
  return true
end

function test_Banded()
  @testset "Banded" begin
    @test_throws ArgumentErrorODE  BandedMatrix{Float64}(
                                   5,4,1,2,zeros(Float64,(4,5)))
    @test_throws ArgumentErrorODE  BandedMatrix{Float64}(
                                   5,4,1,2,zeros(Float64,(5,4)))
    @test_throws ArgumentErrorODE  BandedMatrix{Float64}(
                                   5,4,5,1,zeros(Float64,(7,5)))
    
    bm = BandedMatrix(5,4, 1,2, NaN)
    @test   (bm[1,1] = 1.0) == 1.0
    @test   (bm[2,1] = 2  ) == 2.0
    @test   (bm[2,4] = 7  ) == 7.0
    @test   (bm[3,2] = 8  ) == 8.0
    @test_throws BoundsError  bm[1,4] = 1.1
    @test_throws BoundsError  bm[3,1] = 1.1
    @test_throws BoundsError  bm[4,1] = 1.1
    @test_throws BoundsError  bm[4,2] = 1.1
    @test_throws BoundsError  bm[5,3] = 1.1
 
    bm = BandedMatrix(5,4, 1,2, NaN)
    @test (bm[1:2,1]=[1,5]) == [1,5]
    @test (bm[1:3,2]=[4,2,4]) == [4,2,4]
    @test (bm[1:4,3]=[2,3,3,3]) == [2,3,3,3]
    @test (bm[2:5,4]=[1,0,4,2]) == [1,0,4,2]
    @test full(bm) == [1 4 2 0; 5 2 3 1; 0 4 3 0; 0 0 3 4; 0 0 0 2]
    bm_test = createBandedMatrix(
       [1.0 4 2 0; 5.0 2 3 1; 0 4.0 3 0; 0 0 3 4.0; 0 0 0 2.0]) 
    @test bm_test == bm

    bm = BandedMatrix(5,4, 1,2, NaN)
    @test setdiagonal!(bm,2,[2,1]) == [2,1]
    @test setdiagonal!(bm,1,[4,3,0]) == [4,3,0]
    @test setdiagonal!(bm,0,[1,2,3,4]) == [1,2,3,4]
    @test setdiagonal!(bm,-1,[5,4,3,2]) == [5,4,3,2]
    @test full(bm) == [1 4 2 0; 5 2 3 1; 0 4 3 0; 0 0 3 4; 0 0 0 2]

    bm = BandedMatrix(5,4, 1,2, NaN)
    diagonals = Any[ [2 1], [4 3 0], [1,2,3,4], [5 4 3 2]  ]
    @test setdiagonals!(bm,diagonals) == diagonals
    @test full(bm) == [1 4 2 0; 5 2 3 1; 0 4 3 0; 0 0 3 4; 0 0 0 2]
    
    bm2 = BandedMatrix(5,4, 1,2, NaN)
    setdiagonals!(bm2,bm)
    @test bm ≢ bm2
    @test bm == bm2
  end
end

function test_Options()
  @testset "Options" begin
    opt1 = OptionsODE("test1");
    @test isa(opt1, OptionsODE)
    @test isa(opt1, ODEInterface.AbstractOptionsODE)
    @test setOption!(opt1,"test_key",56) == nothing
    @test setOption!(opt1,"test_key",82) == 56
    @test getOption(opt1,"test_key",0) == 82
    @test getOption(opt1,"nokey",nothing) == nothing
    @test getOption(opt1,"nokey","none") == "none"
    @test setOptions!(opt1, "test_key" => 100, "new_key" => "bla") ==
          [82,nothing]

    opt2 = OptionsODE("test2",opt1)
    @test getOption(opt1,"test_key",0) == 100
  end
end

function test_DLSolvers()
  @testset "DLSolvers" begin
    result = loadODESolvers()
    @testloop for dl in dl_solvers
      @test result[dl].error == nothing
      @test result[dl].libhandle ≠ C_NULL
    end
    
    @testloop for dl in dl_solvers 
      @testloop for method in result[dl].methods
        @test method.error == nothing
        @test method.method_ptr ≠ C_NULL
        @test method.generic_name ≠ ""
        @test method.methodname_found ≠ ""
      end
    end
  end
end

function test_vanilla()
  problems = (test_ode1,test_ode2,test_ode3,)
  @testset "solvers" begin
    @testloop for solver in solvers,
                  problem in problems
      @test problem(solver)
    end
  end
end

function test_mas_solvers()
  problems = (test_massode1,test_massode2,test_massode3,test_massode4)
  @testset "mas-solvers" begin
    @testloop for solver in solvers_mas,
                  problem in problems
      @test problem(solver)
    end
  end
end

function test_jac_solvers()
  problems = (test_jacode1,test_jacode2,test_jacode3,test_jacode4)
  @testset "jac-solvers" begin
    @testloop for solver in solvers_jac,
                  problem in problems
      @test problem(solver)
    end
  end
end

function test_rhsdt_solvers()
  problems = (test_rhstimederiv1,)
  @testset "rhs_dt-sol." begin
    @testloop for solver in solvers_rhsdt,
                  problem in problems
      @test problem(solver)
    end
  end
end

function test_solvers()
  test_vanilla()
  test_mas_solvers()
  test_jac_solvers()
  test_rhsdt_solvers()
end

function test_odecall()
  problems = (test_odecall1,test_odecall2,)
  @testset "odecall" begin
    @testloop for solver in solvers,
                  problem in problems
      @test problem(solver)
    end
  end
end

function test_bvp()
  problems = (test_bvp1,test_bvp2,test_bvp3,)
  @testset "bvp" begin
    @testloop for solver in solvers_bv,
                  problem in problems
      @test problem(solver)
    end
  end
end

function test_all()
  test_Banded()
  test_Options()
  test_DLSolvers()
  test_solvers()
  test_odecall()
  test_bvp()
end

test_all()

end

# vim:syn=julia:cc=79:fdm=indent:
