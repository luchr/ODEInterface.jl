# Functions used for several Hairer-Wanner Solvers

"""
       function hw1rhs(n,t,x,f,cbi::ODEinternalCallInfos)
  
  This function calls `rhs` saved in `...InternalCallInfos`.
  """
function hw1rhs(n,t,x,f,cbi::ODEinternalCallInfos)
  lprefix = cbi.rhs_lprefix

  (lio,l)=(cbi.logio,cbi.loglevel)
  l_rhs = l & LOG_RHS>0

  l_rhs && println(lio,lprefix,"called with n=",n," t=",t," x=",x)
  if cbi.rhs_mode == RHS_CALL_INSITU
    cbi.rhs(t,x,f)
  elseif cbi.rhs_mode == RHS_CALL_RETURNS_ARRAY
    result=cbi.rhs(t,x)
    try
      f[:] = getVectorCheckLength(result,Float64,n,false)
    catch e
      throw(OutputErrorODE(
        string("Cannot convert result $result of ",cbi.rhs," to ",
               "Vector{Float64} of length $n"),cbi.rhs,e))
    end
  else
    throw(InternalErrorODE("Unknown cbi.rhs_mode"))
  end
  l_rhs && println(lio,lprefix,"rhs result=",f)
  return nothing
end

"""
       function unsafe_HW1RHSCallback{FInt}(n_::Ptr{FInt}, t_::Ptr{Float64},
         x_::Ptr{Float64}, f_::Ptr{Float64}, rpar_::Ptr{Float64}, 
         ipar_::Ptr{FInt})
  
  This is the right-hand side given as callback to several Fortran-solvers,
  e.g. dopri5, dop853, odex.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  
  Uses hw1rhs.
  """
function unsafe_HW1RHSCallback{FInt}(n_::Ptr{FInt}, t_::Ptr{Float64},
  x_::Ptr{Float64}, f_::Ptr{Float64}, rpar_::Ptr{Float64}, 
  ipar_::Ptr{FInt})

  n = unsafe_load(n_); t = unsafe_load(t_)
  x = unsafe_wrap(Array, x_,(n,),false)
  f = unsafe_wrap(Array, f_,(n,),false)
  cid = unpackUInt64FromPtr(ipar_)
  cbi = get(GlobalCallInfoDict,cid,nothing)
  cbi == nothing && throw(InternalErrorODE(
      string("Cannot find call-id ",int2logstr(cid[1]),
             " in GlobalCallInfoDict")))
  
  hw1rhs(n,t,x,f,cbi)
  return nothing
end

"""
  `cfunction` pointer for unsafe_HW1RHSCallback with 64bit integers.
  """
const unsafe_HW1RHSCallback_c = cfunction(
  unsafe_HW1RHSCallback, Void, (Ptr{Int64},Ptr{Float64},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}))

"""
  `cfunction` pointer for unsafe_HW1RHSCallback with 32bit integers.
  """
const unsafe_HW1RHSCallbacki32_c = cfunction(
  unsafe_HW1RHSCallback, Void, (Ptr{Int32},Ptr{Float64},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}))

"""
       function hw2rhs(n,t,x,f,cbi::ODEinternalCallInfos)
  
  This function calls `rhs` saved in `...InternalCallInfos`.
  """
function hw2rhs(n,t,x,f,cbi::ODEinternalCallInfos)
  lprefix = cbi.rhs_lprefix

  (lio,l)=(cbi.logio,cbi.loglevel)
  l_rhs = l & LOG_RHS>0

  l_rhs && println(lio,lprefix,"called with n=",n," t=",t," x=",x)
  if cbi.M1 > 0
    # Problem with special structure, fill in special structure in f
    for k in 1:cbi.M1
      f[k] = x[k+cbi.M2]
    end
  end
  if cbi.rhs_mode == RHS_CALL_INSITU
    cbi.rhs(t,x,f)
  elseif cbi.rhs_mode == RHS_CALL_RETURNS_ARRAY
    result=cbi.rhs(t,x)
    try
      f[(cbi.M1+1):n] = getVectorCheckLength(result,Float64,n-cbi.M1,false)
    catch e
      throw(OutputErrorODE(
        string("Cannot convert result $result of ",cbi.rhs," to ",
               "Vector{Float64} of length d-M1=",n,"-",cbi.M1,"=",n-cbi.M1),
        cbi.rhs,e))
    end
  else
    throw(InternalErrorODE("Unknown cbi.rhs_mode"))
  end
  l_rhs && println(lio,lprefix,"rhs result=",f)
  return nothing
end

"""
       function unsafe_HW2RHSCallback{FInt}(n_::Ptr{FInt}, t_::Ptr{Float64},
         x_::Ptr{Float64}, f_::Ptr{Float64}, rpar_::Ptr{Float64}, 
         ipar_::Ptr{FInt})
  
  This is the right-hand side given as callback to Fortran-solvers
  (e.g. radau5 and radau) that can handle problems with "special structure", 
  see `help_specialstructure`.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  
  Uses hw2rhs.
  """
function unsafe_HW2RHSCallback{FInt}(n_::Ptr{FInt}, t_::Ptr{Float64},
  x_::Ptr{Float64}, f_::Ptr{Float64}, rpar_::Ptr{Float64}, 
  ipar_::Ptr{FInt})

  n = unsafe_load(n_); t = unsafe_load(t_)
  x = unsafe_wrap(Array, x_,(n,),false)
  f = unsafe_wrap(Array, f_,(n,),false)
  cid = unpackUInt64FromPtr(ipar_)
  cbi = get(GlobalCallInfoDict,cid,nothing)
  cbi == nothing && throw(InternalErrorODE(
      string("Cannot find call-id ",int2logstr(cid[1]),
             " in GlobalCallInfoDict")))

  hw2rhs(n,t,x,f,cbi)
  return nothing
end

"""
  `cfunction` pointer for unsafe_HW2RHSCallback with 64bit integers.
  """
const unsafe_HW2RHSCallback_c = cfunction(
  unsafe_HW2RHSCallback, Void, (Ptr{Int64},Ptr{Float64},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int64}))

"""
  `cfunction` pointer for unsafe_HW2RHSCallback with 32bit integers.
  """
const unsafe_HW2RHSCallbacki32_c = cfunction(
  unsafe_HW2RHSCallback, Void, (Ptr{Int32},Ptr{Float64},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}))

"""
       function unsafe_HW1MassCallback{FInt}(n_::Ptr{FInt}, 
                   am_::Ptr{Float64}, lmas_::Ptr{FInt}, 
                   rpar_::Ptr{Float64}, ipar_::Ptr{FInt})
  
  This is the MAS callback given to radau5, radau and seulex.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.
  
  This function takes the values of  the mass matrix saved in 
  the InternalCallInfos.
  """
function unsafe_HW1MassCallback{FInt}(n_::Ptr{FInt}, am_::Ptr{Float64},
            lmas_::Ptr{FInt}, rpar_::Ptr{Float64}, ipar_::Ptr{FInt})
  n = unsafe_load(n_)
  lmas = unsafe_load(lmas_)
  am = unsafe_wrap(Array, am_,(lmas,n,),false)
  cid = unpackUInt64FromPtr(ipar_)
  lprefix = string(int2logstr(cid),"unsafe_HW1MassCallback: ")
  cbi = get(GlobalCallInfoDict,cid,nothing)
  cbi==nothing && throw(InternalErrorODE(
      string("Cannot find call-id ",int2logstr(cid[1]),
             " in GlobalCallInfoDict")))

  (lio,l)=(cbi.logio,cbi.loglevel)
  l_mas = l & LOG_MASS>0
  
  l_mas && println(lio,lprefix,"called with n=",n," lmas=",lmas)

  mas = cbi.massmatrix

  if isa(mas,BandedMatrix)
    @assert n == size(mas,2)
    @assert lmas == 1+mas.l+mas.u
    bm = BandedMatrix{Float64}(mas.m,mas.n, mas.l, mas.u, am::Matrix{Float64})
    setdiagonals!(bm,mas)
  else
    @assert (lmas,n) == size(mas)
    am[:] = mas
  end
  
  l_mas && println(lio,lprefix,"am=",am)
  return nothing
end

"""
  `cfunction` pointer for unsafe_HW1MassCallback with 32bit integers.
  """
const unsafe_HW1MassCallbacki32_c = cfunction(
  unsafe_HW1MassCallback, Void, (Ptr{Int32},
    Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}))

"""
  `cfunction` pointer for unsafe_HW1MassCallback with 64bit integers.
  """
const unsafe_HW1MassCallback_c = cfunction(
  unsafe_HW1MassCallback, Void, (Ptr{Int64},
    Ptr{Float64}, Ptr{Int64}, Ptr{Float64}, Ptr{Int64}))


"""
       function extractSpecialStructureOpt{FInt}(d::FInt,
                  opt::AbstractOptionsODE) -> (M1,M2,NM1)

  extracts parameters for special structure (M1, M2).

  reads options: `OPT_M1`, `OPT_M2`
  """
function extractSpecialStructureOpt{FInt}(d::FInt,opt::AbstractOptionsODE)
  OPT = nothing
  M1 = 0; M2 = 0; NM1 = 0
  try
    OPT=OPT_M1; M1 = convert(FInt,getOption(opt,OPT,0))
    @assert M1 ≥ 0
    OPT=OPT_M2; M2 = convert(FInt,getOption(opt,OPT,M1))
    @assert M2 ≥ 0

    OPT=string(OPT_M1," & ",OPT_M2)
    @assert M1+M2 ≤ d
    @assert (M1==M2==0) || (M1≠0≠M2)
    @assert (M1==0) || (0 == M1 % M2)
    NM1 = d - M1
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end
  return (M1,M2,NM1)
end

"""
       function extractMassMatrix{FInt}(M1::FInt, M2::FInt, NM1::FInt,
           args::AbstractArgumentsODESolver,opt::AbstractOptionsODE)
  
  extracts mass matrix and fills `MAS`, `IMAS`, `MLMAS` und `MUMAS` in args.

  reads options: `OPT_MASSMATRIX`
  """
function extractMassMatrix{FInt}(M1::FInt, M2::FInt, NM1::FInt,
    args::AbstractArgumentsODESolver{FInt},opt::AbstractOptionsODE)
  OPT = nothing
  massmatrix = nothing
  try
    OPT = OPT_MASSMATRIX
    massmatrix = getOption(opt,OPT,nothing)
    if massmatrix == nothing
      args.IMAS = [0]; args.MLMAS = [0]; args.MUMAS=[0]
    else
      @assert 2==ndims(massmatrix)
      @assert (NM1,NM1,) == size(massmatrix)
      massmatrix = deepcopy(massmatrix)
      args.IMAS = [1]; 
      # A BandedMatrix with lower bandwidth == NM1 is treated as full!
      if isa(massmatrix,BandedMatrix) && massmatrix.l == NM1
        massmatrix=full(massmatrix)
      end
      if isa(massmatrix,BandedMatrix)
        @assert massmatrix.l < NM1
        args.MLMAS = [ massmatrix.l ]
        args.MUMAS = [ massmatrix.u ]
      else
        massmatrix = convert(Matrix{Float64},massmatrix)
        args.MLMAS = [ NM1 ]; args.MUMAS = [ NM1 ]
      end
    end
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end
  args.MAS = (FInt == Int64)? unsafe_HW1MassCallback_c:
                              unsafe_HW1MassCallbacki32_c
  return massmatrix
end

"""
       function unsafe_HW1JacCallback{FInt}(n_::Ptr{FInt},
               t_::Ptr{Float64},x_::Ptr{Float64},dfx_::Ptr{Float64},
               ldfx_::Ptr{FInt}, rpar_::Ptr{Float64}, ipar_::Ptr{FInt})
  
  This is the JAC callback given to radau5, radau and seulex.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.
  
  This function calls the user-given Julia function cbi.jacobimatrix
  with the appropriate arguments (depending on M1 and jacobibandstruct).
  """
function unsafe_HW1JacCallback{FInt}(n_::Ptr{FInt},
        t_::Ptr{Float64},x_::Ptr{Float64},dfx_::Ptr{Float64},
        ldfx_::Ptr{FInt}, rpar_::Ptr{Float64}, ipar_::Ptr{FInt})
  n = unsafe_load(n_)
  t = unsafe_load(t_)
  x = unsafe_wrap(Array, x_,(n,),false)
  ldfx = unsafe_load(ldfx_)
  cid = unpackUInt64FromPtr(ipar_)
  cbi = get(GlobalCallInfoDict,cid,nothing)
  cbi==nothing && throw(InternalErrorODE(
      string("Cannot find call-id ",int2logstr(cid[1]),
             " in GlobalCallInfoDict")))
  lprefix = cbi.jac_lprefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_jac = l & LOG_JAC>0
  
  l_jac && println(lio,lprefix,"called with n=",n," ldfx=",ldfx)
  jac = cbi.jacobimatrix
  jb = cbi.jacobibandstruct
  if jb == nothing
    @assert ldfx==n-cbi.M1
    J = unsafe_wrap(Array, dfx_,(ldfx,n,),false)
    jac(t,x,J)
  else
    @assert ldfx==1+jb[1]+jb[2]
    if cbi.M1==0
      J = BandedMatrix{Float64}(n,n, jb[1],jb[2], 
            unsafe_wrap(Array, dfx_,(ldfx,n,),false)::Matrix{Float64})
      jac(t,x,J)
    else
      no = Int(cbi.M1/cbi.M2+1)
      J = Vector{BandedMatrix{Float64}}(no)
      for k in 1:no
        ptr = dfx_+(k-1)*ldfx*cbi.M2*sizeof(Float64)
        darr =  unsafe_wrap(Array, ptr, (ldfx,cbi.M2,), false)
        J[k] = BandedMatrix{Float64}(cbi.M2,cbi.M2, jb[1],jb[2],darr)
      end
      jac(t,x,J...)
    end
  end
  
  l_jac && println(lio,lprefix,"dfx=",unsafe_wrap(Array, dfx_,(ldfx,n,),false))
  return nothing
end

"""
  `cfunction` pointer for unsafe_HW1JacCallback with 32bit integers.
  """
const unsafe_HW1JacCallbacki32_c = cfunction(
  unsafe_HW1JacCallback, Void, (Ptr{Int32},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
    Ptr{Int32}, Ptr{Float64}, Ptr{Int32}))

"""
  `cfunction` pointer for unsafe_HW1JacCallback with 64bit integers.
  """
const unsafe_HW1JacCallback_c = cfunction(
  unsafe_HW1JacCallback, Void, (Ptr{Int64},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, 
    Ptr{Int64}, Ptr{Float64}, Ptr{Int64}))

"""
        function unsafe_HWRhsTimeDerivCallback{FInt}(n_::Ptr{FInt},
                t_::Ptr{Float64},x_::Ptr{Float64},dft_::Ptr{Float64},
                rpar_::Ptr{Float64},ipar_::Ptr{FInt})
  
  This is the DFX callback given to rodas.

  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-pointers.
  
  This function calls the user-given Julia function cbi.rhstimederiv
  with the appropriate arguments.
  """
function unsafe_HWRhsTimeDerivCallback{FInt}(n_::Ptr{FInt},
        t_::Ptr{Float64},x_::Ptr{Float64},dfdt_::Ptr{Float64},
        rpar_::Ptr{Float64},ipar_::Ptr{FInt})
  n = unsafe_load(n_)
  t = unsafe_load(t_)
  x = unsafe_wrap(Array, x_,(n,),false)
  dfdt = unsafe_wrap(Array, dfdt_,(n,),false)
  cid = unpackUInt64FromPtr(ipar_)
  cbi = get(GlobalCallInfoDict,cid,nothing)
  cbi == nothing && throw(InternalErrorODE(
      string("Cannot find call-id ",int2logstr(cid[1]),
             " in GlobalCallInfoDict")))
  lprefix = cbi.rhsdt_prefix
  (lio,l)=(cbi.logio,cbi.loglevel)
  l_rhsdt = l & LOG_RHSDT>0

  l_rhsdt && println(lio,lprefix,"called with n=",n," t=",t)
  cbi.rhsdt(t,x,dfdt)
  l_rhsdt && println(lio,lprefix,"dfdt=",dfdt)

  return nothing
end

"""
  `cfunction` pointer for unsafe_HWRhsTimeDerivCallback with 32bit integers.
  """
const unsafe_HWRhsTimeDerivCallbacki32_c = cfunction(
  unsafe_HWRhsTimeDerivCallback, Void, (Ptr{Int32},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    Ptr{Float64}, Ptr{Int32}))

"""
  `cfunction` pointer for unsafe_HWRhsTimeDerivCallback with 64bit integers.
  """
const unsafe_HWRhsTimeDerivCallback_c = cfunction(
  unsafe_HWRhsTimeDerivCallback, Void, (Ptr{Int64},
    Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
    Ptr{Float64}, Ptr{Int64}))

"""
       function extractJacobiOpt{FInt}(d::FInt,M1::FInt, M2::FInt, NM1::FInt,
            cid_str, args::AbstractArgumentsODESolver,opt::AbstractOptionsODE)
  
  extracts jacobi options and
  fills `JAC`, `IJAC`, `MLJAC` and `MUJAC` in args.

  reads options: `OPT_JACOBIMATRIX`, `OPT_JACOBIBANDSTRUCT`
  """
function extractJacobiOpt{FInt}(d::FInt,M1::FInt, M2::FInt, NM1::FInt,
    cid_str, args::AbstractArgumentsODESolver{FInt},opt::AbstractOptionsODE)
  OPT = nothing
  jacobimatrix = nothing
  jacobibandstruct = nothing
  try
    OPT = OPT_JACOBIMATRIX
    jacobimatrix = getOption(opt,OPT,nothing)
    @assert (jacobimatrix == nothing) || isa(jacobimatrix,Function)
    
    if jacobimatrix ≠ nothing
      OPT = OPT_JACOBIBANDSTRUCT
      bs = getOption(opt,OPT,nothing)
      
      if bs ≠ nothing
        jacobibandstruct = ( convert(FInt,bs[1]), convert(FInt,bs[2]) )
        if jacobibandstruct[1] == NM1 
          # A BandedMatrix with lower bandwidth == NM1 is treated as full!
          jacobibandstruct = nothing
        end
      end
      if jacobibandstruct ≠ nothing
        @assert (M1==0) || (M1+M2==d)
        @assert 0 ≤ jacobibandstruct[1] < NM1
        @assert  (M1==0 && 0 ≤ jacobibandstruct[2] ≤ d)  ||
                 (M1>0  && 0 ≤ jacobibandstruct[2] ≤ M2)
      end
    end
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  args.IJAC = [ jacobimatrix==nothing? 0 : 1] 
  args.MLJAC = [ jacobibandstruct==nothing? d : jacobibandstruct[1]  ];
  args.MUJAC=[ jacobibandstruct==nothing? d : jacobibandstruct[2] ]
  args.JAC = (FInt == Int64)? unsafe_HW1JacCallback_c:
                              unsafe_HW1JacCallbacki32_c
  jac_lprefix = string(cid_str,"unsafe_HW1JacCallback: ")
  return (jacobimatrix,jacobibandstruct,jac_lprefix)
end

"""
        function extractRhsTimeDerivOpt{FInt}(cid_str,
                args::AbstractArgumentsODESolver{FInt}, 
                opt::AbstractOptionsODE)
  
  extracts options for callback function for time-derivatives 
  of the right-hand-side and
  fills `DFX` and `IDFX` in args.
  """
function extractRhsTimeDerivOpt{FInt}(cid_str,
        args::AbstractArgumentsODESolver{FInt}, opt::AbstractOptionsODE)
  OPT = OPT_RHSTIMEDERIV
  rhstimederiv = nothing
  try
    rhstimederiv = getOption(opt,OPT,nothing)
    @assert (rhstimederiv == nothing) || isa(rhstimederiv,Function)
  catch e
    throw(ArgumentErrorODE("Option '$OPT': Not valid",:opt,e))
  end

  args.IDFX = [rhstimederiv==nothing? 0 : 1]
  args.DFX = (FInt == Int64)? unsafe_HWRhsTimeDerivCallback_c:
                              unsafe_HWRhsTimeDerivCallbacki32_c
  rhsdt_prefix = string(cid_str,"unsafe_HWRhsTimeDerivCallback: ")
  return (rhstimederiv,rhsdt_prefix)
end

"""
  ## License
  
  This is the license text, which can also be fount at
  
       http://www.unige.ch/~hairer/prog/licence.txt
  
  Copyright (c) 2004, Ernst Hairer
  
  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are
  met:
  
  - Redistributions of source code must retain the above copyright 
  notice, this list of conditions and the following disclaimer.
  
  - Redistributions in binary form must reproduce the above copyright 
  notice, this list of conditions and the following disclaimer in the 
  documentation and/or other materials provided with the distribution.
  
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS 
  IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR 
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  """
const hw_license = nothing

# vim:syn=julia:cc=79:fdm=indent:
