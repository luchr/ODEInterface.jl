# Functions common for serveral Hairer-Wanner Solvers

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
  x = pointer_to_array(x_,(n,),false)
  f = pointer_to_array(f_,(n,),false)
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
  
  This is the right-hand side given as callback to Fortran-solvers that
  can handle problems with "special structure", see `OPT_M1`, e.g. radau5.
  
  The `unsafe` prefix in the name indicates that no validations are 
  performed on the `Ptr`-arguments.
  
  Uses hw2rhs.
  """
function unsafe_HW2RHSCallback{FInt}(n_::Ptr{FInt}, t_::Ptr{Float64},
  x_::Ptr{Float64}, f_::Ptr{Float64}, rpar_::Ptr{Float64}, 
  ipar_::Ptr{FInt})

  n = unsafe_load(n_); t = unsafe_load(t_)
  x = pointer_to_array(x_,(n,),false)
  f = pointer_to_array(f_,(n,),false)
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

# vim:syn=julia:cc=79:fdm=indent:
