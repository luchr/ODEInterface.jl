# Dynamic loading of ODE-Solvers for ODEInterface

"""
  macro for importing (un-)load functions.
  """
macro import_dynamicload()
  :(
    using ODEInterface: loadODESolvers, unloadODESolvers
  )
end

"""
  Type describing a dynamically loaded method.
  """
immutable MethodDLinfo
  generic_name     :: AbstractString     # input for trytoloadmethod
  methodname_found :: AbstractString     # variant of generic_name found
  method_ptr       :: Ptr{Void}          # pointer to method/code
  error                                  # error if there was one
end

"""
  Type describing the "dynamic parts" of a solver.
  """
immutable SolverDLinfo
  libname          :: AbstractString     # the name of the dynamic library
  libfilepath      :: AbstractString     # path to lib, result of find_library
  libhandle        :: Ptr{Void}          # handle of lib, result of dlopen
                                         # an immutable Dict is missing ...
  methods          :: Tuple              #   ... Tuple of MethodDLinfo
  error                                  # error if there was one
end

"""
  global Dict saving informations of all loaded solvers.
  """
const dlSolversInfo = Dict{AbstractString,SolverDLinfo}()

"""
       function trytoloadlib(name::AbstractString,extrapaths::Vector)
  
  tries to (dynamically) load the given shared library given by name.
  
  if `ame="name"` then all the following variants will be tried:
  `"name"`, `"NAME"`, `"Name"`
  
  returns `(ptr,filepath)`
  """
function trytoloadlib(name::AbstractString,extrapaths::Vector)
  ptr = C_NULL

  trylist =  [name, uppercase(name), ucfirst(name)]
  filepath = Libdl.find_library(trylist, extrapaths)
  if isempty(filepath) 
    throw(ErrorException(
      "Cannot find one of $trylist in libpaths or in $extrapaths"))
  end
  ptr = Libdl.dlopen(filepath)
  return (ptr,filepath)
end

"""
       function trytoloadmethod(libhandle::Ptr{Void},
             method_name::AbstractString) -> (ptr,namefound)

  tries to find the given method by name in a dynamically loaded library.
  
  if `method_name="name"` then the following variants will be tried:
  `"name"`, `"NAME"`, `"Name"`,
  `"name_"`, `"NAME_"`, `"Name_"`,
  `"_name"`, `"_NAME"`, `"_Name"`,
  `"_name_"`, `"_NAME_"`, `"_Name_"`
  """
function trytoloadmethod(libhandle::Ptr{Void},method_name::AbstractString)
  namefound = ""
  name_vars = ( method_name, uppercase(method_name), ucfirst(method_name) )
  trylist = tuple( 
    name_vars...,
    map(x->string(x,"_"),name_vars)...,
    map(x->string("_",x),name_vars)...)
  ptr = C_NULL
  for name in trylist
    ptr = Libdl.dlsym_e(libhandle,name)
    if (ptr != C_NULL) 
      namefound = name
      break
    end
  end
  ptr == C_NULL &&
    throw(ErrorException("Cannot find one of the methods $trylist"))
  return (ptr,namefound)
end

"""
       function loadODESolvers(extrapaths::Vector=AbstractString[],
                 loadlibnames::Tuple=() )
                
  tries to (dynamically) load the solvers.
  
  additional locations/paths to look at can be given as argument.
  
  If the 1st argument is an empty Vector, then the method tries to
  find the path of the ODEInterface module and (if successfull)
  uses this path as `extrapaths`.
  
  The 2nd argument is a `Tuple` with libnames of solvers to load. 
  If it is an empty tuple, then all known solvers will be tried.
  
  If an solver is already successfully loaded, then it will *not* be
  loaded again.
  
  returns `Dict` with informations about the loaded solvers (and errors).
  
  If a solver cannot be found (or needed methods inside a dynmic library
  cannot be found) then the errors are not propagated to the caller. The
  errors and expections are saved in the returned `Dict`. Why? Using this
  way, it is possible to see with one call (and try to load all solvers)
  which solvers are found.
  
  You can simply `dump` the values of the output dict to get a human-readable
  form of the result or call `help_solversupport()`.
  
       for k in keys(res); dump(res[k]); end
       ODEInterface.help_solversupport()
  """
function loadODESolvers(extrapaths::Vector=AbstractString[],
          loadlibnames::Tuple=() )
  path_sep = ""
  if isdefined(Base,:path_separator)
    path_sep = Base.path_separator
  end
  if isdefined(Base,:Filesystem) && isdefined(Base.Filesystem,:path_separator)
    path_sep = Base.Filesystem.path_separator
  end
  if isempty(extrapaths)
    try
      path_to_module = Base.find_in_path("ODEInterface")
      if path_to_module ≠ nothing
        if endswith(path_to_module,".jl")
          path_to_module = path_to_module[1:end-3]
        end
        if endswith(path_to_module,"ODEInterface")
          path_to_module = path_to_module[1:end-12]
        end
        if !endswith(path_to_module,path_sep)
          path_to_module = string(path_to_module,path_sep)
        end
        extrapaths = [ path_to_module ]
      end
    catch e
      # At least we tried
    end
  end
  for solver in solverInfo
    for variant in solver.variants
      libname = variant.libname
      if isempty(loadlibnames) || libname ∈ loadlibnames
        if !haskey(dlSolversInfo,libname) || 
           nothing ≠ dlSolversInfo[libname].error 
           
           libhandle = C_NULL; filepath =""; err = nothing
           mArray = Vector{MethodDLinfo}()
           try
             (libhandle,filepath) = trytoloadlib(libname,extrapaths)
             for generic_name in variant.methods
               methodsname_found = ""; method_ptr = C_NULL; merr = nothing
               try
                 (method_ptr, methodsname_found) = 
                   trytoloadmethod(libhandle, generic_name)
               catch e
                 merr = e
               end
               push!(mArray,MethodDLinfo(
                     generic_name, methodsname_found, method_ptr, merr))
             end
           catch e
             err = e
           end
           dlSolversInfo[libname] = SolverDLinfo(
             libname, filepath, libhandle, tuple(mArray...), err)
        end
      end
    end
  end
  
  return copy(dlSolversInfo)
end

"""
       function unloadODESolvers()

  unload all (loaded) solvers.
  """
function unloadODESolvers()
  all_keys = collect( keys(dlSolversInfo) )
  for key in all_keys
    libhandle = dlSolversInfo[key].libhandle
    libhandle ≠ C_NULL && Libdl.dlclose(libhandle)
    delete!(dlSolversInfo,key)
  end
  return nothing
end

"""
  return all method-pointers for a solver.
  
  tries to return all `method_ptr`s for all methods of a solver.
  This method checks if the `method_ptr`s are existent and different
  from `C_NULL`. If not then this method tries to load the 
  `dlname` ODE-Solver with the `loadODESolvers` method and checks again. 
  If even after this the `method_ptr`s are not found a exception is thrown.
  
  see `loadODESolvers`.
  """
function getAllMethodPtrs(dlname::AbstractString)
  load_tried = false
  ret = Vector{Ptr{Void}}()
  while true
    try
      empty!(ret)
      @assert dlSolversInfo[dlname].error == nothing
      for method in dlSolversInfo[dlname].methods
        @assert method.method_ptr ≠ C_NULL
        push!(ret,method.method_ptr)
      end
      break
    catch e
      if load_tried
        throw(SolverODEnotLoaded(string(
          "Cannot find method(s) for $(dlname)! ",
          "I've tried to loadODESolvers(), but it didn't work. ",
          "Please check ODEInterface.help_solversupport() and ",
          "call loadODESolvers and check also this output. ",
          "For further information see also ODEInterface.help_install.")))
      else
        loadODESolvers([],(dlname,)); load_tried = true
      end
    end
  end
  return ret
end

# vim:syn=julia:cc=79:fdm=indent:
