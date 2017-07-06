# backward compatible: no use of @static
macro check_for_windows()
  if isdefined(Base, :is_windows)
    return is_windows()? true : false 
  else
    return :( @windows ? true : false )
  end
end

macro check_for_apple()
  if isdefined(Base, :is_apple)
    return is_apple()? true : false
  else
    return :( @osx ? true : false )
  end
end


windows_flag = @check_for_windows()
apple_flag = @check_for_apple()
file_extension = apple_flag ? ".dylib" : windows_flag ? ".dll" : ".so" 
read_from_cmd = isdefined(Base,:readstring) ? Base.readstring : Base.readall
obj_files = []
verbose_key = "ODEINTERFACE_VERBOSE"
verbose = haskey(ENV, verbose_key) && 
          length(ENV[verbose_key])>0 ? true : false
gfortran = nothing

function search_prog(progname::AbstractString)
  output = ""
  env_key = string("ODEINTERFACE_",uppercase(progname))
  if haskey(ENV,env_key)
    output = ENV[env_key]
  else
    search_cmd = windows_flag ? `where "$progname"` : `which $progname`
    try
      output = rstrip(read_from_cmd(search_cmd))
    catch e
    end
  end
  return output
end

function compile_gfortran(path::AbstractString, basename::AbstractString)
  ffile = joinpath(path,string(basename,".f"))
  flags_i64 = [ "-fdefault-integer-8", "-fdefault-real-8",
                "-fdefault-double-8" ]
  flags_i32 = [ "-fdefault-real-8", "-fdefault-double-8" ]
  comp_flags = windows_flag ? [ "-c" ] : [ "-c", "-fPIC" ]

  ofile = joinpath(path,string(basename,".o"))
  if windows_flag
    cmd_i64 = `"$gfortran"  $comp_flags $flags_i64 -o "$ofile"  "$ffile"`
  else
    cmd_i64 = `"$gfortran"  $comp_flags $flags_i64 -o $ofile  $ffile` 
  end
  verbose && println(cmd_i64)
  run(cmd_i64)
  push!(obj_files,ofile)

  ofile = joinpath(path,string(basename,"_i32.o"))
  if windows_flag
    cmd_i32 = `"$gfortran"  $comp_flags $flags_i32 -o "$ofile"  "$ffile"`
  else
    cmd_i32 = `"$gfortran"  $comp_flags $flags_i32 -o $ofile  $ffile`
  end
  verbose && println(cmd_i32)
  run(cmd_i32)
  push!(obj_files,ofile)

  return nothing
end

function link_gfortran(path::AbstractString, basenames)
  link_flags = windows_flag ? [ "-shared" ] : [ "-shared", "-fPIC" ]
  i64_obj = map( name -> joinpath(path,string(name,".o")), basenames )
  sofile = joinpath(path,string(basenames[1],file_extension))
  if windows_flag
    cmd_i64 = `"$gfortran" $link_flags -o "$sofile" "$i64_obj"`
  else
    cmd_i64 = `"$gfortran" $link_flags -o $sofile $i64_obj`
  end
  verbose && println(cmd_i64)
  run(cmd_i64)

  i32_obj = map( name -> joinpath(path,string(name,"_i32.o")), basenames )
  sofile = joinpath(path,string(basenames[1],"_i32",file_extension))
  if windows_flag
    cmd_i32 = `"$gfortran" $link_flags -o "$sofile" "$i32_obj"`
  else
    cmd_i32 = `"$gfortran" $link_flags -o $sofile $i32_obj`
  end
  verbose && println(cmd_i32)
  run(cmd_i32)
  return nothing
end

function del_obj_files()
  for name in obj_files
    rm(name)
  end
  return nothing
end

function build_dopri(path::AbstractString)
  compile_gfortran(path,"dopri5")
  link_gfortran(path,["dopri5"])

  compile_gfortran(path,"dop853")
  link_gfortran(path,["dop853"])
  return nothing
end

function build_odex(path::AbstractString)
  compile_gfortran(path,"odex")
  link_gfortran(path,["odex"])
  return nothing
end

function compile_lapack(path::AbstractString)
  compile_gfortran(path,"dc_lapack")
  compile_gfortran(path,"lapack")
  compile_gfortran(path,"lapackc")
  return nothing
end

function build_radau(path::AbstractString)
  compile_gfortran(path,"radau5")
  link_gfortran(path,["radau5","dc_lapack","lapack","lapackc"])

  compile_gfortran(path,"radau")
  link_gfortran(path,["radau","dc_lapack","lapack","lapackc"])
  return nothing
end

function build_seulex(path::AbstractString)
  compile_gfortran(path,"seulex")
  link_gfortran(path,["seulex","dc_lapack","lapack","lapackc"])
  return nothing
end

function build_rodas(path::AbstractString)
  compile_gfortran(path,"rodas")
  link_gfortran(path,["rodas","dc_lapack","lapack","lapackc"])
  return nothing
end

function build_bvpsol(path::AbstractString)
  compile_gfortran(path,"bvpsol")
  compile_gfortran(path,"linalg_bvpsol")
  compile_gfortran(path,"zibconst")
  compile_gfortran(path,"ma28_bvpsol")
  link_gfortran(path,
    ["bvpsol","linalg_bvpsol","zibconst","ma28_bvpsol"])
  println("\n\n!!! bvpsol: only non commercial use !!!")
  println("Please note: bvpsol's license only covers non commercial use!")
  println("see using ODEInterface; help_bvpsol_license() for the complete")
  println("license text.")
  return nothing
end

function compile_slatec(path::AbstractString)
  compile_gfortran(path,"slatec")
  # compile_gfortran(path,"d1mach")
  # compile_gfortran(path,"daxpy")
  # compile_gfortran(path,"dcfod")
  # compile_gfortran(path,"ddes")
  # compile_gfortran(path,"dgbfa")
  # compile_gfortran(path,"dgbsl")
  # compile_gfortran(path,"dgefa")
  # compile_gfortran(path,"dgesl")
  # compile_gfortran(path,"dhstrt")
  # compile_gfortran(path,"dhvnrm")
  # compile_gfortran(path,"dintp")
  # compile_gfortran(path,"dintyd")
  # compile_gfortran(path,"dlsod")
  # compile_gfortran(path,"dpjac")
  # compile_gfortran(path,"dscal")
  # compile_gfortran(path,"dslvs")
  # compile_gfortran(path,"dsteps")
  # compile_gfortran(path,"dstod")
  # compile_gfortran(path,"dvnrms")
  # compile_gfortran(path,"fdump")
  # compile_gfortran(path,"i1mach")
  # compile_gfortran(path,"idamax")
  # compile_gfortran(path,"j4save")
  # compile_gfortran(path,"xercnt")
  # compile_gfortran(path,"xerhlt")
  # compile_gfortran(path,"xermsg")
  # compile_gfortran(path,"xerprn")
  # compile_gfortran(path,"xersve")
  # compile_gfortran(path,"xgetua")
  return nothing
end

function build_ddeabm(path::AbstractString)
  compile_gfortran(path,"ddeabm")

  link_gfortran(path,["ddeabm", "slatec"])
  # link_gfortran(path,["ddeabm",
  #   "d1mach", "ddes", "dhstrt", "dhvnrm", "dintp", "dsteps", "fdump",
  #   "i1mach", "j4save", "xercnt", "xerhlt", "xermsg", "xerprn", "xersve",
  #   "xgetua"])
end

function build_ddebdf(path::AbstractString)
  compile_gfortran(path,"ddebdf")

  link_gfortran(path,["ddebdf", "slatec"])
  # link_gfortran(path,["ddebdf",
  #   "d1mach", "daxpy", "dcfod", "dgbfa", "dgbsl", "dgefa", "dgesl", "dhstrt",
  #   "dhvnrm", "dintyd", "dlsod", "dpjac", "dslvs", "dscal", "dstod",
  #   "dvnrms", "fdump", "i1mach", "idamax", "j4save", "xercnt", "xerhlt",
  #   "xermsg", "xerprn", "xersve", "xgetua"])
end

# test for gfortran
gfortran = search_prog("gfortran")
if isempty(gfortran)
  error("Currently only gfortran is supported.")
end

dir_of_this_file = dirname(@__FILE__)
dir_of_src = normpath(joinpath(dir_of_this_file,"..","src"))
if !isdir(dir_of_src)
  error(string("Cannot find src directory. I tried: ",dir_of_src))
end

build_dopri(dir_of_src)
build_odex(dir_of_src)
compile_lapack(dir_of_src)
build_radau(dir_of_src)
build_seulex(dir_of_src)
build_rodas(dir_of_src)
build_bvpsol(dir_of_src)

compile_slatec(dir_of_src)
build_ddeabm(dir_of_src)
build_ddebdf(dir_of_src)

del_obj_files()

# vim:syn=julia:cc=79:fdm=indent:
