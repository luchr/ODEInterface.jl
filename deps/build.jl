unix_flag = @unix ? true : false

if !unix_flag
  error(string("Sorry. ",OS_NAME," is not yet supported."))
end

file_extension = @osx ? ".dylib" : ".so"

obj_files = []

function search_prog(progname::AbstractString)
  output = ""
  search_cmd = `which $progname`
  try
    output = readall(search_cmd)
  catch e
  end
  return output
end

function compile_gfortran(path::AbstractString, basename::AbstractString)
  ffile = joinpath(path,string(basename,".f"))
  flags_i64 = [ "-fdefault-integer-8", "-fdefault-real-8",
                "-fdefault-double-8" ]
  flags_i32 = [ "-fdefault-real-8", "-fdefault-double-8" ]

  ofile = joinpath(path,string(basename,".o"))
  cmd_i64=`gfortran  -c  -fPIC $flags_i64 -o $ofile  $ffile`
  run(cmd_i64)
  push!(obj_files,ofile)

  ofile = joinpath(path,string(basename,"_i32.o"))
  cmd_i32=`gfortran  -c  -fPIC $flags_i32 -o $ofile  $ffile`
  run(cmd_i32)
  push!(obj_files,ofile)

  return nothing
end

function link_gfortran(path::AbstractString, basenames)
  i64_obj = map( name -> joinpath(path,string(name,".o")), basenames )
  sofile = joinpath(path,string(basenames[1],file_extension))
  cmd_i64 = `gfortran -shared -fPIC -o $sofile $i64_obj`
  run(cmd_i64)

  i32_obj = map( name -> joinpath(path,string(name,"_i32.o")), basenames )
  sofile = joinpath(path,string(basenames[1],"_i32",file_extension))
  cmd_i32 = `gfortran -shared -fPIC -o $sofile $i32_obj`
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
  compile_gfortran(dir_of_src,"dopri5")
  link_gfortran(dir_of_src,["dopri5"])

  compile_gfortran(dir_of_src,"dop853")
  link_gfortran(dir_of_src,["dop853"])
  return nothing
end

function build_odex(path::AbstractString)
  compile_gfortran(dir_of_src,"odex")
  link_gfortran(dir_of_src,["odex"])
  return nothing
end

function compile_lapack(path::AbstractString)
  compile_gfortran(dir_of_src,"dc_lapack")
  compile_gfortran(dir_of_src,"lapack")
  compile_gfortran(dir_of_src,"lapackc")
  return nothing
end

function build_radau(path::AbstractString)
  compile_gfortran(dir_of_src,"radau5")
  link_gfortran(dir_of_src,["radau5","dc_lapack","lapack","lapackc"])

  compile_gfortran(dir_of_src,"radau")
  link_gfortran(dir_of_src,["radau","dc_lapack","lapack","lapackc"])
  return nothing
end

function build_seulex(path::AbstractString)
  compile_gfortran(dir_of_src,"seulex")
  link_gfortran(dir_of_src,["seulex","dc_lapack","lapack","lapackc"])
  return nothing
end

# test for gfortran
if isempty(search_prog("gfortran"))
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

del_obj_files()

# vim:syn=julia:cc=79:fdm=indent:
