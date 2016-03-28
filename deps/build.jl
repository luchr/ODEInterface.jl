windows_flag = @windows? true : false
file_extension = @osx ? ".dylib" : @windows ? ".dll" : ".so" 
read_from_cmd = isdefined(Base,:readstring) ? Base.readstring : Base.readall
obj_files = []
verbose = false
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

function build_bvpsol(path::AbstractString)
  compile_gfortran(dir_of_src,"bvpsol")
  compile_gfortran(dir_of_src,"linalg_bvpsol")
  compile_gfortran(dir_of_src,"zibconst")
  compile_gfortran(dir_of_src,"ma28_bvpsol")
  link_gfortran(dir_of_src,
    ["bvpsol","linalg_bvpsol","zibconst","ma28_bvpsol"])
  println("\n\n!!! bvpsol: only non commercial use !!!")
  println("Please note: bvpsol's license only covers non commercial use!")
  println("see using ODEInterface; help_bvpsol_license() for the complete")
  println("license text.")
  return nothing
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
build_bvpsol(dir_of_src)

del_obj_files()

# vim:syn=julia:cc=79:fdm=indent:
