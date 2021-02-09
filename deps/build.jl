try
  using Unicode
catch e
end

windows_flag = Sys.iswindows()
apple_flag = Sys.isapple()
file_extension = nothing
obj_files = []
verbose_key = "ODEINTERFACE_VERBOSE"
verbose = haskey(ENV, verbose_key) && 
          length(ENV[verbose_key])>0 ? true : false
build_script_key = "ODEINTERFACE_BUILD_SCRIPT"
build_script = get(ENV, build_script_key, nothing)
if build_script == ""
  build_script = nothing
end
build_script_io = nothing

gfortran = nothing

function adapt_to_os()
  global file_extension = apple_flag ? ".dylib" : windows_flag ? ".dll" : ".so" 
  return nothing
end

function search_prog(progname::AbstractString)
  output = ""
  env_key = string("ODEINTERFACE_",uppercase(progname))
  if haskey(ENV,env_key)
    output = ENV[env_key]
  else
    search_cmd = windows_flag ? `where "$progname"` : `which $progname`
    try
      output = rstrip(read(search_cmd, String))
    catch e
    end
  end
  return output
end

function build_script_write_cmd(cmd)
  if build_script_io != nothing
    cmd_to_str = string(cmd)
    write(build_script_io, cmd_to_str[2:end-1])
    write(build_script_io, "\n");
    return true
  else
    return false
  end
end

function build_script_write_headline(headline)
  if build_script_io != nothing
    write(build_script_io, "\n# ")
    write(build_script_io, headline)
    write(build_script_io, "\n")
    return true
  else
    return false
  end
end

function compile_gfortran(path::AbstractString, basename::AbstractString,
         options::Dict=Dict())
  fext = get(options, "file_extension", ".f")
  ffile = joinpath(path,string(basename, fext))
  flags_i64 = get(options, "flags_i64",
              [ "-fdefault-integer-8", "-fdefault-real-8",
                "-fdefault-double-8" ])
  append!(flags_i64, get(options, "add_flags_i64", []))
  flags_i32 = get(options, "flags_i32",
              [ "-fdefault-real-8", "-fdefault-double-8" ])
  append!(flags_i32, get(options, "add_flags_i32", []))
  comp_flags = windows_flag ? [ "-c" ] : [ "-c", "-fPIC" ]

  if get(options, "build_i64", true)
    ofile = joinpath(path,string(basename,".o"))
    if windows_flag
      cmd_i64 = `"$gfortran"  $comp_flags $flags_i64 -o "$ofile"  "$ffile"`
    else
      cmd_i64 = `"$gfortran"  $comp_flags $flags_i64 -o $ofile  $ffile` 
    end
    verbose && println(cmd_i64)
    if ! build_script_write_cmd(cmd_i64)
      run(cmd_i64)
    end
    push!(obj_files,ofile)
  end

  if get(options, "build_i32", true)
    ofile = joinpath(path,string(basename,"_i32.o"))
    if windows_flag
      cmd_i32 = `"$gfortran"  $comp_flags $flags_i32 -o "$ofile"  "$ffile"`
    else
      cmd_i32 = `"$gfortran"  $comp_flags $flags_i32 -o $ofile  $ffile`
    end
    verbose && println(cmd_i32)
    if ! build_script_write_cmd(cmd_i32)
      run(cmd_i32)
    end
    push!(obj_files,ofile)
  end

  return nothing
end

function link_gfortran(path::AbstractString, basenames, options::Dict=Dict())
  link_flags = windows_flag ? [ "-shared" ] : [ "-shared", "-fPIC" ]
  
  if get(options, "build_i64", true)
    i64_obj = map( name -> joinpath(path,string(name,".o")), basenames )
    sofile = joinpath(path,string(basenames[1],file_extension))
    if windows_flag
      cmd_i64 = `"$gfortran" $link_flags -o "$sofile" "$i64_obj"`
    else
      cmd_i64 = `"$gfortran" $link_flags -o $sofile $i64_obj`
    end
    verbose && println(cmd_i64)
    if ! build_script_write_cmd(cmd_i64)
      run(cmd_i64)
    end
  end

  if get(options, "build_i32", true)
    i32_obj = map( name -> joinpath(path,string(name,"_i32.o")), basenames )
    sofile = joinpath(path,string(basenames[1],"_i32",file_extension))
    if windows_flag
      cmd_i32 = `"$gfortran" $link_flags -o "$sofile" "$i32_obj"`
    else
      cmd_i32 = `"$gfortran" $link_flags -o $sofile $i32_obj`
    end
    verbose && println(cmd_i32)
    if ! build_script_write_cmd(cmd_i32)
      run(cmd_i32)
    end
  end
  return nothing
end

function del_obj_files()
  for name in obj_files
    rm(name; force=true)
  end
  return nothing
end

function build_dopri(path::AbstractString)
  build_script_write_headline("dopri")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"dopri5", options)
  link_gfortran(path,["dopri5"])

  compile_gfortran(path,"dop853", options)
  link_gfortran(path,["dop853"])
  return nothing
end

function build_odex(path::AbstractString)
  build_script_write_headline("odex")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"odex", options)
  link_gfortran(path,["odex"])
  return nothing
end

function compile_lapack(path::AbstractString)
  build_script_write_headline("lapack")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"dc_lapack", options)
  compile_gfortran(path,"lapack", options)
  compile_gfortran(path,"lapackc", options)
  return nothing
end

function build_radau(path::AbstractString)
  build_script_write_headline("radau")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"radau5", options)
  link_gfortran(path,["radau5","dc_lapack","lapack","lapackc"])

  compile_gfortran(path,"radau", options)
  link_gfortran(path,["radau","dc_lapack","lapack","lapackc"])
  return nothing
end

function build_seulex(path::AbstractString)
  build_script_write_headline("seulex")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"seulex", options)
  link_gfortran(path,["seulex","dc_lapack","lapack","lapackc"])
  return nothing
end

function build_rodas(path::AbstractString)
  build_script_write_headline("rodas")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"rodas", options)
  link_gfortran(path,["rodas","dc_lapack","lapack","lapackc"])
  return nothing
end

function build_bvpsol(path::AbstractString)
  build_script_write_headline("bvpsol")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"bvpsol", options)
  compile_gfortran(path,"linalg_bvpsol", options)
  compile_gfortran(path,"zibconst", options)
  compile_gfortran(path,"ma28_bvpsol", options)
  link_gfortran(path,
    ["bvpsol","linalg_bvpsol","zibconst","ma28_bvpsol"])
  println("\n\n!!! bvpsol: only non commercial use !!!")
  println("Please note: bvpsol's license only covers non commercial use!")
  println("see using ODEInterface; help_bvpsol_license() for the complete")
  println("license text.")
  return nothing
end

function compile_slatec(path::AbstractString)
  build_script_write_headline("slatec")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"slatec", options)
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
  build_script_write_headline("ddeabm")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"ddeabm", options)

  link_gfortran(path,["ddeabm", "slatec"])
  # link_gfortran(path,["ddeabm",
  #   "d1mach", "ddes", "dhstrt", "dhvnrm", "dintp", "dsteps", "fdump",
  #   "i1mach", "j4save", "xercnt", "xerhlt", "xermsg", "xerprn", "xersve",
  #   "xgetua"])
end

function build_ddebdf(path::AbstractString)
  build_script_write_headline("ddebdf")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path,"ddebdf", options)

  link_gfortran(path,["ddebdf", "slatec"])
  # link_gfortran(path,["ddebdf",
  #   "d1mach", "daxpy", "dcfod", "dgbfa", "dgbsl", "dgefa", "dgesl", "dhstrt",
  #   "dhvnrm", "dintyd", "dlsod", "dpjac", "dslvs", "dscal", "dstod",
  #   "dvnrms", "fdump", "i1mach", "idamax", "j4save", "xercnt", "xerhlt",
  #   "xermsg", "xerprn", "xersve", "xgetua"])
end

function build_colnew(path::AbstractString)
  build_script_write_headline("colnew")
  options = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "add_flags_i32" => ["-w", "-std=legacy"],
  )
  compile_gfortran(path, "colnew", options)
  link_gfortran(path, ["colnew",])
  return nothing
end

function build_bvpm2(path::AbstractString)
  build_script_write_headline("bvpm2")
  opt = Dict(
    "add_flags_i64" => ["-w", "-std=legacy"],
    "build_i32"      => false)
  proxy_options = Dict(
    "file_extension" => ".f90",
    "build_i32"      => false,
    "add_flags_i64"  => [
      "-Wall", "-Wextra", "-Wimplicit-interface", "-std=f2008ts"],
  )
  compile_gfortran(path, "bvp_la-2", opt)
  compile_gfortran(path, "bvp_m-2", Dict(
    "file_extension" => ".f90",
    "add_flags_i64"  => [ "-std=f2008" ],
    "build_i32"      => false,))
  compile_gfortran(path, "bvp_m_proxy", proxy_options)
  link_gfortran(path, ["bvp_m_proxy", "bvp_m-2", "bvp_la-2"], opt)
end

function build_all(dir_of_src::AbstractString)
  build_dopri(dir_of_src)
  build_odex(dir_of_src)
  compile_lapack(dir_of_src)
  build_radau(dir_of_src)
  build_seulex(dir_of_src)
  build_rodas(dir_of_src)

  compile_slatec(dir_of_src)
  build_ddeabm(dir_of_src)
  build_ddebdf(dir_of_src)

  build_bvpsol(dir_of_src)
  build_colnew(dir_of_src)
  build_bvpm2(dir_of_src)

  del_obj_files()
  return nothing
end

function build_with_gfortran()
  global gfortran = search_prog("gfortran")
  if isempty(gfortran)
    error("Currently only gfortran is supported.")
  end
  dir_of_this_file = dirname(@__FILE__)
  dir_of_src = normpath(joinpath(dir_of_this_file, "..", "src"))
  if !isdir(dir_of_src)
    error(string("Cannot find src directory. I tried: ", dir_of_src))
  end
  build_all(dir_of_src)
  return nothing
end

function create_build_script()
  global build_script_io = open(build_script, "w")
  global gfortran = "gfortran"
  global windows_flag
  global apple_flag
  dir_of_src = "./"

  try
    write(build_script_io, raw"""if [[ $target == *mingw* ]] ; then""" * "\n")
    (windows_flag, apple_flag) = (true, false); adapt_to_os()
    build_all(dir_of_src)
    write(build_script_io, "\n\n\n")

    write(build_script_io, raw"""elif [[ $target == *apple* ]] ; then""" * "\n")
    (windows_flag, apple_flag) = (false, true); adapt_to_os()
    build_all(dir_of_src)
    write(build_script_io, "\n\n\n")

    write(build_script_io, "else\n")
    (windows_flag, apple_flag) = (false, false); adapt_to_os()
    build_all(dir_of_src)
    write(build_script_io, "\n\n\n")

    write(build_script_io, "fi\n\n")
  finally
    if build_script_io != nothing
      close(build_script_io)
    end
  end
end

if build_script != nothing
  create_build_script()
else
  adapt_to_os()
  build_with_gfortran()
end

# vim:syn=julia:cc=79:fdm=indent:
