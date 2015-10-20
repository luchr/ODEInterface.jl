# ODE-Options for ODEInterface

import Base: show

"""macro for importing OptionsODE and option handling."""
macro import_options()
  :(
    using ODEInterface: OptionsODE, getOption, setOption!, setOptions!, 
                        copyOptions!
  )
end

"""
  Ancestor for all types storing options for ODE-solvers.
  
  ODE-solvers often have serveral parameters for fine-tuning them.
  In this ODEInterface this parameters are called 'options' and 
  they are stored in key/value paris. For the key a 
  `ASCIIString` is used. The value can be `Any`-thing.
  The key is often called the option-name.
  
  All types for this purpose have this abstract type as super-type.
  
  Required fields are: `name`, `lastchanged`, `options`
  """
abstract AbstractOptionsODE <: Any

"""
  Stores options for ODE-Solver(s) together with a name.
  Additionally the time of the last change is saved.
  
  Options can be set at construction time, e.g.
  
       opt=OptionsODE("test",
                      "loglevel" => ODEInterface.LOG_ALL,
                      "logio"    => STDERR)
  
  or later. For changing single options 
  
       oldValue = setOption!(opt,"myopt","new value")
       oldValue = setOption!(opt,"myopt" => "new value")
  
  and for changing many options at once:
  
       oldValues = setOption!(opt,
                   "myopt" => "new value",
                   "oldopt" => 56)
  
  see also: `setOption!`, `setOptions!`
  """
type OptionsODE <: AbstractOptionsODE
  name        :: AbstractString
  lastchanged :: DateTime
  options     :: Dict{ASCIIString,Any}

  function OptionsODE(name::AbstractString="")
    obj = new(name,now(),Dict{ASCIIString,Any}())
    return obj
  end
end

function OptionsODE(name::AbstractString,copyOptionsFrom::AbstractOptionsODE) 
  opt = OptionsODE(name)
  copyOptions!(opt,copyOptionsFrom)
  return opt
end

function OptionsODE(copyOptionsFrom::AbstractOptionsODE) 
  return OptionsODE("",copyOptionsFrom)
end

function OptionsODE(name::AbstractString,pairs::Pair...)
  opt = OptionsODE(name)
  setOptions!(opt,pairs...)
  return opt
end

function OptionsODE(pairs::Pair...)
  return OptionsODE("",pairs...)
end

"""
     function getOption(opt::AbstractOptionsODE,name::ASCIIString,
                        default::Any=nothing)

  get saved value of option with given `name` or `default` 
  if option is unknown.
  """
function getOption(opt::AbstractOptionsODE,name::ASCIIString,
                   default::Any=nothing)
  return haskey(opt.options,name)?opt.options[name]:default
end

"""
     function setOption!(opt::AbstractOptionsODE,name::ASCIIString,value::Any)

  set ODE-Option with given `name` and return old value 
  (or `nothing` if there was no old value).
  """
function setOption!(opt::AbstractOptionsODE,name::ASCIIString,value::Any)
  oldValue = getOption(opt,name)
  opt.options[name]=value
  opt.lastchanged=now()
  return oldValue
end

"""
     function setOption!(opt::AbstractOptionsODE,pair::Pair)

  set ODE-Option with given (`name`,`value`) pair and return old value 
  (or `nothing` if there was no old value).
  """
function setOption!(opt::AbstractOptionsODE,pair::Pair)
  return setOption!(opt,pair.first,pair.second)
end

"""
     function setOptions!(opt::AbstractOptionsODE,pairs::Pair...)

  set many ODE-Options and return an array with the old option values.
  """
function setOptions!(opt::AbstractOptionsODE,pairs::Pair...)
  oldValues=Any[]
  for (name,value) in pairs
    push!(oldValues, setOption!(opt,name,value))
  end
  return oldValues
end

"""
     function copyOptions!(dest::AbstractOptionsODE,source::AbstractOptionsODE)

  copy all options from other ODE-Option object.
  """
function copyOptions!(dest::AbstractOptionsODE,source::AbstractOptionsODE)
  merge!(dest.options,source.options)
  dest.lastchanged=now()
  return dest
end

function show(io::IO, opt::AbstractOptionsODE)
  print(io,typeof(opt)," ")
  isempty(opt.name) || print(io,"'",opt.name,"' ")
  len=length(opt.options)
  print(io,"with ",len," option",len!=1?"s":"",len>0?":":"."); println(io)
  if len>0
    maxLen=2+max( 0,map(length,keys(opt.options))... )
    for key in sort(collect(keys(opt.options)))
      print(io,lpad(key,maxLen),": ")
      show(io,opt.options[key]); println(io)
    end
  end
  print(io,"lasted changed "); show(io,opt.lastchanged); 
end

# vim:syn=julia:cc=79:fdm=indent:
