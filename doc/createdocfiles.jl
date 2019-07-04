# This is a "hack" to automatically produce some doc files for github 
# using julia's doc-strings.
using Markdown

using ODEInterface
@ODEInterface.import_huge

const NL = "\n"

# takebuf_string deprecation: 358c4419
if VERSION < v"0.6.0-dev+1254"
  buf2str(buf) = takebuf_string(buf)
else
  buf2str(buf) = String(take!(buf))
end

"""
  escapes characters with "&#...;" HTML-notation.
  """
function escapeChars(s::AbstractString,toreplace=r"([^a-zA-Z0-9 \n])")
  return replace(s,toreplace => c -> string("&#",Int(c[1]),";"))
end

formatTable_new_row_for_nl = false

"""
  try to convert the Unicode doc-tables to HTML-tables.
  """
function formatTable(io,s::AbstractString)
  table_head = true
  std_line_content = r"^\s*║([^║]+)║\s*$"m
  head_stop  = r"^\s*╠[═╪]+╣\s*$"m
  row_stop = r"^\s*(╟[─┼]+╢|╚[═╧]+╝)\s*"m

  write(io,"<table>",NL)

  lines = split(s,"\n")
  columns = Vector{IOBuffer}(undef, length(split(lines[1],"╤")))
  for k in 1:length(columns)
    columns[k] = IOBuffer()
  end
  write_row = false
  for line in lines[2:end]
    mo = match(std_line_content,line)
    if mo !== nothing
      contents = split(mo[1],"│")
      @assert length(contents)==length(columns)
      for (col,content) in zip(columns,contents)
        data = rstrip(content)
        if length(data)>0
          write(col,escapeChars(data),NL)
        end
      end
    end
    mo = match(head_stop,line)
    if mo !== nothing
      @assert table_head
      write_row = true
    end
    mo = match(row_stop,line)
    if mo !== nothing
      write_row = true
    end
    if !table_head && formatTable_new_row_for_nl
      write_row = true
    end
    if write_row
      write(io,"<tr>")
      for col in columns
        write(io,table_head ? "<th>" : "<td>",
          "<pre>",buf2str(col),"</pre>",
          table_head ? "</th>" : "</td>",NL)
      end
      write(io,"</tr>",NL)
      write_row = false; table_head = false
    end
  end
  write(io,"</table>",NL)
  return nothing
end

function formatMDelement(io,e::Markdown.MD)
  for element in e.content
    formatMDelement(io,element)
    write(io,"\n")
  end
  return nothing
end

function formatMDelement(io,code::Markdown.Code)
  table_check = r"^\s*[╔══╤╗]+\s*$"m
  if occursin(table_check,code.code)
    formatTable(io,code.code)
  else
    Markdown.plain(io,code)
  end
  return nothing
end

function formatMDelement(io,rest::Any)
  Markdown.plain(io,rest)
  return nothing
end

function introHeader(io)
  write(io,"[ This file was auto-generated from the module's documentation ",
           "included in the doc-strings. Use julia's help system to get ",
           "these informations in a nicer output format. ]",NL,NL)
  return nothing
end

function docSolverOptions(filename)
  io = open(filename,"w")
  introHeader(io)
  namewomodule = r"\.([^.]+)$"

  for solver in (dopri5,dop853,odex,seulex,rodas,ddeabm,ddebdf,)
    solvername = string(solver)
    mo = match(namewomodule,solvername)
    if mo !== nothing
      solvername = mo[1]
    end
    write(io,"# ",solvername,NL,NL)
    formatMDelement(io,Base.Docs.doc(solver))
  end
  # radau5 and radau share same documentation
  write(io,"# radau and radau5",NL,NL)
  formatMDelement(io,Base.Docs.doc(radau))
  
  # bvpsol
  write(io,"# bvpsol",NL,NL)
  formatMDelement(io,Base.Docs.doc(bvpsol))

  # colnew
  write(io,"# colnew",NL,NL)
  formatMDelement(io,Base.Docs.doc(colnew))

  # bvpm2
  write(io,"# bvpm2",NL,NL)
  formatMDelement(io,Base.Docs.doc(Bvpm2))
  formatMDelement(io,Base.Docs.doc(bvpm2_init))
  formatMDelement(io,Base.Docs.doc(bvpm2_solve))

  close(io)
  return nothing
end

function docstringToFile(filename,docobjs)
  io = open(filename,"w")
  introHeader(io)
  for docobj in docobjs
    md_elem = isa(docobj,Markdown.MD) ? docobj : Base.Docs.doc(docobj)
    formatMDelement(io,md_elem)
    write(io,NL)
  end
  close(io)
  return nothing
end


docSolverOptions("./SolverOptions.md")

formatTable_new_row_for_nl = true
docstringToFile("./OptionOverview.md",[
   ODEInterface.help_options, ODEInterface.help_options()])
formatTable_new_row_for_nl = false

docstringToFile("./OutputFunction.md",[ODEInterface.help_outputfcn])
docstringToFile("./SpecialStructure.md",[ODEInterface.help_specialstructure])
docstringToFile("./CallSolvers.md",[
   ODEInterface.help_callsolvers,ODEInterface.odecall])

# vim:syn=julia:cc=79:fdm=indent:
