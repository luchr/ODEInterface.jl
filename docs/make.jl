using Pkg

Pkg.develop(path = joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Documenter, ODEInterface

makedocs(;
    modules = [ODEInterface],
    sitename = "ODEInterface.jl",
    repo = "github.com/luchr/ODEInterface.jl",
    checkdocs = :public,
    pages = ["Public API" => "index.md"],
)
