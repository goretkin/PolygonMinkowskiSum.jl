using Documenter, PolygonMinkowskiSum

makedocs(
    modules = [PolygonMinkowskiSum],
    format = Documenter.HTML(),
    checkdocs = :exports,
    sitename = "PolygonMinkowskiSum.jl",
    pages = Any["index.md"]
)

deploydocs(
    repo = "github.com/goretkin/PolygonMinkowskiSum.jl.git",
)
