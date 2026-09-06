using Documenter
using METTSLibrary

makedocs(
    sitename = "METTSLibrary.jl",
    modules  = [METTSLibrary],
    format   = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        edit_link  = "main",
        canonical  = "https://awietek.github.io/METTSLibrary.jl",
        assets     = String[],
    ),
    pages = [
        "Overview"                => "index.md",
        "Getting started"         => "getting_started.md",
        "Using samples"           => "using.md",
        "Adding samples"          => "adding.md",
        "Lattice files"           => "lattices.md",
        "Converting legacy output" => "converting.md",
        "Sharing"                 => "sharing.md",
        "Hosting and copies"      => "hosting.md",
        "File format"             => "schema.md",
        "API"                     => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
    repo     = Remotes.GitHub("awietek", "METTSLibrary.jl"),
)

deploydocs(repo = "github.com/awietek/METTSLibrary.jl.git", devbranch = "main")
