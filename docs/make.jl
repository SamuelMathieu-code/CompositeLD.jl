push!(LOAD_PATH, "../src/")
using Documenter, CompositeLD


makedocs(sitename="CompositeLD.jl", format = Documenter.HTML(prettyurls = false))