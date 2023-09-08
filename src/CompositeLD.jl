module CompositeLD


export maf, ld_r2, getLDmat, formatSnpData!
export clump, tclump
export getStrongLD


include("ld.jl")
include("clump.jl")
include("getStrongLD.jl")

end
