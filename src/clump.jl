using SnpArrays
using Base.Threads

"""
Implementation of the clumping algorithm prioritising first snps in given Vector and formated Genotypes SnpData (see formatSnpData!)
    returns a vector of booean indication if each given snp is kept

**arguments :**

`ref_genotypes::SnpData` : reference genotypes (formated by [`formatSnpData!`](@ref) or not)
`snps::Union{AbstractVector{<:Tuple{integer, Integer}}, AbstractVector{<:AbstractString}}` : snps to clump


**options :** 

`formated::Bool` : indicates if ref SnpData is already formated according to :chr_pos (chr, pos) or :snpid (id as string)\\
`r2_tresh::AbstractFloat` : minimal r² for 2 snps to be considered strongly correlated.


## Examples :

```julia-repl
julia> ref = SnpData(datadir("some/data"));

julia> kept_v_b::Vector{Bool} = clump([(1, 123), (1, 456), (1, 789)], ref)
3-element Vector{Int64}:
 1
 0
 1

julia> formatSnpData!(ref, :snpid);

julia> kept_v_b::Vector{Bool} = clump(["rs123", "rs456", "rs789"], ref, formated = true)
3-element Vector{Int64}:
 1
 0
 1
```

If formatSnpData has already been called on good snp info type (`:chr_pos` or `:snpid`), `formated = true` option does not verify or modify formating.
See [`formatSnpData!`](@ref).
"""
function clump(ref_genotypes::SnpData, 
               snps::AbstractVector{<:Tuple{Integer, Integer}}; 
               r2_tresh::AbstractFloat = 0.1,
               formated = false,
               min_maf::Real = 0
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    indx_v_b = Vector{Bool}(undef, size(snps, 1))
    for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsorted(ref_genotypes.snp_info.chr_pos, chr_pos_sing) # range valid if of length 1 (> 1 => triallelic, < 1 => not found)
        if length(j) == 1
            snps_indx[i] = ref_genotypes.snp_info.idx[j[]]

            # if maf insufficient : discard variant
            if maf(ref_genotypes.snparray[:, snps_indx[i]]) > min_maf
                indx_v_b[i] = true
            else
                indx_v_b[i] = false
            end
        else
            # if variant not biallelic : discard
            indx_v_b[i] = false
        end
    end
    
    for i in 1:(lastindex(snps)-1)
        if indx_v_b[i]
            for j in i+1:lastindex(snps)
                if indx_v_b[j]
                    s1 = @view ref_genotypes.snparray[:, snps_indx[i]]
                    s2 = @view ref_genotypes.snparray[:, snps_indx[j]]
                    if _ld_r²_(s1, s2) > r2_tresh
                        indx_v_b[j] = false
                    end
                end
            end
        end
    end

    return indx_v_b
end


function clump(ref_genotypes::SnpData, 
               snps::AbstractVector{<:AbstractString}; 
               r2_tresh::AbstractFloat = 0.1,
               formated = false,
               min_maf::Real = 0
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes, :snpid)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    indx_v_b = Vector{Bool}(undef, size(snps, 1))
    for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsorted(ref_genotypes.snp_info.snpid, chr_pos_sing)
        if length(j) == 1
            snps_indx[i] = ref_genotypes.snp_info.idx[j[]]

            # if maf insufficient : discard variant
            if maf(ref_genotypes.snparray[:, snps_indx[i]]) > min_maf
                indx_v_b[i] = true
            else
                indx_v_b[i] = false
            end
        else
            # if variant not biallelic : discard
            indx_v_b[i] = false
        end
    end
    
    for i in 1:(lastindex(snps)-1)
        if indx_v_b[i]
            for j in i+1:lastindex(snps)
                if indx_v_b[j]
                    s1 = @view ref_genotypes.snparray[:, snps_indx[i]]
                    s2 = @view ref_genotypes.snparray[:, snps_indx[j]]
                    if _ld_r²_(s1, s2) > r2_tresh
                        indx_v_b[j] = false
                    end
                end
            end
        end
    end

    return indx_v_b
end

"""
threaded implementation of the clumping algorithm prioritising first snps in given Vector and formated Genotypes SnpData (see formatSnpData!)
    returns a vector of booean indication if each given snp is kept

**arguments :**

`ref_genotypes::SnpData` : reference genotypes (formated by [`formatSnpData!`](@ref) or not)
`snps::Union{AbstractVector{<:Tuple{integer, Integer}}, AbstractVector{<:AbstractString}}` : snps to clump


**options :** 

`formated::Bool` : indicates if ref SnpData is already formated according to :chr_pos (chr, pos) or :snpid (id as string)\\
`r2_tresh::AbstractFloat` : minimal r² for 2 snps to be considered strongly correlated.
        
## Examples :

```julia
ref = SnpData(datadir("some/data"))

kept_v_b::Vector{Bool} = tclump([(1, 123), (1, 456), (1, 789)], ref)

formatSnpData!(ref, :snpid)

kept_v_b::Vector{Bool} = tclump(["rs123", "rs456", "rs789"], ref, formated = true)
```

If formatSnpData has already been called on good snp info type (`:chr_pos` or `:snpid`), `formated = true` option does not verify or modify formating.
See [`formatSnpData!`](@ref).
"""
function tclump(ref_genotypes::SnpData, 
               snps::AbstractVector{<:Tuple{Integer, Integer}}; 
               r2_tresh::AbstractFloat = 0.1,
               formated = false,
               min_maf::Real = 0
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1)) # indices 
    indx_v_b = Vector{Bool}(undef, size(snps, 1)) # found or not
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsorted(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        @inbounds if length(j) == 1
            snps_indx[i] = ref_genotypes.snp_info.idx[j[]]

            # if maf insufficient : discard variant
            if maf(ref_genotypes.snparray[:, snps_indx[i]]) > min_maf
                indx_v_b[i] = true
            else
                indx_v_b[i] = false
            end
        else
            # if variant not biallelic : discard
            @inbounds indx_v_b[i] = false
        end
    end

    
    for i in 1:(lastindex(snps)-1)
        if indx_v_b[i]
            @threads for j in i+1:lastindex(snps)
                if indx_v_b[j]
                    s1 = @view ref_genotypes.snparray[:, snps_indx[i]]
                    s2 = @view ref_genotypes.snparray[:, snps_indx[j]]
                    if _ld_r²_(s1, s2) > r2_tresh
                            indx_v_b[j] = false
                    end
                end
            end
        end
    end

    return indx_v_b
end


function tclump(ref_genotypes::SnpData, 
               snps::AbstractVector{<:AbstractString}; 
               r2_tresh::AbstractFloat = 0.1,
               formated::Bool = false,
               min_maf::Real = 0
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes, :snpid)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    indx_v_b = Vector{Bool}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsorted(ref_genotypes.snp_info.snpid, chr_pos_sing)
        @inbounds if length(j) == 1
            snps_indx[i] = ref_genotypes.snp_info.idx[j[]]

            # if maf insufficient : discard variant
            if maf(ref_genotypes.snparray[:, snps_indx[i]]) > min_maf
                indx_v_b[i] = true
            else
                indx_v_b[i] = false
            end
        else
            # if variant not biallelic : discard
            indx_v_b[i] = false
        end
    end

    
    for i in 1:(lastindex(snps)-1)
        if indx_v_b[i]
            @threads for j in i+1:lastindex(snps)
                if indx_v_b[j]
                    s1 = @view ref_genotypes.snparray[:, snps_indx[i]]
                    s2 = @view ref_genotypes.snparray[:, snps_indx[j]]
                    if _ld_r²_(s1, s2) > r2_tresh
                            indx_v_b[j] = false
                    end
                end
            end
        end
    end

    return indx_v_b
end