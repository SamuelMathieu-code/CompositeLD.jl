using SnpArrays
using Folds

"""
Function : Get strongly correlated variants from list.
"""
function getStrongLD(ref_genotypes::SnpData, 
                     snps::AbstractVector{<:Tuple{Integer, Integer}}; 
                     r2_tresh::Float64 = 0.8,
                     window::Int64 = 500000,
                     formated::Bool = false)
    if !formated
        formatSnpData!(ref_genotypes)
    end

    snps_indx = Vector{Int}(undef, size(snps, 1))
    kept_input = Vector{Bool}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        @inbounds begin 
            kept_input[i] = firstindex(ref_genotypes.snp_info.chr_pos) ≤ j ≤ lastindex(ref_genotypes.snp_info.chr_pos) && 
                          ref_genotypes.snp_info.chr_pos[j] == chr_pos_sing
            snps_indx[i] = j
        end
    end
    snps_indx = snps_indx[kept_input]
    kept_idx = ref_genotypes.snp_info.idx[snps_indx]
    kept_chr = ref_genotypes.snp_info.chromosome[snps_indx]
    kept_pos = ref_genotypes.snp_info.position[snps_indx]

    res = Vector{Vector{Int}}(undef, length(kept_idx))
    @threads for i in eachindex(kept_idx)
        t1 = ref_genotypes.snp_info.chromosome .== kept_chr[i]
        t2 = abs.(ref_genotypes.snp_info.position .- kept_pos[i]) .< window
        idx_v_b = t1 .& t2
        indexes = ref_genotypes.snp_info.idx[idx_v_b]
        kept = findall(idx_v_b)
        corr = [ld_r²(s, kept_idx[i], ref_genotypes.snparray) for s in indexes]
        res[i] = kept[corr .> r2_tresh]
    end
    res2 = Folds.reduce(vcat, res, init = [])

    return ref_genotypes.snp_info.chr_pos[res2]
    
end


function getStrongLD(ref_genotypes::AbstractVector{SnpData}, 
                     snps::AbstractVector{<:Tuple{Integer, Integer}}; 
                     r2_tresh::Float64 = 0.8,
                     window::Int64 = 500000,
                     formated::Bool = false)
    
    a = Vector{Vector{Tuple{Int8, Int}}}(repeat([[]], 22))
    for c in snps
        push!(a[c[1]], c)
    end
    r = Vector{Tuple{Int8, Int64}}([])
    for i in 1:22
        append!(r, getStrongLD(ref_genotypes[i], a[i], r2_tresh = r2_tresh, window = window, formated = formated))
    end
    return r
end


function getStrongLD(ref_genotypes::SnpData, 
                     snps::AbstractVector{<:AbstractString};
                     r2_tresh::Float64 = 0.8,
                     window::Int64 = 500000,
                     formated::Bool = false)
    if !formated
        formatSnpData!(ref_genotypes, :snpid)
    end

    snps_indx = Vector{Int}(undef, size(snps, 1))
    kept_input = Vector{Bool}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.snpid, chr_pos_sing)
        @inbounds begin 
            kept_input[i] = firstindex(ref_genotypes.snp_info.snpid) ≤ j ≤ lastindex(ref_genotypes.snp_info.snpid) && 
                          ref_genotypes.snp_info.snpid[j] == chr_pos_sing
            snps_indx[i] = j
        end
    end
    snps_indx = snps_indx[kept_input]
    kept_idx = ref_genotypes.snp_info.idx[snps_indx]
    kept_chr = ref_genotypes.snp_info.chromosome[snps_indx]
    kept_pos = ref_genotypes.snp_info.position[snps_indx]


    res = Vector{Vector{Int}}(undef, length(kept_idx))
    @threads for i in eachindex(kept_idx)
        chr, pos = kept_chr[i], kept_pos[i]
        t1 = ref_genotypes.snp_info.chromosome .== chr
        t2 = abs.(ref_genotypes.snp_info.position .- pos) .< window
        idx_v_b = t1 .& t2
        indexes = ref_genotypes.snp_info.idx[idx_v_b]
        kept = findall(idx_v_b)
        corr = [ld_r²(s, kept_idx[i], ref_genotypes.snparray) for s in indexes]
        res[i] = kept[corr .> r2_tresh]
    end
    res2 = Folds.reduce(vcat, res, init = [])

    return ref_genotypes.snp_info.snpid[res2]
    
end


function getStrongLD(ref_genotypes::AbstractVector{SnpData}, 
    snps::AbstractVector{<:AbstractString}; 
    r2_tresh::Float64 = 0.8,
    window::Int64 = 500000,
    formated::Bool = false)

    res = Vector{String}([])
    for i in eachindex(ref_genotypes)
        append!(res, getStrongLD(ref_genotypes[i], snps, r2_tresh = r2_tresh, window = window, formated = formated))
    end

    return unique(res)
end

