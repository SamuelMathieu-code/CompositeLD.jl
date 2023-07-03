using SnpArrays
using Folds


function getStrongLD(ref_genotypes::SnpData, 
                     snps::AbstractVector{<:Tuple{Integer, Integer}}, 
                     r2_tresh::Float64 = 0.8,
                     window::Int64 = 500000,
                     snpDataFormated::Bool = true)
    if !snpDataFormated
        formatSnpData!(ref_genotypes)
    end

    snps_indx = Vector{Union{Int}}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        local j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        if j > lastindex(ref_genotypes.snp_info.chr_pos) || ref_genotypes.snp_info.chr_pos[j] != chr_pos_sing
            j = -1
        end
        snps_indx[i] = j
    end
    kept_input = snps_indx .> 0
    kept_indx = snps_indx[kept_input]

    res = Vector{Vector{Int}}(undef, length(kept_indx))
    for (i, j) in enumerate(kept_indx)
        idx_v_b = ref_genotypes.snp_info.chromosome .== snps[i][1] .& abs.(ref_genotypes.snp_info.position .- snps[i][2]) .< window
        indexes = ref_genotypes.idx[idx_v_b]
        kept = ld_r².(indexes, j, ref_genotypes.snparray) .> r2_tresh
        res[i] = indexes[kept]
    end
    res2 = Folds.reduce(vcat, res)

    return ref_genotypes.snp_info.chr_pos[res2]
    
end

function getStrongLD(ref_genotypes::AbstractVector{SnpData}, 
                     snps::AbstractVector{<:Tuple{Integer, Integer}}; 
                     r2_tresh::Float64 = 0.8,
                     window::Int64 = 500000,
                     snpDataFormated::Bool = true)
    if !snpDataFormated
    formatSnpData!.(ref_genotypes)
    end

    snps_indx = Vector{Int}(undef, size(snps, 1))
    snps_chr = Vector{Int}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        local j = searchsortedfirst(ref_genotypes[chr_pos_sing[1]].snp_info.chr_pos, chr_pos_sing)
        if j > lastindex(ref_genotypes[chr_pos_sing[1]].snp_info.chr_pos) || ref_genotypes[chr_pos_sing[1]].snp_info.chr_pos[j] != chr_pos_sing
            j = -1
        end
        snps_indx[i] = j
        snps_chr[i] = snps[i][1]
    end
    kept_input = snps_indx .> 0
    kept_snps = snps[kept_input]
    kept_indx = snps_indx[kept_input]
    kept_chr = snps_chr[kept_input]
  

    res = Vector{Vector{Int}}(undef, length(kept_indx))
    res_chrpos = Vector{Vector{Tuple{Int8, Int32}}}(undef, length(kept_indx))
    @threads for i in eachindex(kept_indx, kept_chr, res)
        idx_v_b = abs.(ref_genotypes[kept_chr[i]].snp_info.position .- kept_snps[i][2]) .< window
        indexes = ref_genotypes[kept_chr[i]].snp_info.idx[idx_v_b]
        corr = [ld_r²(s, kept_indx[i], ref_genotypes[kept_chr[i]].snparray) for s in indexes]
        kept = corr .> r2_tresh
        res[i] = indexes[kept]
    end
    
    @threads for i in eachindex(kept_chr, res, kept_indx)
        res_chrpos[i] = ref_genotypes[kept_chr[i]].snp_info.chr_pos[res[i]]
    end

    return kept_chr, Folds.reduce(vcat, res_chrpos, init = [])

end