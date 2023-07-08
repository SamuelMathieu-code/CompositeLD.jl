using SnpArrays
using Base.Threads

"""
Implementation of the clumping algorithm prioritising first snps in given Vector and formated Genotypes SnpData (see formatSnpData!)
    returns a vector of booean indication if each given snp is kept
"""
function clump(ref_genotypes::SnpData, 
               snps::AbstractVector{<:Tuple{Integer, Integer}}; 
               r2_tresh::AbstractFloat = 0.1,
               formated = false
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    indx_v_b = Vector{Bool}(undef, size(snps, 1))
    for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        @inbounds begin 
            indx_v_b[i] = firstindex(ref_genotypes.snp_info.chr_pos) ≤ j ≤ lastindex(ref_genotypes.snp_info.chr_pos) && 
                          ref_genotypes.snp_info.chr_pos[j] == chr_pos_sing
            if indx_v_b[i]
                snps_indx[i] = ref_genotypes.snp_info.idx[j]
            end
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
               formated = false
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes, :snpid)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    indx_v_b = Vector{Bool}(undef, size(snps, 1))
    for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.snpid, chr_pos_sing)
        @inbounds begin 
            indx_v_b[i] = firstindex(ref_genotypes.snp_info.snpid) ≤ j ≤ lastindex(ref_genotypes.snp_info.snpid) && 
                          ref_genotypes.snp_info.snpid[j] == chr_pos_sing
            if indx_v_b[i]
                snps_indx[i] = ref_genotypes.snp_info.idx[j]
            end
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
"""
function tclump(ref_genotypes::SnpData, 
               snps::AbstractVector{<:Tuple{Integer, Integer}}; 
               r2_tresh::AbstractFloat = 0.1,
               formated = false
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1)) # indices 
    indx_v_b = Vector{Bool}(undef, size(snps, 1)) # found or not
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        @inbounds begin 
            indx_v_b[i] = firstindex(ref_genotypes.snp_info.chr_pos) ≤ j ≤ lastindex(ref_genotypes.snp_info.chr_pos) && 
                          ref_genotypes.snp_info.chr_pos[j] == chr_pos_sing
            if indx_v_b[i]
                snps_indx[i] = ref_genotypes.snp_info.idx[j]
            end
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
               formated::Bool = false
               )::Vector{Bool}
    
    if !formated
        formatSnpData!(ref_genotypes, :snpid)
    end

    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    indx_v_b = Vector{Bool}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.snpid, chr_pos_sing)
        @inbounds begin 
            indx_v_b[i] = firstindex(ref_genotypes.snp_info.snpid) ≤ j ≤ lastindex(ref_genotypes.snp_info.snpid) && 
                          ref_genotypes.snp_info.snpid[j] == chr_pos_sing
            if indx_v_b[i]
                snps_indx[i] = ref_genotypes.snp_info.idx[j]
            end
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