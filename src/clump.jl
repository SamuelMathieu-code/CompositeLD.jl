using SnpArrays
using Base.Threads

"""
Implementation of the clumping algoritm prioritising first snps in given Vector and formated Genotypes SnpData (see formatSnpData!)
    returns a vector of booean indication if each given snp is kept
"""
function clump(ref_genotypes::SnpData, 
               snps::Union{AbstractVector{<:Tuple{Integer, Integer}}, AbstractVector{<:AbstractString}}; 
               r2_tresh::Float64 = 0.1,
               formated = false
               )::Vector{Bool}
    
    r2_mat, indx_v_b = getLDmat(ref_genotypes, snps, formated)
    idx_on_mat = accumulate(+, indx_v_b)
    
    for i in 1:(lastindex(snps)-1)
        if indx_v_b[i]
            @threads for j in (i+1):lastindex(snps)
                if r2_mat[idx_on_mat[i], idx_on_mat[j]] > r2_tresh
                    indx_v_b[j] = false
                end
            end
        end
    end

    return indx_v_b
end
