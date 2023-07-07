using SnpArrays
import IterTools.subsets
using LinearAlgebra
using Base.Threads

########## Genotypes in Plink .bed format #########

macro HOMO_M()  # AA
    return :(0x00)
end
macro MISSING()
    return :(0x01)
end
macro HETER()   # aA / Aa
    return :(0x02)
end
macro HOMO_m()  # aa
    return :(0x03)
end

############ LD calculations ############

"""
LD r² composite of pair of SNPs given two vectors of genotypes
"""
function ld_r²(snp1::Vector{UInt8}, snp2::Vector{UInt8})::Float64
    # drop missing values
    idx = (snp1 .!= @MISSING) .& (snp2 .!= @MISSING)
    s1, s2 = snp1[idx], snp2[idx]
    
    # number of not missing samples
    n = length(idx)
    
    # calculate needed frequencies
    naa = naA = nAA = nbb = nbB = nBB =
		nAABB = naabb = naaBB = nAAbb = 0

    for (g1, g2) in zip(s1, s2)
        if  g1 == @HOMO_M
            nAA +=1
            if g2 == @HOMO_M # AABB
                nBB +=1
                nAABB +=1
            elseif g2 == @HOMO_m # AAbb
                nbb +=1
                nAAbb +=1
            else    #HETER -> AAbB
                nbB +=1
            end
        elseif g1 == @HOMO_m
            naa +=1
            if g2 == @HOMO_M # aaBB
                nBB +=1
                naaBB +=1
            elseif g2 == @HOMO_m  # aabb
                nbb +=1
                naabb +=1
            else    # HETER -> aabB
                nbB +=1
            end
        else # @HETER
            naA +=1
            if g2 == @HOMO_M # aABB
                nBB +=1
            elseif g2 == @HOMO_m # aAbb
                nbb+=1
            else    # HETER -> aAbB
                nbB +=1
            end
        end
    end

    # final calculations
    Δ = (nAABB + naabb - naaBB - nAAbb) / (2*n) - 
        (naa-nAA)*(nbb-nBB) / (2*n*n)
    
    pa = (2*naa + naA) / (2*n)
    pb = (2*nbb + nbB) / (2*n)
    pA = 1 - pa
    pB = 1 - pb
    pAA = nAA / n
    pBB = nBB / n
    DA = pAA - pA*pA
    DB = pBB - pB*pB
    t = (pA*pa + DA) * (pB*pb + DB)

    return Δ^2 / t

end


"""
LD r² composite of pair of SNPs given genotype matrix and indices of snps in matrix
"""
function ld_r²(snp1::Integer, snp2::Integer, ref::AbstractSnpArray)::Float64
    return ld_r²(ref[:, snp1], ref[:, snp2])    
end


"""
LD r² composite of pair of SNPs given snp ids and refernce data
"""
function ld_r²(snp1::AbstractString, snp2::AbstractString, ref::SnpData; formated = false)::Float64
    if !formated
        formatSnpData!(ref, :snpid)
    end

    id1 = searchsortedfirst(ref.snp_info.snpid, snp1)
    id2 = searchsortedfirst(ref.snp_info.snpid, snp2)

    if (id1 > lastindex(ref.snp_info.snpid) || 
        id2 > lastindex(ref.snp_info.snpid) ||
        ref.snp_info.snpid[id1] != snp1 ||
        ref.snp_info.snpid[id2] != snp2)

        return NaN
    else
        return ld_r²(ref.snp_info.idx[id1], ref.snp_info.idx[id2], ref.snparray)
    end
end


"""
LD r² composite of pair of SNPs given (chromosome, position) and reference data
"""
function ld_r²(snp1::Tuple{Integer, Integer}, snp2::Tuple{Integer, Integer}, ref::SnpData, formated = false)
    if !formated
        formatSnpData!(ref)
    end

    id1 = searchsortedfirst(ref.snp_info.chr_pos, snp1)
    id2 = searchsortedfirst(ref.snp_info.chr_pos, snp2)

    if (id1 > lastindex(ref.snp_info.chr_pos) || 
        id2 > lastindex(ref.snp_info.chr_pos) ||
        ref.snp_info.chr_pos[id1] != snp1 ||
        ref.snp_info.chr_pos[id2] != snp2)
        
        return NaN
    else
        return ld_r²(ref.snp_info.idx[id1], ref.snp_info.idx[id2], ref.snparray)
    end
end


"""
LD r² composite matrix of n SNPs from index in SnpArray
"""
function getLDmat(arr::SnpArray, idx::AbstractVector{Int})::Matrix{Float64}
    M_corr = Matrix{Float64}(I, length(idx), length(idx))
    
    d = Dict(zip(idx, 1:lastindex(idx)))
    @threads for (i, j) in collect(subsets(idx, 2))
        snp1 = arr[:,i]
        snp2 = arr[:,j]
        M_corr[d[i], d[j]] = M_corr[d[j], d[i]] = ld_r²(snp1, snp2) #function implemented following paper pmid :18757931, for r² type r\^2
    end
    return M_corr
end


"""
Get correlation Matrix for specifies snps (tuple of the form (chr, pos) 
    given reference genotype SnpData.
"""
function getLDmat(ref_genotypes::SnpData, 
                  snps::AbstractVector{<:Tuple{Integer, Integer}},
                  formated::Bool = false
                  )::Tuple{Matrix{Float64}, Vector{Bool}}

    if !formated
        formatSnpData!(ref_genotypes)
    end

    snps_indx = Vector{Union{Int}}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        local j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        if j > length(ref_genotypes.snp_info.chr_pos) || ref_genotypes.snp_info.chr_pos[j] != chr_pos_sing
            j = -1
        end
        snps_indx[i] = (j < 0) ? j : ref_genotypes.snp_info.idx[j]
    end
    kept_indx = snps_indx[snps_indx .> 0]

    return getLDmat(ref_genotypes.snparray, kept_indx), snps_indx .> 0
end

"""
Get correlation Matrix for specifies snps (snpid) 
    given reference genotype SnpData.
"""
function getLDmat(ref_genotypes::SnpData,
    snps::AbstractVector{<:AbstractString},
    formated::Bool = false
    )::Tuple{Matrix{Float64}, Vector{Bool}}
    
    if !formated
        formatSnpData!(ref_genotypes, :snpid)
    end

    snps_indx = Vector{Union{Int}}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        local j = searchsortedfirst(ref_genotypes.snp_info.snpid, chr_pos_sing)
        if j > lastindex(ref_genotypes.snp_info.snpid) || ref_genotypes.snp_info.snpid[j] != chr_pos_sing
            j = -1
        end
        snps_indx[i] = (j < 0) ? j : ref_genotypes.snp_info.idx[j]
    end
    kept_indx = snps_indx[snps_indx .> 0]

    return getLDmat(ref_genotypes.snparray, kept_indx), snps_indx .> 0
end


"""
format Genotype information contained in SnpData for snp sorted search based on chromosome and position (snpid).
Adds a column chr_pos (chr::Int8, pos::Int) in snp_info and sorts snp_info by chr_pos (snpid). Adds a column containing original indices.
    returns nothing
"""
function formatSnpData!(Genotypes::SnpData, sort_by::Symbol = :chr_pos)
    if !hasproperty(Genotypes.snp_info, :idx)
        Genotypes.snp_info.idx = collect(1:size(Genotypes.snp_info, 1))
    end
    if !hasproperty(Genotypes.snp_info, :chr_pos)
        Genotypes.snp_info.chr_pos = collect(
                zip(parse.(Int8, Genotypes.snp_info.chromosome), 
                    Genotypes.snp_info.position)
            )
    end
    if !issorted(getproperty(Genotypes.snp_info, sort_by))
        sort!(Genotypes.snp_info, sort_by)
    end
end

