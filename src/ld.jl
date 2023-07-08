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
@inline function _ld_r²_(s1, s2)::Float64
    n = 0
    
    # calculate needed frequencies
    naa = naA = nAA = nbb = nbB = nBB =
		nAABB = naabb = naaBB = nAAbb = 0

    @inbounds @simd for i in eachindex(s1, s2)
        g1, g2 = s1[i], s2[i]
        if  g1 == @HOMO_M
            if g2 == @HOMO_M # AABB
                nBB +=1
                nAABB +=1
                n += 1
                nAA +=1
            elseif g2 == @HOMO_m # AAbb
                nbb +=1
                nAAbb +=1
                n += 1
                nAA +=1
            elseif g2 == @HETER    #HETER -> AAbB
                nbB +=1
                n += 1
                nAA +=1
            end
        elseif g1 == @HOMO_m
            if g2 == @HOMO_M # aaBB
                nBB +=1
                naaBB +=1
                n += 1
                naa +=1
            elseif g2 == @HOMO_m  # aabb
                nbb +=1
                naabb +=1
                n += 1
                naa +=1
            elseif g2 == @HETER    # HETER -> aabB
                nbB +=1
                n += 1
                naa +=1
            end
        elseif g1 == @HETER # @HETER
            if g2 == @HOMO_M # aABB
                nBB +=1
                n += 1
                naA +=1
            elseif g2 == @HOMO_m # aAbb
                nbb+=1
                n += 1
                naA +=1
            elseif g2 == @HETER    # HETER -> aAbB
                nbB +=1
                n += 1
                naA +=1
            end
        end
    end

    # final calculations
    @fastmath begin
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

        r = Δ^2 / t
    end

    return r
end


function ld_r²(s1, s2)::Float64
    n = 0
    
    # calculate needed frequencies
    naa = naA = nAA = nbb = nbB = nBB =
		nAABB = naabb = naaBB = nAAbb = 0

    for i in eachindex(s1, s2)
        g1, g2 = s1[i], s2[i]
        if  g1 == @HOMO_M
            if g2 == @HOMO_M # AABB
                nBB +=1
                nAABB +=1
                n += 1
                nAA +=1
            elseif g2 == @HOMO_m # AAbb
                nbb +=1
                nAAbb +=1
                n += 1
                nAA +=1
            elseif g2 == @HETER    #HETER -> AAbB
                nbB +=1
                n += 1
                nAA +=1
            end
        elseif g1 == @HOMO_m
            if g2 == @HOMO_M # aaBB
                nBB +=1
                naaBB +=1
                n += 1
                naa +=1
            elseif g2 == @HOMO_m  # aabb
                nbb +=1
                naabb +=1
                n += 1
                naa +=1
            elseif g2 == @HETER    # HETER -> aabB
                nbB +=1
                n += 1
                naa +=1
            end
        elseif g1 == @HETER # @HETER
            if g2 == @HOMO_M # aABB
                nBB +=1
                n += 1
                naA +=1
            elseif g2 == @HOMO_m # aAbb
                nbb+=1
                n += 1
                naA +=1
            elseif g2 == @HETER    # HETER -> aAbB
                nbB +=1
                n += 1
                naA +=1
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
    return @views _ld_r²_(ref[:, snp1], ref[:, snp2])
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
        return @views _ld_r²_(ref.snparray[:, ref.snp_info.idx[id1]], ref.snparray[:, ref.snp_info.idx[id2]])
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
        return @views _ld_r²_(ref.snparray[:, ref.snp_info.idx[id1]], ref.snparray[:, ref.snp_info.idx[id2]])
    end
end


"""
LD r² composite matrix of n SNPs from index in given SnpArray
"""
function getLDmat(arr::SnpArray, idx::AbstractVector{<:Integer})::Matrix{Float64}
    M_corr = Matrix{Float64}(I, length(idx), length(idx))
    
    @threads for (i, j) in collect(subsets(eachindex(idx), 2))
        @inbounds M_corr[i, j] = M_corr[j, i] = @views _ld_r²_(arr[:, idx[i]], arr[:, idx[j]]) #function implemented following paper pmid :18757931, for r² type r\^2
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

    # snps_indx = Vector{Union{Int}}(undef, size(snps, 1))
    # @threads for (i, chr_pos_sing) in collect(enumerate(snps))
    #     local j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
    #     if j > length(ref_genotypes.snp_info.chr_pos) || ref_genotypes.snp_info.chr_pos[j] != chr_pos_sing
    #         j = -1
    #     end
    #     snps_indx[i] = (j < 0) ? j : ref_genotypes.snp_info.idx[j]
    # end
    # kept_indx = snps_indx[snps_indx .> 0]
    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    kept_indx = Vector{Bool}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.chr_pos, chr_pos_sing)
        @inbounds begin 
            kept_indx[i] = firstindex(ref_genotypes.snp_info.chr_pos) ≤ j ≤ lastindex(ref_genotypes.snp_info.chr_pos) && 
                          ref_genotypes.snp_info.chr_pos[j] == chr_pos_sing
            if kept_indx[i]
                snps_indx[i] = ref_genotypes.snp_info.idx[j]
            end
        end
    end

    return getLDmat(ref_genotypes.snparray, snps_indx[kept_indx]), kept_indx
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

    # snps_indx = Vector{Union{Int}}(undef, size(snps, 1))
    # @threads for (i, chr_pos_sing) in collect(enumerate(snps))
    #     local j = searchsortedfirst(ref_genotypes.snp_info.snpid, chr_pos_sing)
    #     if j > lastindex(ref_genotypes.snp_info.snpid) || ref_genotypes.snp_info.snpid[j] != chr_pos_sing
    #         j = -1
    #     end
    #     snps_indx[i] = (j < 0) ? j : ref_genotypes.snp_info.idx[j]
    # end
    # kept_indx = snps_indx[snps_indx .> 0]
    # search for indices of given snps
    snps_indx = Vector{Int}(undef, size(snps, 1))
    kept_indx = Vector{Bool}(undef, size(snps, 1))
    @threads for (i, chr_pos_sing) in collect(enumerate(snps))
        j = searchsortedfirst(ref_genotypes.snp_info.snpid, chr_pos_sing)
        @inbounds begin 
            kept_indx[i] = firstindex(ref_genotypes.snp_info.snpid) ≤ j ≤ lastindex(ref_genotypes.snp_info.snpid) && 
                          ref_genotypes.snp_info.snpid[j] == chr_pos_sing
            if kept_indx[i]
                snps_indx[i] = ref_genotypes.snp_info.idx[j]
            end
        end
    end

    return getLDmat(ref_genotypes.snparray, snps_indx[kept_indx]), kept_indx
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

