
data_dir_plink_files = joinpath(pwd(), "data", "sample_chosen")

chosen_sample = [(1, 87917746),
                 (1, 100046246),
                 (1, 170645774),
                 (1, 201746768)]
chosen_sample_ids = ["1:87917746",
                     "1:100046246",
                     "1:170645774",
                     "1:201746768"]

@testset "formatSnpData!" begin
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))
    formatSnpData!(data)

    @test issorted(data.snp_info.chr_pos)
    @test hasproperty(data.snp_info, :idx)
    
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))
    formatSnpData!(data, :snpid)

    @test issorted(data.snp_info.snpid)
    @test hasproperty(data.snp_info, :idx)
end

@testset "ld_r²" begin
    resp = 0.00124
    tol = 1e-5
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))

    @test isapprox(ld_r²(data.snparray[:, 1], data.snparray[:, 4]), resp, atol = tol)
    @test isapprox(ld_r²(1, 4, data.snparray), resp, atol = tol)

    @test isapprox(ld_r²(chosen_sample[1], chosen_sample[4], data), resp, atol = tol)
    @test isapprox(ld_r²(chosen_sample_ids[1], chosen_sample_ids[4], data), resp, atol = tol)

end