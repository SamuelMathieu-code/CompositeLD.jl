
data_dir_plink_files = joinpath(pwd(), "data", "sample_chosen")

chosen_sample = [(1, 87917746),
                 (1, 100046246),
                 (1, 170645774),
                 (1, 201746768),
                 (0, 0)]
chosen_sample_ids = ["1:87917746",
                     "1:100046246",
                     "1:170645774",
                     "1:201746768",
                     "0:0"]

@testset "formatSnpData!" begin
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))
    formatSnpData!(data)

    @test issorted(data.snp_info.chr_pos)
    @test hasproperty(data.snp_info, :idx)
    @test hasproperty(data.snp_info, :chr_pos)
    
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))
    formatSnpData!(data, :snpid)

    @test issorted(data.snp_info.snpid)
    @test hasproperty(data.snp_info, :idx)
    @test hasproperty(data.snp_info, :chr_pos)
end

@testset "ld_r2" begin
    resp = 0.00124
    tol = 1e-5
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))

    @test isapprox(ld_r2(data.snparray[:, 1], data.snparray[:, 4]), resp, atol = tol)
    @test isapprox(ld_r2(1, 4, data.snparray), resp, atol = tol)

    @test isapprox(ld_r2(chosen_sample[1], chosen_sample[4], data), resp, atol = tol)
    @test isapprox(ld_r2(chosen_sample_ids[1], chosen_sample_ids[4], data), resp, atol = tol)

end

@testset "getLDmat" begin
    resp = [1.0         0.00213594  0.00109997   0.00124035;
            0.00213594  1.0         0.00090874   0.00271792;
            0.00109997  0.00090874  1.0          0.000142889;
            0.00124035  0.00271792  0.000142889  1.0]
    tol = 1e-5
    data = SnpData(SnpArrays.datadir(data_dir_plink_files))

    @test isapprox(resp, getLDmat(data.snparray, [1, 2, 3, 4]))
    formatSnpData!(data, :chr_pos)
    m, _ = getLDmat(data, chosen_sample, true)
    @test isapprox(resp, m)
    formatSnpData!(data, :snpid)
    m, _ = getLDmat(data, chosen_sample_ids, true)
    @test isapprox(resp, m)

    m, _ = getLDmat(data, chosen_sample)
    @test isapprox(resp, m)
    m, _ = getLDmat(data, chosen_sample)
    @test isapprox(resp, m) # test for fromatSnpData! run twice :chr_pos
    m, _ = getLDmat(data, chosen_sample_ids)
    @test isapprox(resp, m)
    m, _ = getLDmat(data, chosen_sample_ids)
    @test isapprox(resp, m) # test for formatSnpData! run twice :snpid
    
end