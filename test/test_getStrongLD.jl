data_dir_plink_files = joinpath(pwd(), "data", "sample_chosen")

chosen_sample = [(1, 87917746),
                 (1, 100046246),
                 (1, 170645774),
                 (1, 201746768)]
chosen_sample_ids = ["1:87917746",
                     "1:100046246",
                     "1:170645774",
                     "1:201746768"]

@testset "getStrongLD tuple" begin
    snp = (1, 100046246)
    data = data = SnpData(SnpArrays.datadir(data_dir_plink_files))
    r1 = getStrongLD(data, [snp])
    r2 = getStrongLD(repeat([data], 22), [snp])

    @test r1 == r2 == [snp]
    
    r1 = getStrongLD(data, [snp], window = 50000000, r2_tresh = 0.0)
    r2 = getStrongLD(repeat([data], 22), [snp], window = 50000000, r2_tresh = 0.0)

    @test Set(r1) == Set(r2) == Set([(1, 87917746), snp])

    r1 = getStrongLD(data, [snp], window = Int(1e12), r2_tresh = 0.002)
    r2 = getStrongLD(repeat([data], 22), [snp], window = Int(1e12), r2_tresh = 0.002)

    @test Set(r1) == Set(r2) == Set([(1, 87917746),
                                     (1, 100046246),
                                     (1, 201746768)])

    r1 = getStrongLD(data, [snp], window = Int(1e12), r2_tresh = 0.0)
    r2 = getStrongLD(repeat([data], 22), [snp], window = Int(1e12), r2_tresh = 0.0)

    @test Set(r1) == Set(r2) == Set(chosen_sample)

    r = getStrongLD(repeat([data], 22), [snp, (1, 1234), (2, 2345)])
    
    @test r == [snp]
end


@testset "getStrongLD id" begin
    snp = "1:100046246"
    data = data = SnpData(SnpArrays.datadir(data_dir_plink_files))
    r1 = getStrongLD(data, [snp])
    r2 = getStrongLD(repeat([data], 22), [snp])

    @test r1 == r2 == [snp]
    
    r1 = getStrongLD(data, [snp], window = 50000000, r2_tresh = 0.0)
    r2 = getStrongLD(repeat([data], 22), [snp], window = 50000000, r2_tresh = 0.0)

    @test Set(r1) == Set(r2) == Set(["1:87917746", snp])

    r1 = getStrongLD(data, [snp], window = Int(1e12), r2_tresh = 0.002)
    r2 = getStrongLD(repeat([data], 22), [snp], window = Int(1e12), r2_tresh = 0.002)

    @test Set(r1) == Set(r2) == Set(["1:87917746",
                                     "1:100046246",
                                     "1:201746768"])

    r1 = getStrongLD(data, [snp], window = Int(1e12), r2_tresh = 0.0)
    r2 = getStrongLD(repeat([data], 22), [snp], window = Int(1e12), r2_tresh = 0.0)

    @test Set(r1) == Set(r2) == Set(chosen_sample_ids)

    r = getStrongLD(repeat([data], 22), [snp, "1:1234", "2:2345"])
    
    @test r == [snp]
end