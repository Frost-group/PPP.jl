using Test
using PPP

@testset "Test futzing Bedogni2024ModelParams ... " begin
    @testset "Read and write params..." begin
        params = Bedogni2024ModelParams()
        @test params.T == -2.4
        @test params.UC == 11.26
        @test params.cutoff == 1.4
        
        # Test custom parameters
        custom_params = Bedogni2024ModelParams(
            T = -3.0,
            UC = 12.0,
            cutoff = 1.5
        )
        @test custom_params.T == -3.0
        @test custom_params.UC == 12.0
        @test custom_params.cutoff == 1.5
    end
    
    @testset "Transfer integral sensitivity..." begin
        # Test with different hopping integrals
        params1 = Bedogni2024ModelParams(T = -2.0)
        params2 = Bedogni2024ModelParams(T = -3.0)
        
        _, _, scf1 = run_ppp_calculation("../molecules/2T-7N.xyz", params1)
        _, _, scf2 = run_ppp_calculation("../molecules/2T-7N.xyz", params2)
        
        gap1 = scf1.energies[8] - scf1.energies[7]
        gap2 = scf2.energies[8] - scf2.energies[7]
        
        # Different parameters should give different results
        @test abs(gap1 - gap2) > 1e-2
    end
    
    @testset "Futz Bedogni2024ModelParams and see if it converges..." begin
        # Test with modified Hubbard U for carbon
        params_uc = Bedogni2024ModelParams(UC = 15.0)
        _, _, scf_uc = run_ppp_calculation("../molecules/2T-7N.xyz", params_uc)
        
        # Test with modified Hubbard U for nitrogen
        params_un = Bedogni2024ModelParams(UN = 18.0, UN_AZA = 18.5)
        _, _, scf_un = run_ppp_calculation("../molecules/2T-7N.xyz", params_un)
        
        # Test with modified site energy for nitrogen
        params_es = Bedogni2024ModelParams(ES_NPY = -15.0, ES_NAZA = -7.0)
        _, _, scf_es = run_ppp_calculation("../molecules/2T-7N.xyz", params_es)
        
        # Test with modified bond cutoff: though if you start actually getting extra atoms, the whole thing will break
        params_cutoff = Bedogni2024ModelParams(cutoff = 1.6)
        _, _, scf_cutoff = run_ppp_calculation("../molecules/2T-7N.xyz", params_cutoff)
        
        # All should converge
        @test scf_uc.converged
        @test scf_un.converged
        @test scf_es.converged
        @test scf_cutoff.converged
    end
    

 
end 