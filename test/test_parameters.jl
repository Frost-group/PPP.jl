using Test
using PPP
using PPP: BEDOGNI_HUBBARD_U, BEDOGNI_SITE_ENERGY

@testset "Test futzing Bedogni2024Model ... " begin
    @testset "Read and write params..." begin
        params = Bedogni2024Model()
        @test params.T == -2.4
        @test params.HUBBARD_U[:C] == 11.26
        # Note: cutoff is now hardcoded in calculate_connectivity, not a model parameter
        
        # Test custom parameters with modified dictionaries
        custom_hubbard = copy(BEDOGNI_HUBBARD_U)
        custom_hubbard[:C] = 12.0
        custom_params = Bedogni2024Model(
            T = -3.0,
            HUBBARD_U = custom_hubbard
        )
        @test custom_params.T == -3.0
        @test custom_params.HUBBARD_U[:C] == 12.0
    end
    
    @testset "Transfer integral sensitivity..." begin
        # Test with different hopping integrals
        params1 = Bedogni2024Model(T = -2.0)
        params2 = Bedogni2024Model(T = -3.0)
        
        _, _, scf1 = run_ppp_calculation("../molecules/2T-7N.xyz", params1)
        _, _, scf2 = run_ppp_calculation("../molecules/2T-7N.xyz", params2)
        
        gap1 = scf1.energies[8] - scf1.energies[7]
        gap2 = scf2.energies[8] - scf2.energies[7]
        
        # Different parameters should give different results
        @test abs(gap1 - gap2) > 1e-2
    end
    
    @testset "Futz Bedogni2024Model and see if it converges..." begin
        # Test with modified Hubbard U for carbon
        hubbard_uc = copy(BEDOGNI_HUBBARD_U)
        hubbard_uc[:C] = 15.0
        params_uc = Bedogni2024Model(HUBBARD_U = hubbard_uc)
        _, _, scf_uc = run_ppp_calculation("../molecules/2T-7N.xyz", params_uc)
        
        # Test with modified Hubbard U for nitrogen
        hubbard_un = copy(BEDOGNI_HUBBARD_U)
        hubbard_un[:N2] = 18.0  # Pyrrole
        hubbard_un[:N1] = 18.5  # Aza
        params_un = Bedogni2024Model(HUBBARD_U = hubbard_un)
        _, _, scf_un = run_ppp_calculation("../molecules/2T-7N.xyz", params_un)
        
        # Test with modified site energy for nitrogen
        site_es = copy(BEDOGNI_SITE_ENERGY)
        site_es[:N2] = -15.0  # Pyrrole
        site_es[:N1] = -7.0   # Aza
        params_es = Bedogni2024Model(SITE_ENERGY = site_es)
        _, _, scf_es = run_ppp_calculation("../molecules/2T-7N.xyz", params_es)
        
        # Note: cutoff is now hardcoded in calculate_connectivity, not a model parameter
        
        # All should converge
        @test scf_uc.converged
        @test scf_un.converged
        @test scf_es.converged
    end
    

 
end 