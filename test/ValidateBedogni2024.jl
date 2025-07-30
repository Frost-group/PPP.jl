using PPP
using Test

@testset "Validate Bedogni2024ModelParams with 2T-7N" begin
    # Test that default parameters work without explicit parameter passing
    system, huckel_result, scf_result = run_ppp_calculation("../molecules/2T-7N.xyz", Bedogni2024ModelParams())
    
    @test scf_result.converged
    @test length(system.atoms) == 13
    @test system.n_electrons == 14        
    # Verify atom composition
    carbon_count = count(atom -> atom.symbol == :C, system.atoms)
    nitrogen_count = count(atom -> atom.symbol == :N, system.atoms)
    @test carbon_count == 6
    @test nitrogen_count == 7

# From raw output on GitHub of Bedogni2024 
# Huckel MO energies (eV)
# 7    -5.00000
# 8     0.54229
    @test huckel_result.energies[7] ≈ -5.000 atol=1e-3
    @test huckel_result.energies[8] ≈ 0.542 atol=1e-3
    @test huckel_result.energies[8] - huckel_result.energies[7] ≈ 5.542 atol=1e-3
    @test huckel_result.total_energy ≈ -107.101913355153 atol=1e-3

# SCF convergence reached after     62  iterations, with threshold=   0.10000E-11
# GS energy=  -47.0167692550647      eV
# HOMO energy=  2.751288921716742E-003 eV
# LUMO energy=   10.0450931042920      eV
# HOMO-LUMO gap=   10.0423418153703      eV
# HOMO-LUMO exchange integral= -0.216408175656806      eV

    homo_lumo_gap = scf_result.energies[8] - scf_result.energies[7]
    @test homo_lumo_gap ≈ 10.0423418153703 atol=1e-2

    @test PPP.calculate_HOMO_LUMO_exchange(scf_result, system) ≈ -0.216408175656806 atol=1e-2

    @test scf_result.total_energy ≈ -47.0167692550647 atol=1e-2

end
