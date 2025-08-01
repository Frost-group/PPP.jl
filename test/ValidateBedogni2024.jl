using PPP
using Test
using Printf

@testset "Validate Bedogni2024ModelParams with 2T-7N c.f. Bedogni2024 published values" begin
    system, huckel_result, scf_result = run_ppp_calculation("../molecules/2T-7N.xyz", Bedogni2024ModelParams())
    
    @test scf_result.converged
    @test length(system.atoms) == 13
    @test system.n_electrons == 14        
    # Verify atom composition
    carbon_count = count(atom -> atom.symbol == :C, system.atoms)
    nitrogen_count = count(atom -> atom.symbol == :N, system.atoms)
    @test carbon_count == 6
    @test nitrogen_count == 7

# From raw output on GitHub of Bedogni2024 : https://github.com/francescodimaiolo/Hartree-Fock_PPP_tool/blob/main/output.out
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


# Validate CIS calculation
@testset "Validate CIS calculation with 2T-7N c.f. Bedogni2024 published values" begin
    system, huckel_result, scf_result = run_ppp_calculation("../molecules/2T-7N.xyz", Bedogni2024ModelParams())
    singlet, triplet = run_cis_calculation(system, scf_result)

#Reference values courtesy of PlotDigitizer, lightly edited. Fig 4d: 2T-7N at different levels of theory
    @test triplet.energies[1] ≈ 3.428550679531072 atol=1e-2
    @test triplet.energies[2] ≈ 3.947089947089947 atol=1e-2

    @test singlet.energies[1] ≈ 4.052910052910053 atol=1e-2
    @test singlet.energies[2] ≈ 5.216931216931217 atol=1e-2

    # Copy paste so we can see how close we are getting. 
    @printf("T1: %.4f T1(Bedogni2024): %.4f diff: %.4f\n", triplet.energies[1], 3.428550679531072, triplet.energies[1] - 3.428550679531072)
    @printf("S1: %.4f S1(Bedogni2024): %.4f diff: %.4f\n", singlet.energies[1], 4.052910052910053, singlet.energies[1] - 4.052910052910053)
    @printf("T2: %.4f T2(Bedogni2024): %.4f diff: %.4f\n", triplet.energies[2], 3.947089947089947, triplet.energies[2] - 3.947089947089947)
    @printf("S2: %.4f S2(Bedogni2024): %.4f diff: %.4f\n", singlet.energies[2], 5.216931216931217, singlet.energies[2] - 5.216931216931217)
end