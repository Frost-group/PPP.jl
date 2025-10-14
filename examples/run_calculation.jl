using PPP
using Printf

# Run basic PPP calculation
println("\n" * "="^50)
println("Running PPP SCF Calculation")
println("="^50)
system, huckel_result, scf_result = PPP.run_ppp_calculation("molecules/2T-N.xyz", Bedogni2024Model())

# Run CIS calculation for singlet
println("\n" * "="^50)
println("Running CIS Calculation (singlet)")
println("="^50)
singlet, triplet = PPP.run_cis_calculation(system, scf_result)
n_singlets = length(singlet.energies)
n_triplets = length(triplet.energies)

# Print CIS analysis
println("\nCIS Singlet Analysis:")
println("----")
for i in 1:min(n_singlets,5)  # Shows up to 5 singlet states
    config = singlet.configurations[singlet.dominant_configurations[i]]
    @printf("State %d: ΔE = %.3f eV, f = %.4f, %d→%d excitation\n",
            i, singlet.energies[i], singlet.oscillator_strengths[i],
            config.from_orbitals[1], config.to_orbitals[1])
end

println("\nCIS Triplet Analysis:")
println("----")
for i in 1:min(n_triplets,5) # Shows up to 5 triplet states
    config = triplet.configurations[triplet.dominant_configurations[i]]
    @printf("State %d: ΔE = %.3f eV, f = %.4f, %d→%d excitation\n",
            i, triplet.energies[i], triplet.oscillator_strengths[i],
            config.from_orbitals[1], config.to_orbitals[1])
end

println("\nCIS Gap Analysis:")
println("----")
for i in 1:min(n_singlets,5)  # Shows up to 5 gaps
    @printf("ΔE = %.3f eV\n",
            singlet.energies[i]-triplet.energies[i])
end

dsp_correction = calculate_dsp(system, scf_result, singlet, triplet)

println("\n" * "="^50)
println("DSP Correction = ", dsp_correction)
println("="^50)

println("\nDSP corrected CIS Gap Analysis:")
println("----")
for i in 1:min(n_singlets,5)  # Shows up to 5 gaps
    @printf("ΔE = %.3f eV\n",
            singlet.energies[i]-triplet.energies[i]+dsp_correction)
end


# Run CISD calculation
println("\n" * "="^50)
println("Running CISD Calculation")
println("="^50)
cisd_result = PPP.run_cisd_calculation(system, scf_result)
n_cisd = length(cisd_result.energies)

# Print CISD analysis
println("\nCISD Analysis:")
println("-----")
for i in 1:min(n_cisd,5)  # Show up to 5 states
    config = cisd_result.configurations[cisd_result.dominant_configurations[i]]
    if config.type isa SingleExcitation
        @printf("State %d: ΔE = %.3f eV, f = %.4f, %d→%d excitation\n",
                i, cisd_result.energies[i], cisd_result.oscillator_strengths[i],
                config.from_orbitals[1], config.to_orbitals[1])
    else
        @printf("State %d: ΔE = %.3f eV, f = %.4f, (%d,%d)→(%d,%d) double excitation\n",
                i, cisd_result.energies[i], cisd_result.oscillator_strengths[i],
                config.from_orbitals[1], config.from_orbitals[2],
                config.to_orbitals[1], config.to_orbitals[2])
    end
end

