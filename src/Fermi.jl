# Interface with Fermi.jl for PPP -> Fermi wavefunction handoff

using PPP
using Fermi

using LinearAlgebra
using Printf

"""
    to_fermi_rhf(system::PPP.MolecularSystem, scf_result::PPP.SCFResult)

Construct a `Fermi.HartreeFock.RHF` wavefunction from PPP SCF results with
consistent dimensions and electron count.

Notes:
- Uses an identity AO overlap (orthonormal π-site basis).
- Sets the molecule charge so that total electrons match 2×ndocc.
- Uses multiplicity 1 for closed-shell (even-electron) cases, 2 otherwise.
"""
function to_fermi_rhf(system::PPP.MolecularSystem, scf_result::PPP.SCFResult)
    # Dimensions
    n_sites = size(scf_result.eigenvectors, 1)
    ndocc = scf_result.n_occupied
    nvir = n_sites - ndocc

    # OK, this is utterly horrific, but I got it working. 
#   Essentially we put a Hydrogen everywhere there is a PPP site, and then request a STO-3G basis, which tricks all of the Fermi machinery which assumes Gaussian orbitals into treatin the PPP orbitals correctly. I am so so sorry ~ Jarv.
# We use element H for all sites; the charge is adjusted to match electrons = 2*ndocc.
    coords = [a.position for a in system.atoms]
    @assert length(coords) == n_sites
    # Compose XYZ-like molstring for Fermi.Molecule
    buf = IOBuffer()
    for r in coords
        # element x y z
        println(buf, "H $(r[1]) $(r[2]) $(r[3])")
    end
    molstring = String(take!(buf))

    total_Z = n_sites                    # Z=1 per H
    total_electrons = 2 * ndocc          # closed-shell electrons for PPP
    charge = total_Z - total_electrons   # electrons = Z - charge
    multiplicity = iseven(total_electrons) ? 1 : 2
    # Ensure Fermi global options reflect this molecule, so downstream integral helpers
    # use consistent AO basis and electron count
    Fermi.Options.set("molstring", molstring)
    Fermi.Options.set("unit", "angstrom")
    Fermi.Options.set("charge", charge)
    Fermi.Options.set("multiplicity", multiplicity)
    Fermi.Options.set("basis", "sto-3g")

    molecule = Fermi.Molecule(
        molstring=molstring,
        unit=:angstrom,
        charge=charge,
        multiplicity=multiplicity,
    )

    # MO data from PPP SCF (convert eV → Hartree for Fermi)
    ev2au = 1 / 27.21138505
    mo_energies = scf_result.energies .* ev2au
    mo_coefficients = scf_result.eigenvectors

    # Fermi.Orbitals.RHFOrbitals expects:
    #   (molecule::Molecule, basis::String, eps::Vector{Float64}, sd_energy::Float64, C::Matrix{Float64})
    orbitals = Fermi.Orbitals.RHFOrbitals(
        molecule,
        "sto-3g",
        mo_energies,
        scf_result.total_energy * ev2au,
        mo_coefficients,
    )

    hf_energy = scf_result.total_energy * ev2au
    return Fermi.HartreeFock.RHF(molecule, hf_energy, ndocc, nvir, orbitals, 0.0, 0.0)
end

# Example usage (uncomment to run manually):
using PPP
sys, huckel, scf = PPP.run_ppp_calculation("molecules/ethene-2C.xyz", PPP.Bedogni2024ModelParams())
wfn = to_fermi_rhf(sys, scf)
mp2 = @energy wfn => mp2

# CCSD: create integral helpers explicitly (does not dispatch on RHF)
moints = Fermi.Integrals.IntegralHelper{Float64}(orbitals=wfn.orbitals)
aoints = Fermi.Integrals.IntegralHelper{Float64}(molecule=wfn.molecule)
ccsd = Fermi.CoupledCluster.RCCSD(moints, aoints, Fermi.CoupledCluster.get_rccsd_alg())
display(ccsd)

# FCI: API expects (aoints, rhf)
fci = Fermi.ConfigurationInteraction.RFCI(aoints, wfn, Fermi.ConfigurationInteraction.get_rfci_alg())
display(fci)

Ha2eV = 27.21138505
# Compare total and correlation energies across methods
println("\nTotal Energies (eV):")
println("Huckel: ", huckel.total_energy)
println("PPP: ", scf.total_energy)

println("SCF:  ", wfn.energy * Ha2eV)
println("MP2:  ", mp2.energy * Ha2eV)
println("CCSD: ", ccsd.energy * Ha2eV)
println("FCI:  ", fci.energy * Ha2eV)

println("\nCorrelation Energies (eV):")
println("MP2:  ", mp2.energy * Ha2eV - wfn.energy * Ha2eV)
println("CCSD: ", ccsd.energy * Ha2eV - wfn.energy * Ha2eV) 
println("FCI:  ", fci.energy * Ha2eV - wfn.energy * Ha2eV)


#Reference values courtesy of PlotDigitizer, lightly edited. Fig 4d: 2T-7N at different levels of theory

#singlets=cis_tda_singlet(wfn; nroots=5)
#triplets=cis_tda_triplet(wfn; nroots=5)
#println("TDA-CIS singlet roots (eV): ", singlets)
#println("TDA-CIS triplet roots (eV): ", triplets)
   
# Copy paste so we can see how close we are getting. 
#@printf("T1: %.4f T1(Bedogni2024): %.4f diff: %.4f\n", triplets[1], 3.428550679531072, triplets[1] - 3.428550679531072)
#@printf("S1: %.4f S1(Bedogni2024): %.4f diff: %.4f\n", singlets[1], 4.052910052910053, singlets[1] - 4.052910052910053)
#@printf("T3: %.4f T3(Bedogni2024): %.4f diff: %.4f\n", triplets[3], 3.947089947089947, triplets[3] - 3.947089947089947)
#@printf("S2: %.4f S2(Bedogni2024): %.4f diff: %.4f\n", singlets[2], 5.216931216931217, singlets[2] - 5.216931216931217)

