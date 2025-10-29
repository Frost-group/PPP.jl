# Interface with Fermi.jl for PPP -> Fermi wavefunction handoff

using PPP
using Fermi

using LinearAlgebra
using Printf

 

# ============================================================================ #
# Pack PPP integrals into Fermi IntegralHelpers (PPP → Fermi moints/aoints)
# ============================================================================ #

"""
    to_fermi_ppp_integrals(system::PPP.MolecularSystem,
                           huckel::PPP.HuckelResult,
                           scf::PPP.SCFResult)
        -> (rhf::Fermi.HartreeFock.RHF,
            moints::Fermi.Integrals.IntegralHelper,
            aoints::Fermi.Integrals.IntegralHelper)

Create a Fermi RHF wavefunction and integral helpers backed by PPP data (ZDO/Ohno):
- Orbital coefficients from PPP `scf.eigenvectors`
- Orbital energies and reference energy from PPP in Hartree
- AO ERIs in Chonky 4-index form with ZDO: (μν|ρσ) = δ_{μν} δ_{ρσ} V_raw[μ,ρ]
- AO one-electron matrix set to Hückel Hamiltonian (eV) as `V` and `T = 0`
These allow running Fermi methods (MP2/CCSD/FCI) over PPP integrals.
"""
function to_fermi_ppp_integrals(system::PPP.MolecularSystem,
                                huckel::PPP.HuckelResult,
                                scf::PPP.SCFResult)
    # Build minimal Fermi Molecule carrying the electron count
    n = size(scf.eigenvectors, 1)
    ndocc = scf.n_occupied
    nvir = n - ndocc
    coords = [a.position for a in system.atoms]
    buf = IOBuffer()
    for r in coords
        println(buf, "H $(r[1]) $(r[2]) $(r[3])")
    end
    molstring = String(take!(buf))
    charge = n - 2*ndocc
    Fermi.Options.set("molstring", molstring)
    Fermi.Options.set("unit", "angstrom")
    Fermi.Options.set("charge", charge)
    Fermi.Options.set("multiplicity", iseven(2*ndocc) ? 1 : 2)
    Fermi.Options.set("basis", "sto-3g")
    molecule = Fermi.Molecule()

    # Hartree ↔ eV
    Ha2eV = 27.21138505
    eV2Ha = 1 / Ha2eV

    # Data
    C = scf.eigenvectors                    # AO→MO (columns = MOs)
    eps_Ha = scf.energies .* eV2Ha          # MO energies (Ha)
    Eref_Ha = scf.total_energy * eV2Ha      # Reference energy (Ha)

    # Construct PPP-backed RHFOrbitals and RHF (immutable types)
    orbs_ppp = Fermi.Orbitals.RHFOrbitals(molecule, "sto-3g", eps_Ha, Eref_Ha, C)
    rhf_ppp = Fermi.HartreeFock.RHF(molecule, Eref_Ha, ndocc, nvir, orbs_ppp, 0.0, 0.0)

    # Integral helpers
    moints = Fermi.Integrals.IntegralHelper{Float64}(orbitals=rhf_ppp.orbitals)
    aoints = Fermi.Integrals.IntegralHelper{Float64}(molecule=rhf_ppp.molecule, eri_type=Fermi.Integrals.Chonky())

    # AO ERIs in ZDO: (μν|ρσ) = δ_{μν} δ_{ρσ} V_raw[μ,ρ]
    V_eV = scf.V_raw
    AOERI = Array{Float64,4}(undef, n, n, n, n)
    @inbounds for μ in 1:n, ν in 1:n, ρ in 1:n, σ in 1:n
        AOERI[μ,ν,ρ,σ] = (μ == ν && ρ == σ) ? (V_eV[μ, ρ] * eV2Ha) : 0.0
    end
    # One-electron AO: set T=0, V = Hückel AO Hamiltonian (eV) in Hartree
    T_AO = zeros(Float64, n, n)
    V_AO = huckel.Hamiltonian .* eV2Ha

    # Populate AO integral cache
    aoints["ERI"] = AOERI
    aoints["T"] = T_AO
    aoints["V"] = V_AO
    # Identity overlap (orthonormal PPP site basis)
    aoints["S"] = Matrix{Float64}(I, n, n)

    return rhf_ppp, moints, aoints
end

# Example usage (uncomment to run manually):
using PPP
sys, huckel, scf = PPP.run_ppp_calculation("molecules/ethene-2C.xyz", PPP.Zhang2011Model())
#sys, huckel, scf = PPP.run_ppp_calculation("molecules/2T-7N.xyz", PPP.Bedogni2024ModelParams())

Ha2eV = 27.21138505

ppp_rhf, ppp_moints, ppp_aoints = to_fermi_ppp_integrals(sys, huckel, scf)

mp2 = Fermi.MollerPlesset.RMP2(ppp_moints, ppp_aoints, Fermi.MollerPlesset.get_rmp2_alg())

ccsd = Fermi.CoupledCluster.RCCSD(ppp_moints, ppp_aoints, Fermi.CoupledCluster.get_rccsd_alg())

fci = Fermi.ConfigurationInteraction.RFCI(ppp_aoints, ppp_rhf, Fermi.ConfigurationInteraction.get_rfci_alg())

# Trying to deduct the nuclear repulsion from FCI, but still massive disagreement with MP2 and CCSD. 
enuc = Fermi.Molecules.nuclear_repulsion(ppp_rhf.molecule.atoms)
println("Nuclear repulsion (FCI correction) (eV): ", enuc * Ha2eV)
fci_energy_vnuc0 = fci.energy - enuc


println("\nTotal Energies (eV):")
println("Huckel: ", huckel.total_energy)
println("PPP: ", scf.total_energy)
println("SCF:  ", ppp_rhf.energy * Ha2eV)
println("MP2:  ", mp2.energy * Ha2eV)
println("CCSD: ", ccsd.energy * Ha2eV)
println("FCI (Vnuc=0):  ", fci_energy_vnuc0 * Ha2eV)

println("\nCorrelation Energies (eV):")
println("MP2:  ", (mp2.energy - ppp_rhf.energy) * Ha2eV)
println("CCSD: ", (ccsd.energy - ppp_rhf.energy) * Ha2eV)
println("FCI:  ", (fci_energy_vnuc0 - ppp_rhf.energy) * Ha2eV)


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

