# Interface with Fermi.jl for PPP -> Fermi wavefunction handoff

using PPP
using Fermi
using LinearAlgebra

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

    # MO data from PPP SCF
    mo_energies = scf_result.energies
    mo_coefficients = scf_result.eigenvectors

    # Fermi.Orbitals.RHFOrbitals expects:
    #   (molecule::Molecule, basis::String, eps::Vector{Float64}, sd_energy::Float64, C::Matrix{Float64})
    orbitals = Fermi.Orbitals.RHFOrbitals(
        molecule,
        "sto-3g",
        mo_energies,
        scf_result.total_energy,
        mo_coefficients,
    )

    hf_energy = scf_result.total_energy
    return Fermi.HartreeFock.RHF(molecule, hf_energy, ndocc, nvir, orbitals, 0.0, 0.0)
end

# Example usage (uncomment to run manually):
using PPP
sys, huckel, scf = PPP.run_ppp_calculation("molecules/pyrrole.xyz", PPP.Bedogni2024ModelParams())
wfn = to_fermi_rhf(sys, scf)
@energy wfn => mp2

# CCSD: create integral helpers explicitly (does not dispatch on RHF)
moints = Fermi.Integrals.IntegralHelper{Float64}(orbitals=wfn.orbitals)
aoints = Fermi.Integrals.IntegralHelper{Float64}(molecule=wfn.molecule)
ccsd = Fermi.CoupledCluster.RCCSD(moints, aoints, Fermi.CoupledCluster.get_rccsd_alg())
display(ccsd)

# FCI: API expects (aoints, rhf)
fci = Fermi.ConfigurationInteraction.RFCI(aoints, wfn, Fermi.ConfigurationInteraction.get_rfci_alg())
display(fci)
