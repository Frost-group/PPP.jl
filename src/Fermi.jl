# Interface with Fermi.jl for all that CI goodness

using PPP, Fermi

system, huckel_result, scf_result = run_ppp_calculation("molecules/2T-7N.xyz", Bedogni2024ModelParams())

function FermiRHF(system, scf_result)
    # 1st argument: molecule object
    molstringA="""
H 0 0 1.0
H 0 0 2.0
H 0 0 3.0
H 0 0 4.0
H 0 0 5.0
H 0 0 6.0
H 0 0 7.0
H 0 0 8.0
H 0 0 9.0
H 0 0 10.0
H 0 0 11.0
H 0 0 12.0
H 0 0 13.0
"""
# I'm trying to hack in a 13 state basis, as I think Fermi.jl uses the molecule.Nα and molecule.Nβ as indices

    molecule = Fermi.Molecule(
        molstring=molstringA,
        unit=:angstrom,
        charge=0,
        multiplicity=2
    )

    # 2nd argument: energy
    energy = scf_result.total_energy
    # 3rd argument: Number of doubly occ orbitals
    ndocc = scf_result.n_occupied
    # 4th argument: Number of virtual orbitals
    nvir = size(scf_result.eigenvectors,1) - ndocc
    # 5th argument: RHFOrbitals object
    orbitals = Fermi.Orbitals.RHFOrbitals(molecule, "PPP", [scf_result.F[a,a] for a in 1:size(scf_result.F,1)], energy, scf_result.K)
    # 6th and 7th, convergency parameters. We will skip those for now.

    return Fermi.HartreeFock.RHF(molecule, energy, ndocc, nvir, orbitals, 0.0, 0.0)
end

# Generate Fermi.jl objects, and try and do a calculation with them 
wfn=FermiRHF(system,scf_result)
@energy wfn => mp2 

# attempts to start calculation, but fails with error: 
# ERROR: LoadError: DimensionMismatch: non-matching sizes in contracted dimensions
# Stacktrace:
#   [1] dimcheck_tensorcontract
#     @ ~/.julia/packages/TensorOperations/RaU9j/src/implementation/abstractarray.jl:167 [inlined]
#   [2] tensorcontract!
#     @ ~/.julia/packages/TensorOperations/RaU9j/src/implementation/abstractarray.jl:145 [inlined]
#   [3] tensorcontract!
#     @ ~/.julia/packages/TensorOperations/RaU9j/src/implementation/abstractarray.jl:145 [inlined]
#   [4] tensorcontract!
#     @ ~/.julia/packages/TensorOperations/RaU9j/src/implementation/abstractarray.jl:145 [inlined]
