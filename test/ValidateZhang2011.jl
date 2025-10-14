using Test
using PPP
using Printf
using StaticArrays

@testset "Validate Zhang2011Model with ethene" begin
    system, huckel_result, scf_result = run_ppp_calculation("../molecules/ethene.xyz", Zhang2011Model())
    
    @test scf_result.converged
    @test length(system.atoms) == 2
    @test system.n_electrons == 2       
    # Verify atom composition
    carbon_count = count(atom -> atom.symbol == :C, system.atoms)
    @test carbon_count == 2

    @test huckel_result.energies[1] ≈ -10.0450931042920 atol=1e-3
    @test huckel_result.energies[2] ≈ 0.0 atol=1e-3

    singlet, triplet = run_cis_calculation(system, scf_result)

    @test singlet.energies[1] - triplet.energies[1] ≈ 3.45 atol=0.1
end

@testset "Repro Fig 1 from Zhang2011Model" begin
    system, huckel_result, scf_result = run_ppp_calculation("../molecules/ethene.xyz", Zhang2011Model())

    # Local helper function to deal with kruft of immutability in PPP
    function update_system_atom_position(system::PPP.MolecularSystem, atom_idx::Int, new_position::SVector{3,Float64})
        new_atoms = copy(system.atoms)
        old_atom = new_atoms[atom_idx]
        new_atoms[atom_idx] = PPP.Atom(
            old_atom.symbol, 
            new_position, 
            old_atom.n_bonds, 
            old_atom.nz, 
            old_atom.site_energy, 
            old_atom.Hubbard_U, 
            old_atom.Z_eff, 
            old_atom.n_number
        )
        
        # Recalculate connectivity based on new positions
        new_positions = [atom.position for atom in new_atoms]
        new_connectivity = PPP.calculate_connectivity(new_positions, Zhang2011Model())
        
        return PPP.MolecularSystem(new_atoms, new_connectivity, system.n_electrons)
    end

    t_iis = []
    t_ijs = []
    gamma_iis = []
    gamma_ijs = []
    rs = 1.30:0.01:1.60

    for r in rs
        new_position = system.atoms[1].position + SVector{3,Float64}(r, 0.0, 0.0)
        system = update_system_atom_position(system, 2, new_position)

        append!(t_iis, t_ii(Zhang2011Model(), system, 1))
        append!(t_ijs, t_ij(Zhang2011Model(), system, 1, 2))
        append!(gamma_iis, γ_ii(Zhang2011Model(), system, 1))
        append!(gamma_ijs, γ_ij(Zhang2011Model(), system, 1, 2))

    end

    for d in zip(rs,t_iis, t_ijs, gamma_iis, gamma_ijs)
        @printf("%f %f %f %f %f\n", d...)
    end
end
