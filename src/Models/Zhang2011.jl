"""
    Zhang2011Model <: AbstractModel

    Carbon only model, fitted to stretched ethene molecules with CASPT2, then
    solved in PPP at full CI as closed solutions (2 sites!). 
    
    Zhang et al, J.Chem.Phys. 134, 024114 (2011). 

All energies are in eV, distances in Angstroms, inferred from figure 1.

Fields (optimizable parameters):
- `t_ii_E`: Base site energy (eV)
- `t_ii_I`: Induction correction denominator offset
- `t_ij_E`: Hopping amplitude (eV)
- `t_ij_B`: Hopping exponential decay (1/Å)
- `γ_ii_E`: Base Hubbard U (eV)
- `γ_ii_I`: Induction amplitude (eV·Å²)
- `γ_ii_I_B`: Induction offset (Å)
- `γ_ij_A`: Ohno numerator (eV·Å)
- `γ_ij_B`: Ohno denominator offset (Å²)
"""
@kwdef struct Zhang2011Model <: AbstractModel
    # t_ii parameters
    t_ii_E::Float64 = -9.13597
    t_ii_I::Float64 = 0.15076
    # t_ij parameters  
    t_ij_E::Float64 = -28.07749
    t_ij_B::Float64 = 1.65878
    # γ_ii parameters
    γ_ii_E::Float64 = 11.56718
    γ_ii_I::Float64 = 0.16640
    γ_ii_I_B::Float64 = 2.30448
    # γ_ij parameters
    γ_ij_E::Float64 = 21.88221
    γ_ij_I::Float64 = 7.05909
    # Connectivity cutoff
    cutoff::Float64 = 1.65
end

function t_ii(model::Zhang2011Model, system, i)
    t_ii = model.t_ii_E
    for j in 1:length(system.atoms)
        if j == i
            continue
        end
        Rij = norm(system.atoms[i].position - system.atoms[j].position)
        t_ii -= (Rij)/sqrt(Rij^2 - model.t_ii_I) * γ_ij(model, system, i, j) 
    end
    return t_ii
end

function t_ij(model::Zhang2011Model, system, i, j)
    Rij = norm(system.atoms[i].position - system.atoms[j].position)
    return model.t_ij_E * exp(-model.t_ij_B * Rij) 
end

function γ_ii(model::Zhang2011Model, system, i)
    atom = system.atoms[i]
    if atom.symbol == :C
        # for all neighbours, calculate an induction contribution
        neighbours = [j for j in 1:length(system.atoms) if system.connectivity[i,j] == 1]
        if isempty(neighbours)
            return model.γ_ii_E
        end
        Rijs = [norm(system.atoms[i].position - system.atoms[j].position) for j in neighbours]
        induction = sum(model.γ_ii_I / (Rij - model.γ_ii_I_B)^2 for Rij in Rijs)
        return model.γ_ii_E - induction
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
end

function γ_ij(model::Zhang2011Model, system, i, j)
    Rij = norm(system.atoms[i].position - system.atoms[j].position)
    return model.γ_ij_E / sqrt(Rij^2 + model.γ_ij_I)
end
