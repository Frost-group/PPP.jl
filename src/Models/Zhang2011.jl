"""
    Zhang2011Model <: AbstractModel

    Carbon only model, fitted to stretched ethene molecules with CASPT2, then
    solved in PPP at full CI as closed solutions (2 sites!). 
    
    Zhang et al, J.Chem.Phys. 134, 024114 (2011). 

All energies are in eV, distances in Angstroms, inferred from figure 1.
"""
@kwdef struct Zhang2011Model <: AbstractModel
    # I've just left this one as magic numbers for now, as I think this keeps it closer to the paper. 
    # I suppose I'll also need a generic, numbers packed away, version of this for fitting. 

    Z_EFF_C=0.0 # OK, this is some weird holdover from the old interface 
    N_C=0.0 
end

function t_ii(model::Zhang2011Model, system, i)
    t_ii = -9.13597
    for j in 1:length(system.atoms)
        if j == i
            continue
        end

        Rij = norm(system.atoms[i].position - system.atoms[j].position)
        t_ii -= (Rij)/sqrt(Rij^2-0.15076) * γ_ij(model, system, i, j) 
    end
    return t_ii
end

function t_ij(model::Zhang2011Model, system, i, j)
    Rij = norm(system.atoms[i].position - system.atoms[j].position)
    return -28.07749 * exp(-1.65878*Rij) 
end

function γ_ii(model::Zhang2011Model, system, i)
    atom = system.atoms[i]
    if atom.symbol == :C
        # for all neighbours, calculate an induction contribution
        neighbours=[j for j in 1:length(system.atoms) if system.connectivity[i,j] == 1]
        if isempty(neighbours)
            return 11.56718  # No neighbours, no induction contribution
        end
        Rijs = [norm(system.atoms[i].position - system.atoms[j].position) for j in neighbours]
        induction=sum(0.16640/(Rij-2.30448)^2 for Rij in Rijs) # could write this iwth a splat operator, but I think clearer this way 
        return 11.56718-induction
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
end

function γ_ij(model::Zhang2011Model, system, i, j)
    Rij = norm(system.atoms[i].position - system.atoms[j].position)
    return 21.88221 / sqrt(Rij^2 + 7.05909)
end
