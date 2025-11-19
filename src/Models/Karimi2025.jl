"""
    Karimi2025Model <: AbstractModel

    Carbon only model, fitted to corone experimental data (?) 
    
    Karimi et al, ArXiv, 2508.18963v1

    Initial stab at parameters: a bit unclear from the paper specifically what is going on. 

    Huckel style (i.e. no distance dependence hopping), but dist. dep for γ_ij. ?
"""
@kwdef struct Karimi2025Model <: AbstractModel
    Z_EFF_C=0.0 # OK, this is some weird holdover from the old interface 
    N_C=0.0 
    cutoff=1.65 # Cut-off for connectivity matrix, in Angstroms
end

function t_ii(model::Karimi2025Model, system, i)
    t_ii = 0.0 
    return t_ii
end

function t_ij(model::Karimi2025Model, system, i, j)
    Rij = norm(system.atoms[i].position - system.atoms[j].position)
    return 2.4 
end

function γ_ii(model::Karimi2025Model, system, i)
    atom = system.atoms[i]
    if atom.symbol == :C
        return 8.0  # only carbon, 
    else
        error("Unsupported atom type: $(atom.symbol)")
    end
end

function γ_ij(model::Karimi2025Model, system, i, j)
    Rij = norm(system.atoms[i].position - system.atoms[j].position)
    r_0 = 0.25 # Angstroms, pass-over from OhNo semi-classical to Hubbard U. 
    return γ_ii(model, system, i) / sqrt(1.0 + (Rij/r_0)^2)
end
