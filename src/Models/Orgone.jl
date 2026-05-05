# Orgone.jl - the pseudoscientific concept of anti-entropic life force
#
# ALSO - a PPP model which calculates the the electron repulsion gamma (γ) via 3 overlapping
# Gaussian pseudo-densities. 
# The most diffuse of these Gaussians then 'breathes' (expand/contract) depending on the
# local induction from neighbouring atoms, dynamically modifying the chemical 'hardness'. 
# The pseudo-densities can become negative, allowing for Löwdin orthogonalisation holes etc. 
# Also added a very basic form for hopping integrals to vary with the normal plane, so that
# this could hopefully be used for a nonadiabatic dynamics simulation. 
#
# So an initial attempt to make pairwise PPP many-body.
#
# Only tested against the ~10 rattled cycobutadiene structures I had on my laptop - may well
# still be horrible bugs! 
# Also apologies for the naming conventions - it was late and I was listening to Kate Bush
# on loop. - Jarv

using TOML
using SpecialFunctions
using LinearAlgebra

# We assume q, ε₀ are available (from PPP.jl exports/includes)

# #I still dream of Orgonon 
@kwdef struct OrgoneAtomParams
    ζ::Vector{Float64} # Exponents (inverse width) of Gaussians for PPP 
    c::Vector{Float64} # Weights
    χ::Float64         # Polarizing power
    α::Float64         # Susceptibility
    Z_eff::Float64 = 3.25 # Slater Z_eff for hopping
    n_star::Float64 = 2.0 # Principal quantum number effective
    k_t::Float64 = 1.0    # Hopping scaling factor
    IP::Float64 = 11.16   # Ionization potential (site energy)
    pi_electrons::Int = 1 # Number of pi electrons donated
end

@kwdef struct OrgoneModel <: AbstractModel
    params::Dict{Symbol,OrgoneAtomParams}
    ohnoconstant::Float64 = q / (4π * ε₀ * 1e-10) # e²/4πε₀ in eV⋅Å
    cutoff::Float64 = 1.4 # Connectivity cutoff
end

# Save and load; TOML
# #I can't hide you from the government
function save_orgone(model::OrgoneModel, path::String)
    dict_out = Dict{String,Any}()
    for (k, v) in model.params
        dict_out[string(k)] = Dict(
            "zeta" => v.ζ,
            "c" => v.c,
            "chi" => v.χ,
            "alpha" => v.α,
            "Z_eff" => v.Z_eff,
            "n_star" => v.n_star,
            "k_t" => v.k_t,
            "IP" => v.IP,
            "pi_electrons" => v.pi_electrons
        )
    end
    open(path, "w") do io
        TOML.print(io, dict_out)
    end
end

function load_orgone(path::String)::OrgoneModel
    dict_in = TOML.parsefile(path)
    params = Dict{Symbol,OrgoneAtomParams}()
    for (k, v) in dict_in
        params[Symbol(k)] = OrgoneAtomParams(
            ζ=Vector{Float64}(v["zeta"]),
            c=Vector{Float64}(v["c"]),
            χ=Float64(v["chi"]),
            α=Float64(v["alpha"]),
            Z_eff=Float64(get(v, "Z_eff", 3.25)),
            n_star=Float64(get(v, "n_star", 2.0)),
            k_t=Float64(get(v, "k_t", 1.0)),
            IP=Float64(get(v, "IP", 11.16)),
            pi_electrons=Int(get(v, "pi_electrons", 1))
        )
    end
    return OrgoneModel(params=params)
end

# #Every time it rains
function default_orgone_model(N::Int=3)::OrgoneModel
    params = Dict{Symbol,OrgoneAtomParams}()

    # Carbon: from SAMIN optimisation against rattled cyclobutadiene (EOM-CCSD(T) reference)
    params[:C] = OrgoneAtomParams(
        ζ=[3.42, 0.74, 0.79][1:N],
        c=[0.01, 0.52, 0.47][1:N],
        χ=3.87,
        α=0.20,
        Z_eff=3.15,
        n_star=2.0,
        k_t=0.59,
        IP=11.16,
        pi_electrons=1
    )

    # Heteroatoms: ζ scaled by (Z_eff_X/Z_eff_C), χ scaled by electronegativity,
    # IPs from standard PPP valence-state values (Hinze-Jaffé)

    # Pyridine-like N (2 bonds, lone pair in plane, 1 π electron)
    params[:N1] = OrgoneAtomParams(
        ζ=[4.24, 0.92, 0.98][1:N],
        c=[0.01, 0.52, 0.47][1:N],
        χ=4.8,
        α=0.20,
        Z_eff=3.90,
        n_star=2.0,
        k_t=0.59,
        IP=14.12,
        pi_electrons=1
    )

    # Pyrrole-like N (3 bonds, lone pair in π system, 2 π electrons)
    params[:N2] = OrgoneAtomParams(
        ζ=[4.02, 0.87, 0.93][1:N],
        c=[0.01, 0.52, 0.47][1:N],
        χ=4.8,
        α=0.20,
        Z_eff=3.90,
        n_star=2.0,
        k_t=0.59,
        IP=28.71,
        pi_electrons=2
    )

    # Carbonyl O (1 bond, 1 π electron)
    params[:O1] = OrgoneAtomParams(
        ζ=[4.94, 1.07, 1.14][1:N],
        c=[0.01, 0.52, 0.47][1:N],
        χ=5.6,
        α=0.20,
        Z_eff=4.55,
        n_star=2.0,
        k_t=0.59,
        IP=17.70,
        pi_electrons=1
    )

    # Ether/Furan O (2 bonds, lone pair in π system, 2 π electrons)
    params[:O2] = OrgoneAtomParams(
        ζ=[4.70, 1.02, 1.09][1:N],
        c=[0.01, 0.52, 0.47][1:N],
        χ=5.6,
        α=0.20,
        Z_eff=4.55,
        n_star=2.0,
        k_t=0.59,
        IP=34.12,
        pi_electrons=2
    )

    # Normalize c to sum to 1
    for (k, v) in params
        if N > 1
            v.c[end] = 1.0 - sum(v.c[1:end-1])
        else
            v.c[1] = 1.0
        end
    end

    return OrgoneModel(params=params)
end

# I'm not sure if it's just the IPA, but I'm starting to think like a chemist...
# #Ooh, I just know that something good is gonna happen
function get_breathed_zeta(model::OrgoneModel, system, i::Int)
    atom = system.atoms[i]
    type = get_atomtype(atom)

    if !haskey(model.params, type)
        error("Atom type $type not found in Orgone parameters.")
    end

    base_params = model.params[type]
    ζ_breathed = copy(base_params.ζ)

    # Calculate induction from local cluster
    env_sum = 0.0
    for j in 1:length(system.atoms)
        if system.connectivity[i, j] == 1
            neighbor_type = get_atomtype(system.atoms[j])
            if haskey(model.params, neighbor_type)
                env_sum += model.params[neighbor_type].χ
            end
        end
    end

    # Find the most diffuse Gaussian (the one with the smallest exponent ζ)
    # The diffuse valence tail is the part of the cloud most easily deformed by neighbors.
    diffuse_idx = argmin(ζ_breathed)

    # Apply the local cluster induction to contract/expand the diffuse tail
    #  thus induction --> contract the diffuse tail --> increase the Hubbard U (chemical hardness)
    ζ_breathed[diffuse_idx] += base_params.α * env_sum

    # and now transformed by the local chemical environment, it will be used to calculate gamma...
    return ζ_breathed, base_params.c
end

# #what made it special, made it dangerous
function γ_ij(model::OrgoneModel, system, i::Int, j::Int)
    r_ij = norm(system.atoms[i].position - system.atoms[j].position)

    ζ_i, c_i = get_breathed_zeta(model, system, i) # induction effect...
    ζ_j, c_j = get_breathed_zeta(model, system, j)

    val = 0.0
    N = length(ζ_i) # Assuming same N for all atoms

    for k in 1:N # sum over overlapping Gaussians in aux basis
        for l in 1:N
            α_kl = sqrt((ζ_i[k] * ζ_j[l]) / (ζ_i[k] + ζ_j[l]))
            if r_ij < 1e-8 # limit r -> 0
                val += c_i[k] * c_j[l] * 2 * α_kl / sqrt(π)
            else
                val += c_i[k] * c_j[l] * erf(α_kl * r_ij) / r_ij
            end
        end
    end

    return val * model.ohnoconstant
end

function γ_ii(model::OrgoneModel, system, i::Int)
    return γ_ij(model, system, i, i)
end

# Compute local sp2 plane normal for an atom
#  man its a mess, walking neighbour thingy etc. ; not very high performance or robust
function get_normal(system, i::Int)
    atom = system.atoms[i]
    neighbors = []
    for j in 1:length(system.atoms)
        if system.connectivity[i, j] == 1
            push!(neighbors, system.atoms[j].position)
        end
    end

    if length(neighbors) >= 3
        # Symmetric normal for sp2 centers: fixes missing 3rd atom force
        v1, v2, v3 = neighbors[1] - atom.position, neighbors[2] - atom.position, neighbors[3] - atom.position
        n1, n2, n3 = cross(v1, v2), cross(v2, v3), cross(v3, v1)
        # Align them to prevent cancellation due to arbitrary index ordering
        n = n1 + sign(dot(n1, n2))*n2 + sign(dot(n1, n3))*n3
        if norm(n) > 1e-8
            return normalize(n)
        end
    elseif length(neighbors) == 2
        v1 = neighbors[1] - atom.position
        v2 = neighbors[2] - atom.position
        n = cross(v1, v2)
        if norm(n) > 1e-8
            return normalize(n)
        end
    elseif length(neighbors) == 1
        # Terminal atom (e.g., carbonyl oxygen): inherit normal from its neighbor
        j = findfirst(x -> x == 1, system.connectivity[i, :])
        if j !== nothing && j != i
            # To avoid infinite recursion, only do a simple neighbor check here
            neighbor = system.atoms[j]
            n_neighbors = []
            for k in 1:length(system.atoms)
                if system.connectivity[j, k] == 1 && k != i
                    push!(n_neighbors, system.atoms[k].position)
                end
            end
            if length(n_neighbors) >= 1
                v1 = atom.position - neighbor.position
                v2 = n_neighbors[1] - neighbor.position
                n = cross(v1, v2)
                if norm(n) > 1e-8
                    return normalize(n)
                end
            end
        end
    end

    # Fallback if structure is perfectly linear or isolated
    return [0.0, 0.0, 1.0]
end

# hacked off Jorner model to provide twist effect of t_ij
function t_ij(model::OrgoneModel, system, i::Int, j::Int)
    r_ij = norm(system.atoms[i].position - system.atoms[j].position)

    type_i = get_atomtype(system.atoms[i])
    type_j = get_atomtype(system.atoms[j])

    p_i = model.params[type_i]
    p_j = model.params[type_j]

    # Slater exponents Z_eff / n_star 
    exp_i = p_i.Z_eff / p_i.n_star
    exp_j = p_j.Z_eff / p_j.n_star
    k_t = 0.5 * (p_i.k_t + p_j.k_t) # Averaged scaling factor

    r_bohr = Angstrom2Bohr * r_ij

    # Use slater_grad from PPP.jl namespace
    overlap_grad = slater_grad(r_bohr, p_i.n_star, p_j.n_star, exp_i, exp_j, 0.01)
    beta = overlap_grad / r_bohr

    # Directional Cosine Math for twists: cos(theta) = n_i ⋅ n_j
    n_i = get_normal(system, i)
    n_j = get_normal(system, j)
    # Removed abs() to eliminate 90-degree cusp and preserve Berry phase
    dir_cos = dot(n_i, n_j) 

    # Note: dir_cos can now be negative. To prevent Mobius frustration in rings, 
    # a reference sign matrix should ideally be multiplied here for NAMD initialization.

    return k_t * beta * Ha2eV * dir_cos
end

# site energy by IP lookup
function t_ii(model::OrgoneModel, system, i::Int)
    type = get_atomtype(system.atoms[i])
    p = model.params[type]
    
    # Calculate induction to shift IP
    env_sum = 0.0
    for j in 1:length(system.atoms)
        if system.connectivity[i, j] == 1
            neighbor_type = get_atomtype(system.atoms[j])
            if haskey(model.params, neighbor_type)
                env_sum += model.params[neighbor_type].χ
            end
        end
    end
    
    # Contracted cloud (harder) = higher Ionization Potential
    return -(p.IP + p.α * env_sum)
end

# was defaulting to 1, which worekd for carbon...
function get_pi_electrons(model::OrgoneModel, atom)
    type = get_atomtype(atom)
    return model.params[type].pi_electrons
end

function update_atom_params(model::OrgoneModel, atom)
    type = get_atomtype(atom)
    p = model.params[type]
    return Atom(atom.symbol, atom.position, atom.n_bonds,
        p.pi_electrons, -p.IP, 0.0, p.Z_eff, p.n_star)
end

