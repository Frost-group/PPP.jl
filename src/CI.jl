using LinearAlgebra

# Types for excitations
abstract type ExcitationType end
struct SingleExcitation <: ExcitationType end
struct DoubleExcitation <: ExcitationType end


# Configuration represents |Φᵢ⟩ in CI expansion
struct Configuration
    type::ExcitationType
    from_orbitals::Vector{Int}  # i,j indices
    to_orbitals::Vector{Int}    # a,b indices
end

# Slater determinant for FCI
struct SlaterDeterminant
    occupied_orbitals::Vector{Int}
    n_electrons::Int
end

# Results types
struct CIResult
    method::String
    energies::Vector{Float64}        # Excitation energies in eV
    eigenvectors::Matrix{Float64}    # CI coefficients
    configurations::Vector{Configuration}
    dominant_configurations::Vector{Int}
    oscillator_strengths::Vector{Float64}
    characters::Vector{String}
end

struct FCIResult
    energies::Vector{Float64}
    eigenvectors::Matrix{Float64}
    determinants::Vector{SlaterDeterminant}
    oscillator_strengths::Vector{Float64}
    state_compositions::Vector{Dict{String,Float64}}
end

# Generate configurations for CI expansion
function generate_single_excitations(n_orbs::Int, n_occ::Int)
    configs = Configuration[]
    for i in 1:n_occ, a in (n_occ+1):n_orbs
        push!(configs, Configuration(SingleExcitation(), [i], [a]))
    end
    return configs
end

function generate_double_excitations(n_orbs::Int, n_occ::Int)
    configs = Configuration[]
    for i in 1:n_occ, j in (i+1):n_occ
        for a in (n_occ+1):n_orbs, b in (a+1):n_orbs
            push!(configs, Configuration(DoubleExcitation(), [i,j], [a,b]))
        end
    end
    return configs
end

# Matrix elements using PPP Hamiltonian
function calculate_ci_matrix_element(config1::Configuration, config2::Configuration, 
                                   scf_result::SCFResult, system::MolecularSystem)
    if config1.type == SingleExcitation && config2.type == SingleExcitation
        i, a = config1.from_orbitals[1], config1.to_orbitals[1]
        j, b = config2.from_orbitals[1], config2.to_orbitals[1]
        
        if i == j && a == b
            # Diagonal element: εₐ - εᵢ + ⟨ia||ia⟩
            return (scf_result.energies[a] - scf_result.energies[i]) + 
                   (2 * scf_result.K[a,i] - scf_result.K[a,a] - scf_result.K[i,i])
        elseif i == j || a == b
            # Off-diagonal element: ⟨ia||jb⟩
            return scf_result.K[max(a,b), min(i,j)]
        end
    elseif config1.type == DoubleExcitation && config2.type == DoubleExcitation
        i, j = config1.from_orbitals
        a, b = config1.to_orbitals
        k, l = config2.from_orbitals
        c, d = config2.to_orbitals
        
        if i == k && j == l && a == c && b == d
            # Diagonal element
            return (scf_result.energies[a] + scf_result.energies[b] - 
                   scf_result.energies[i] - scf_result.energies[j]) +
                   scf_result.K[a,b] - scf_result.K[i,j]
        end
    end
    return 0.0
end

# Transform PPP two-electron integrals from AO to MO basis
function transform_two_electron_integral(i::Int, a::Int, j::Int, b::Int, scf_result::SCFResult)
    C = scf_result.eigenvectors
    V = scf_result.V_raw
    return sum(C[μ,i] * C[μ,a] * V[μ,ν] * C[ν,j] * C[ν,b]
              for μ in axes(C,1), ν in axes(C,1))
end

function calculate_two_electron_term(i::Int, a::Int, scf_result::SCFResult)
    C = scf_result.eigenvectors
    V = scf_result.V_raw

    n_sites = size(V, 1)
    term = 0.0
    
    # Calculate ⟨ia|ia⟩ - ⟨ii|aa⟩
    for μ in 1:n_sites, ν in 1:n_sites
        for λ in 1:n_sites, σ in 1:n_sites
            # Direct term ⟨ia|ia⟩
            term += C[μ,i] * C[ν,a] * V[μ,ν] * C[λ,i] * C[σ,a]
            
            # Exchange term -⟨ii|aa⟩
            term -= C[μ,i] * C[ν,i] * V[μ,ν] * C[λ,a] * C[σ,a]
        end
    end
    
    return term
end

# Calculate oscillator strengths using transition dipoles
function calculate_oscillator_strength(config::Configuration, scf_result::SCFResult, 
                                     system::MolecularSystem)
    config.type != SingleExcitation && return 0.0
    
    i, a = config.from_orbitals[1], config.to_orbitals[1]
    C = scf_result.eigenvectors
    
    μ = sum(system.atoms[μ].position .* C[μ,a] * C[μ,i] for μ in eachindex(system.atoms))
    ΔE = (scf_result.energies[a] - scf_result.energies[i]) * 0.0367493  # eV to a.u.
    
    return (2/3) * ΔE * sum(abs2, μ)
end

"""
Calculate CIS matrix element ⟨Φᵢᵃ|H|Φⱼᵇ⟩ using PPP Hamiltonian
"""
function calculate_cis_matrix_element(i::Int, a::Int, j::Int, b::Int, scf_result::SCFResult; singlet::Bool=true)
    if singlet==true 
        if i == j && a == b
            # Diagonal element
            orbital_energy_diff = scf_result.energies[a] - scf_result.energies[i]
            return orbital_energy_diff + 2 * transform_two_electron_integral(i,a,j,b, scf_result) - transform_two_electron_integral(i,j,a,b, scf_result)
        else
            return 2 * transform_two_electron_integral(i,a,j,b, scf_result) - transform_two_electron_integral(i,j,a,b, scf_result)
        end
    else
        if i == j && a == b
            # Diagonal element
            orbital_energy_diff = scf_result.energies[a] - scf_result.energies[i]
            return orbital_energy_diff - transform_two_electron_integral(i,j,a,b, scf_result)
        else
            return - transform_two_electron_integral(i,j,a,b, scf_result)
        end
    end
end


"""
Calculate CISD matrix element ⟨Φᵢⱼᵃᵇ|H|Φₖₗᶜᵈ⟩ using PPP Hamiltonian
"""
function calculate_cisd_matrix_element(i::Int, j::Int, a::Int, b::Int, 
                                     k::Int, l::Int, c::Int, d::Int, scf_result::SCFResult)
    if i == k && j == l && a == c && b == d
        # Diagonal: ⟨Φᵢⱼᵃᵇ|H|Φᵢⱼᵃᵇ⟩
        return (scf_result.energies[a] + scf_result.energies[b] - 
                scf_result.energies[i] - scf_result.energies[j]) +
                scf_result.K[a,b] - scf_result.K[i,j]
    end
    return 0.0
end

function run_cis_calculation(system::MolecularSystem, scf_result::SCFResult)
    n_orbs = size(scf_result.eigenvectors, 2)
    n_occ = system.n_electrons ÷ 2
    
    println("\nStarting CIS calculation:")
    println("Number of occupied orbitals: $n_occ")
    println("Total number of orbitals: $n_orbs")
    
    # Print orbital energies
    println("\nOrbital energies:")
    for i in 1:n_orbs
        @printf("ε[%d] = %8.3f eV\n", i, scf_result.energies[i])
    end
    
    # Generate all single excitations
    configs = generate_single_excitations(n_orbs, n_occ)
    n_configs = length(configs)
    println("\nNumber of configurations: $n_configs")
    
    # Build CIS matrix
    H_S = zeros(n_configs, n_configs)
    H_T = zeros(n_configs, n_configs)

    for p in 1:n_configs, q in 1:p
        i, a = configs[p].from_orbitals[1], configs[p].to_orbitals[1]
        j, b = configs[q].from_orbitals[1], configs[q].to_orbitals[1]
        H_S[p,q] = calculate_cis_matrix_element(i,a,j,b,scf_result; singlet=true)
        H_S[q,p] = H_S[p,q]
        H_T[p,q] = calculate_cis_matrix_element(i,a,j,b,scf_result; singlet=false)
        H_T[q,p] = H_T[p,q]
    end
    
    # Print CIS matrix
    println("\nCIS matrix:")
    for i in 1:min(5, n_configs)
        for j in 1:min(5, n_configs)
            @printf(" %8.3f", H_S[i,j])
        end
        println()
    end

    println("\nCIS matrix:")
    for i in 1:min(5, n_configs)
        for j in 1:min(5, n_configs)
            @printf(" %8.3f", H_T[i,j])
        end
        println()
    end
    
    # Solve eigenvalue problem
    E_S, C_S = eigen(Symmetric(H_S))
    E_T, C_T = eigen(Symmetric(H_T))
    
    
    # Calculate oscillator strengths
    f = zeros(n_configs)
    for state in 1:n_configs
        μ = zeros(2)  # transition dipole in x,y
        for (idx, config) in enumerate(configs)
            # Sum over all configurations contributing to this state
            for μ_idx in 1:length(system.atoms)
                pos = system.atoms[μ_idx].position
                i, a = config.from_orbitals[1], config.to_orbitals[1]
                μ .+= C_S[idx,state] * pos .* (scf_result.eigenvectors[μ_idx,a] * 
                                           scf_result.eigenvectors[μ_idx,i])
            end
        end
        ΔE = E_S[state] * 0.0367493  # eV to a.u.
        f[state] = (2/3) * ΔE * sum(abs2, μ)
    end
    
    return (CIResult(
        "CIS",
        E_S,
        C_S,
        configs,
        [argmax(abs.(C_S[:,i])) for i in 1:n_configs],
        f,
        ["$(c.from_orbitals[1])→$(c.to_orbitals[1])" for c in configs]
        ),
        CIResult(
            "CIS",
            E_T,
            C_T,
            configs,
            [argmax(abs.(C_T[:,i])) for i in 1:n_configs],
            f,
            ["$(c.from_orbitals[1])→$(c.to_orbitals[1])" for c in configs]
        ))
end


export run_cis_calculation
export SingleExcitation, DoubleExcitation, Configuration, CIResult