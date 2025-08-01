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

# Generates all single and double excitations
function generate_excitations(n_orbs::Int, n_occ::Int)
    configs = Configuration[]
    singles = generate_single_excitations(n_orbs, n_occ)
    doubles = generate_double_excitations(n_orbs, n_occ)
    configs = vcat(singles,doubles)
end

# Transform PPP two-electron integrals from AO to MO basis
function transform_two_electron_integral(i::Int, a::Int, j::Int, b::Int, scf_result::SCFResult)
    C = scf_result.eigenvectors
    V = scf_result.V_raw
    return sum(C[μ,i] * C[μ,a] * V[μ,ν] * C[ν,j] * C[ν,b]
              for μ in axes(C,1), ν in axes(C,1))
end

"""
Calculate CIS matrix element ⟨Φᵢᵃ|H|Φⱼᵇ⟩ using PPP Hamiltonian
"""
# Matrix elements obtained from Table 4.1 in Szabo & Ostlund: Modern Quantum Chemistry and 10.1016/0584-8539(72)80159-4.
function calculate_cis_matrix_element(i::Int, a::Int, j::Int, b::Int, scf_result::SCFResult; singlet::Bool=true)
    ΔE = scf_result.energies[a] - scf_result.energies[i]
    if singlet==true 
        # Singlet matrix elements: (ε_a-ε_i)*δ_ab*δ_ij + 2*(ia|jb) - (ij|ab)
        if i == j && a == b
            # Singlet Diagonal element
            return ΔE+ 2 * transform_two_electron_integral(i,a,j,b, scf_result) - transform_two_electron_integral(i,j,a,b, scf_result)
        else
            # Singlet Off-Diagonal element (one or more kronecker delta = 0)
            return 2 * transform_two_electron_integral(i,a,j,b, scf_result) - transform_two_electron_integral(i,j,a,b, scf_result)
        end
    else
        # Triplet matrix elements: (ε_a-ε_i)*δ_ab*δ_ij - (ij|ab)
        if i == j && a == b
            # Triplet Diagonal element
            return ΔE - transform_two_electron_integral(i,j,a,b, scf_result)
        else
            # Triplet Off-Diagonal element (one or more kronecker delta = 0)
            return - transform_two_electron_integral(i,j,a,b, scf_result)
        end
    end
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
    println("\nFirst 10-block of Singlet CIS matrix:")
    for i in 1:min(10, n_configs)
        for j in 1:min(10, n_configs)
            @printf(" %8.3f", H_S[i,j])
        end
        println()
    end

    println("\nFirst 10-block of Triplet CIS matrix:")
    for i in 1:min(10, n_configs)
        for j in 1:min(10, n_configs)
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
            i, a = config.from_orbitals[1], config.to_orbitals[1]
            for μ_idx in 1:length(system.atoms)
                pos = system.atoms[μ_idx].position[1:2]
                μ .+= C_S[idx,state] * pos * (scf_result.eigenvectors[μ_idx,a] * 
                                           scf_result.eigenvectors[μ_idx,i])
            end
        end
        ΔE = E_S[state] * 0.0367493  # eV to a.u.
        f[state] = (2/3) * ΔE * sum(abs2, μ)
    end
    
    println("Singlet excitations:")
    for i in 1:n_configs
        @printf("ε[%d] = %8.3f eV f = %8.3f \n", i, E_S[i], f[i])
    end

    println("Triplet excitations:")
    for i in 1:n_configs
        @printf("ε[%d] = %8.3f eV f = %8.3f \n", i, E_T[i], f[i])
    end

    return (
        CIResult(
            "Singlet CIS",
            E_S,
            C_S,
            configs,
            [argmax(abs.(C_S[:,i])) for i in 1:n_configs],
            f,
            ["$(c.from_orbitals[1])→$(c.to_orbitals[1])" for c in configs]
        ),
        CIResult(
            "Triplet CIS",
            E_T,
            C_T,
            configs,
            [argmax(abs.(C_T[:,i])) for i in 1:n_configs],
            zeros(n_configs),
            ["$(c.from_orbitals[1])→$(c.to_orbitals[1])" for c in configs]
        ))
end


"""
Calculate CISD matrix element ⟨Φᵢⱼᵃᵇ|H|Φₖₗᶜᵈ⟩ using PPP Hamiltonian
"""
# FIXME : unfinished, requires all matrix elements to be coded in for double-double and double-single
function calculate_cisd_matrix_element(config1::Configuration, config2::Configuration, 
                                   scf_result::SCFResult; singlet::Bool=true)
    if config1.type == SingleExcitation && config2.type == SingleExcitation
        i, a = config1.from_orbitals[1], config1.to_orbitals[1]
        j, b = config2.from_orbitals[1], config2.to_orbitals[1]
        if singlet==true    
            return calculate_cis_matrix_element(i, a, j, b, scf_result; singlet)
        else
            return calculate_cis_matrix_element(i, a, j, b, scf_result; triplet)
        end
    elseif config1.type == DoubleExcitation && config2.type == DoubleExcitation
        i, j = config1.from_orbitals
        a, b = config1.to_orbitals
        k, l = config2.from_orbitals
        c, d = config2.to_orbitals

        ΔE = scf_result.energies[a] + scf_result.energies[b] - scf_result.energies[i] - scf_result.energies[j]

        if singlet==true 
            if i == k && j == l && a == c && b == d
            # Diagonal element
                ΔE = scf_result.energies[a] + scf_result.energies[b] - scf_result.energies[i] - scf_result.energies[j]
                return ΔE + 2 * transform_two_electron_integral(i,a,j,b, scf_result) - transform_two_electron_integral(i,j,a,b, scf_result)
            end
            if i == k && j == l && a == c && b == d
            # Diagonal element
            return ΔE - transform_two_electron_integral(i,j,a,b, scf_result)
            end            
        end
    else
        return 0.0
    end
end

export run_cis_calculation
export SingleExcitation, DoubleExcitation, Configuration, CIResult

