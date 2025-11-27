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
    singles = generate_single_excitations(n_orbs, n_occ)
    doubles = generate_double_excitations(n_orbs, n_occ)
    return vcat(singles,doubles)
end

# Transform PPP two-electron integrals from AO to MO basis: (ij|ab)
function transform_two_electron_integral(i::Int, j::Int, a::Int, b::Int, scf_result::SCFResult)
    C = scf_result.eigenvectors
    V = scf_result.V_raw
    return sum(C[μ,i] * C[μ,j] * C[ν,a] * C[ν,b] * V[μ,ν]
              for μ in axes(C,1), ν in axes(C,1))
end

# Transform Fock Matrix from AO to MO basis for 1 electron integrals
function transform_fock(scf_result::SCFResult)
    F_AO = scf_result.F
    C = scf_result.eigenvectors
    return C'*F_AO*C
end

"""
Calculate CIS matrix element ⟨Φᵢᵃ|H|Φⱼᵇ⟩ using PPP Hamiltonian

Singlet matrix elements from Table 4.1 in Szabo & Ostlund: Modern Quantum Chemistry (Singlets) 
Triplet matrix elements from: 
    Gamba, A., Tantardini, G.F., Simonetta, M., 1972. 
    A study of ground and excited states of biphenyl by the “Molecules in Molecules” method. Spectrochimica Acta Part A: Molecular Spectroscopy 28, 1877–1888. 
    https://doi.org/10.1016/0584-8539(72)80159-4

"""
function calculate_singlet_cis_matrix_element(i::Int, a::Int, j::Int, b::Int, scf_result::SCFResult)
    ε = @view scf_result.energies[:]
    # Singlet matrix elements: (ε_a-ε_i)*δ_ab*δ_ij + 2*(ai|jb) - (ab|ji)
    return (ε[a]-ε[i])*((i==j) && (a==b)) + 2*transform_two_electron_integral(a,i,j,b, scf_result) - transform_two_electron_integral(a,b,j,i, scf_result)
end 

function calculate_triplet_cis_matrix_element(i::Int, a::Int, j::Int, b::Int, scf_result::SCFResult)
    ε = @view scf_result.energies[:]    
    # Triplet matrix elements: (ε_a-ε_i)*δ_ab*δ_ij - (ab|ji)
    return (ε[a]-ε[i])*((i==j) && (a==b)) - transform_two_electron_integral(a,b,j,i, scf_result)
end

function run_cis_calculation(system::MolecularSystem, scf_result::SCFResult)
    n_orbs = size(scf_result.eigenvectors, 2)
    n_occ = system.n_electrons ÷ 2
    
    @debug "Starting CIS calculation" n_occ n_orbs
    @debug "Orbital energies (eV)" energies=scf_result.energies
    
    # Generate all single excitations
    configs = generate_single_excitations(n_orbs, n_occ)
    n_configs = length(configs)
    @debug "CIS configurations" n_configs
    
    # Build CIS matrix
    H_S = zeros(n_configs, n_configs)
    H_T = zeros(n_configs, n_configs)

    for p in 1:n_configs, q in 1:p
        i, a = configs[p].from_orbitals[1], configs[p].to_orbitals[1]
        j, b = configs[q].from_orbitals[1], configs[q].to_orbitals[1]
        H_S[p,q] = calculate_singlet_cis_matrix_element(i,a,j,b,scf_result)
        H_S[q,p] = H_S[p,q]
        H_T[p,q] = calculate_triplet_cis_matrix_element(i,a,j,b,scf_result)
        H_T[q,p] = H_T[p,q]
    end
    
    @debug "CIS singlet matrix (10x10 block)" singlet=H_S[1:min(10,n_configs), 1:min(10,n_configs)] 
    @debug "CIS triplet matrix (10x10 block)" triplet=H_T[1:min(10,n_configs), 1:min(10,n_configs)]
    
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
        ΔE = E_S[state] * eV2Ha  # eV to a.u.
        f[state] = (2/3) * ΔE * sum(abs2, μ)
    end
    
    @debug "Singlet excitations (eV) oscillator strengths" energies=E_S oscillator_strengths=f
    @debug "Triplet excitations (eV)" energies=E_T

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
Calculate dynamic spin polarisation correction using perturbation theory

Perturbation theory on CIS states according to: https://doi.org/10.1007/BF00549021

"""
function calculate_dsp(system::MolecularSystem, scf_result::SCFResult, CIS_S::CIResult, CIS_T::CIResult)
    n_occ = system.n_electrons ÷ 2
    n_orbs = size(scf_result.eigenvectors, 2)
    homo_idx = n_occ
    lumo_idx = n_occ + 1

    F_MO = transform_fock(scf_result)
    # single excitations excluding HOMO and LUMO
    single_excitations = [(i,a) for i in 1:n_occ-1, a in n_occ+2:n_orbs]
    
    excitations = Dict()
    for (i,a) in single_excitations
        k_x = transform_two_electron_integral(i,homo_idx,homo_idx,a, scf_result)
        k_y = transform_two_electron_integral(i,lumo_idx,lumo_idx,a, scf_result)

        gap_s = F_MO[a,a] + F_MO[lumo_idx,lumo_idx] - F_MO[homo_idx,homo_idx] - F_MO[i,i] - (CIS_S.energies[1])
        gap_t = F_MO[a,a] + F_MO[lumo_idx,lumo_idx] - F_MO[homo_idx,homo_idx] - F_MO[i,i] - (CIS_T.energies[1])

        s_1 = (3/2) * (k_x - k_y)^2 / gap_s
        t_1 = (1/2) * (k_x - k_y)^2 / gap_t
        t_2 = (k_x + k_y)^2 / gap_t

        excitations[(i,a)] = Dict("s_1" => s_1,
            "t_1" => t_1,
            "t_2" => t_2,
            "dsp" => -s_1 + t_1 + t_2)
    end

    dsp = sum(excitation["dsp"] for excitation in values(excitations))
    return dsp
end

"""
Calculate CISD matrix element ⟨Φᵢⱼᵃᵇ|H|Φₖₗᶜᵈ⟩ using PPP Hamiltonian
 Attempted to use matrix elements from: 'A study of ground and excited states of biphenyl by the “Molecules in Molecules” method' https://doi.org/10.1016/0584-8539(72)80159-4
 Notation in this paper doesn't match Szabo and Ostlund. Szabo and Ostlund does not feature matrix elements for <S|H|D> (or vice versa).
"""
# FIXME : not working. Tried implementing using Slater-Condon rules
# <D|H|D> matrix elements for CISD  
function calculate_double_double_matrix_element(p::Configuration, q::Configuration, scf_result::SCFResult)
    i, j = p.from_orbitals
    a, b = p.to_orbitals
    k, l = q.from_orbitals
    c, d = q.to_orbitals
    F = transform_fock(scf_result)
    # Identical determinants
    if i==k && j==l && a==c && b==d 
        return F[a,a] + F[b,b] - F[i,i] - F[j,j] + 2*transform_two_electron_integral(a,i,j,b, scf_result) - transform_two_electron_integral(a,b,j,i, scf_result)
    # Determinants differ by 1 spin orbital
    elseif i!==k && j==l && a==c && b==d
        return -F[i,k] + transform_two_electron_integral(a,i,k,b,scf_result)
    elseif i==k && j!==l && a==c && b==d
        return -F[j,l] + transform_two_electron_integral(a,j,l,b,scf_result)
    elseif i==k && j==l && a!==c && b==d
        return F[a,c] + transform_two_electron_integral(a,i,j,c,scf_result)
    elseif i==k && j==l && a==c && b!==d
        return F[b,d] + transform_two_electron_integral(b,i,j,d,scf_result)
    # Determinants differ by 2 spin orbitals
    elseif i==k && j!==l && a!==c && b==d
        return transform_two_electron_integral(a,j,l,c,scf_result)
    elseif i!==k && j==l && a!==c && b==d
        return transform_two_electron_integral(a,i,k,c,scf_result)
    elseif i!==k && j!==l && a==c && b==d
        return transform_two_electron_integral(i,j,l,k,scf_result)
    elseif i!==k && j==l && a==c && b!==d
        return transform_two_electron_integral(b,i,k,d,scf_result)
    elseif i==k && j!==l && a==c && b!==d
        return transform_two_electron_integral(b,j,l,d,scf_result)
    elseif i==k && j==l && a!==c && b!==d
        return transform_two_electron_integral(a,b,d,c,scf_result)
    # Determinants differ by 3 or more spin orbitals
    else
        return 0.0
    end
end

# <S|H|D> || <D|H|S> matrix elements for CISD
function calculate_single_double_matrix_element(p::Configuration,q::Configuration, scf_result::SCFResult)
    i, a = p.from_orbitals[1], p.to_orbitals[1] #single case
    j, k = q.from_orbitals # double case
    b, c = q.to_orbitals # double case
    F = transform_fock(scf_result)

    # Determinants always differ by 2 spin orbitals
    if i==j && a==b
        return 2*transform_two_electron_integral(a,i,k,c, scf_result) - transform_two_electron_integral(a,c,k,i, scf_result)
    elseif i==j && a==c
        return 2*transform_two_electron_integral(a,i,k,b, scf_result) - transform_two_electron_integral(a,b,k,i, scf_result)
    elseif i==k && a==b
        return 2*transform_two_electron_integral(a,i,j,c, scf_result) - transform_two_electron_integral(a,c,j,i, scf_result)
    elseif i==k && a==c 
        return 2*transform_two_electron_integral(a,i,j,b, scf_result) - transform_two_electron_integral(a,b,j,i, scf_result)
    else
        return 0.0
    end
end

# Builds CISD matrix from Single-Single, Single-Double, Double-Single and Double-Double
function calculate_singlet_cisd_matrix_element(p::Configuration,q::Configuration, scf_result::SCFResult)
    if p.type isa SingleExcitation && q.type isa SingleExcitation
        i, a = p.from_orbitals[1], p.to_orbitals[1]
        j, b = q.from_orbitals[1], q.to_orbitals[1]
        return calculate_singlet_cis_matrix_element(i,a,j,b,scf_result)
    elseif p.type isa DoubleExcitation && q.type isa DoubleExcitation
        return calculate_double_double_matrix_element(p,q,scf_result)
    elseif p.type isa SingleExcitation && q.type isa DoubleExcitation
        return calculate_single_double_matrix_element(p,q,scf_result)
    elseif p.type isa DoubleExcitation && q.type isa SingleExcitation
        return calculate_single_double_matrix_element(q,p,scf_result)
    else
        return 0.0
    end
end

# Nearly identical to run CIS calculation
function run_cisd_calculation(system::MolecularSystem, scf_result::SCFResult)
    n_orbs = size(scf_result.eigenvectors, 2)
    n_occ = system.n_electrons ÷ 2
    
    @debug "Starting CISD calculation" n_occ n_orbs
    @debug "Orbital energies (eV)" energies=scf_result.energies
    
    # Generate all excitations
    configs = generate_excitations(n_orbs, n_occ)
    n_configs = length(configs)
    @debug "CISD configurations: " n_configs
    
    # Build CISD matrix
    H_CISD = zeros(n_configs, n_configs)

    for p in 1:n_configs, q in 1:p
        H_CISD[p,q] = calculate_singlet_cisd_matrix_element(configs[p], configs[q], scf_result)
        H_CISD[q,p] = H_CISD[p,q]
    end
    
    @debug "CISD matrix (10x10 block)" H_CISD=H_CISD[1:min(10,n_configs), 1:min(10,n_configs)]
    
    # Solve eigenvalue problem
    E_CISD, C_CISD = eigen(Symmetric(H_CISD))

    # Calculate oscillator strengths (Copied from CIS)
    f = zeros(n_configs)
    for state in 1:n_configs
        μ = zeros(2)  # transition dipole in x,y
        for (idx, config) in enumerate(configs)
            # Sum over all configurations contributing to this state
            i, a = config.from_orbitals[1], config.to_orbitals[1]
            for μ_idx in 1:length(system.atoms)
                pos = system.atoms[μ_idx].position[1:2]
                μ .+= C_CISD[idx,state] * pos * (scf_result.eigenvectors[μ_idx,a] * 
                                           scf_result.eigenvectors[μ_idx,i])
            end
        end
        ΔE = E_CISD[state] * eV2Ha  # eV to a.u.
        f[state] = (2/3) * ΔE * sum(abs2, μ)
    end

    chars = String[]
    for config in configs
        if config.type isa SingleExcitation
            push!(chars, "$(config.from_orbitals[1])→$(config.to_orbitals[1])")
        elseif config.type isa DoubleExcitation
            push!(chars, "($(config.from_orbitals[1]),$(config.from_orbitals[2]))→($(config.to_orbitals[1]),$(config.to_orbitals[2]))")
        end
    end

    return (CIResult(
        "CISD",
        E_CISD,
        C_CISD,
        configs,
        [argmax(abs.(C_CISD[:,i])) for i in 1:n_configs],
        f,
        chars
    ))
end

export run_cis_calculation, run_cisd_calculation, calculate_dsp
export SingleExcitation, DoubleExcitation, Configuration, CIResult