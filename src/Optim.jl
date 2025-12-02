# Optimise PPP models using Optim.jl
# Hacked together, only works for Zhang currently. 

using Optim

const MoleculeReference = Tuple{String, Float64, Float64}  # (mol_path, S1_ref, T1_ref)

"""
    build_loss(make_model::Function, dataset::Vector{MoleculeReference}) -> Function

This hurts my head. Some deep Julia magic; why not just generate a function with a function? 
What Would Church Do (WWCD)

Arguments:
- `make_model`: Function mapping parameter vector x::Vector{Float64} to an AbstractModel
- `dataset`: Vector of (molecule_path, S1_ref, T1_ref) tuples

Returns a loss function `loss(x)` that computes sum of squared errors.
"""
function build_loss(make_model::Function, dataset::Vector{MoleculeReference})
    function loss(x::Vector{Float64})
        model = make_model(x)
        total = 0.0
        for (mol_path, S1_ref, T1_ref) in dataset
            sys, _, scf = run_ppp_calculation(mol_path, model)
            singlet, triplet = run_cis_calculation(sys, scf)
            S1, T1 = singlet.energies[1], triplet.energies[1]
            total += (S1 - S1_ref)^2 + (T1 - T1_ref)^2
        end
        return total
    end
end

"""
    optimize_ppp(make_model, x0, dataset; lb, ub, options) -> Optim.OptimizationResults

Optimize PPP model parameters to fit reference S1/T1 data.

Arguments:
- `make_model`: Function mapping parameter vector to AbstractModel
- `x0`: Initial parameter guess
- `dataset`: Vector of (molecule_path, S1_ref, T1_ref) tuples
- `lb`: Lower bounds (default: -Inf)
- `ub`: Upper bounds (default: Inf)
- `options`: Optim.Options (default: show_trace=true)

Returns Optim.jl result object. Access optimal params via `Optim.minimizer(result)`.
"""
function optim_ppp(make_model::Function, x0::Vector{Float64}, dataset::Vector{MoleculeReference};
                      lb::Vector{Float64}=fill(-Inf, length(x0)),
                      ub::Vector{Float64}=fill(Inf, length(x0)),
                      options::Optim.Options=Optim.Options(show_trace=true, iterations=100))
    loss = build_loss(make_model, dataset)
#    optimize(loss, lb, ub, x0, Fminbox(NelderMead()), options)
    optimize(loss, lb, ub, x0, SAMIN(), options)


# Currently gradient free; god knows whether autodiff will work. Probably not. 
end

export MoleculeReference, build_loss, optim_ppp

