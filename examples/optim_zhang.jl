# Optim Zhang2011Model parameters to fit initial Orca EOM ref data

using PPP
using Printf
using Optim

using Logging
# Silent mode (default - only warnings/errors - otherwise gets loads of PPP chatter)
global_logger(ConsoleLogger(stderr, Logging.Warn))

# Reference data: (mol_path, S1_ref, T1_ref) from ORCA EOM-CCSD(T) (actually, was it pertubative triplets?)
const DATASET = [
 #   ("molecules/cyclobutadiene.xyz",    3.159, 1.367),
    ("molecules/butadiene.xyz",         6.507, 3.232),
    ("molecules/hexatriene.xyz",        5.442, 2.475),
    ("molecules/pentalene.xyz",         1.909, 0.993),
#    ("molecules/cyclooctatetraene.xyz", 4.123, 2.797),
    ("molecules/azulene.xyz",           1.833, 1.923),
#    ("molecules/octatetraene.xyz",      4.701, 2.045),
#    ("molecules/heptalene.xyz",         2.460, 1.606),
#    ("molecules/phenanthrene.xyz",      3.711, 2.641),
#    ("molecules/decapentaene.xyz",      4.171, 1.759),
#    ("molecules/dodecahexaene.xyz",     3.766, 1.545),
#    ("molecules/tetradecaheptaene.xyz", 3.453, 1.406),
]

# Map parameter vector to Zhang2011Model (all 9 parameters)
const PARAM_NAMES = ["t_ii_E", "t_ii_I", "t_ij_E", "t_ij_B", 
                     "γ_ii_E", "γ_ii_I", "γ_ii_I_B", "γ_ij_E", "γ_ij_I"]

function make_zhang(x)
    Zhang2011Model(
        t_ii_E  = x[1],
        t_ii_I   = x[2],

        t_ij_E     = x[3],
        t_ij_B     = x[4],

        γ_ii_E  = x[5],
        γ_ii_I = x[6],
        γ_ii_I_B = x[7],
        
        γ_ij_E     = x[8],
        γ_ij_I     = x[9],
    )
end

# Initial guess (original Zhang2011 values; eV and Å or other factors therefore)
x0 = [-9.13597, 0.15076, -28.07749, 1.65878, 11.56718, 0.16640, 2.30448, 21.88221, 7.05909]
# Bounds: guess physically reasonable ranges
lb = [-15.0,    0.01,    -50.0,     0.5,      0.0,     0.01,    1.5,     15.0,     4.0]
ub = [-5.0,     1.0,     -10.0,     3.0,     15.0,     1.0,     3.5,     30.0,    12.0]

println("="^60)
println("PPP Zhang2011 Parameter Optimization")
println("="^60)
println("\nInitial parameters:")
for (name, val) in zip(PARAM_NAMES, x0)
    @printf("  %12s = %10.5f\n", name, val)
end

println("\nOptim over $(length(DATASET)) molecules...")
println("-"^60)

# HOLD MY BEER
result = optim_ppp(make_zhang, x0, DATASET; lb=lb, ub=ub)

# Report optimiser summary
println("\nOptimization Summary:")
println("  Converged ?:      ", Optim.converged(result))
println("  Iterations:     ", Optim.iterations(result))
println("  f(x) calls:     ", result.f_calls)
println("  Final loss:     ", Optim.minimum(result))

println("-"^60)
println("\nFinal parameters:")
x_opt = Optim.minimizer(result)
for (i, name) in enumerate(PARAM_NAMES)
    @printf("  %12s = %10.5f (was %10.5f)\n", name, x_opt[i], x0[i])
end

println("\nInitial loss: ", build_loss(make_zhang, DATASET)(x0))
println("\nFinal loss: ", Optim.minimum(result))
println("Converged (?): ", Optim.converged(result))

# Compare predictions with optimized model
println("\n" * "="^60)
println("Comparison: Reference data vs PPP-Zhang2011 vs Optim PPP-Zhang2011")
println("="^60)
model_opt = make_zhang(x_opt)
model_zhang = Zhang2011Model()
println("Molecule Name          | Reference | Optim PPP | Zhang2011 PPP")
println("-"^60)
for (mol_path, S1_ref, T1_ref) in DATASET
    sys, _, scf = run_ppp_calculation(mol_path, model_opt)
    singlet, triplet = run_cis_calculation(sys, scf)
    S1, T1 = singlet.energies[1], triplet.energies[1]

    sys, _, scf = run_ppp_calculation(mol_path, model_zhang)
    singlet, triplet = run_cis_calculation(sys, scf)
    S1_zhang, T1_zhang = singlet.energies[1], triplet.energies[1]

    mol_name = basename(mol_path)[1:end-4]
    @printf("%15s   Singlets %8.3f  %8.3f  %8.3f\n", mol_name, S1_ref, S1, S1_zhang)
    @printf("%15s   Triplets %8.3f  %8.3f  %8.3f\n", mol_name, T1_ref, T1, T1_zhang)
end

