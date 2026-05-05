# Optim OrgoneModel parameters to fit rattled Cyclobutadiene
using PPP
using Printf
using Optim
using Logging

# Silent mode (default - only warnings/errors)
global_logger(ConsoleLogger(stderr, Logging.Warn))

# Reference data: (mol_path, S1_ref, T1_ref) from ORCA EOM-CCSD(T)
const DATASET = [
    ("examples/optim_zhang_rattled/01-rattled01-cyclobutadiene.xyz", 3.234, 1.427),
    ("examples/optim_zhang_rattled/01-rattled02-cyclobutadiene.xyz", 3.145, 1.305),
    ("examples/optim_zhang_rattled/01-rattled03-cyclobutadiene.xyz", 3.055, 1.242),
    ("examples/optim_zhang_rattled/01-rattled04-cyclobutadiene.xyz", 2.831, 0.978),
    ("examples/optim_zhang_rattled/01-rattled05-cyclobutadiene.xyz", 3.149, 1.318),
    ("examples/optim_zhang_rattled/01-rattled06-cyclobutadiene.xyz", 3.026, 1.206),
    ("examples/optim_zhang_rattled/01-rattled07-cyclobutadiene.xyz", 2.904, 1.060),
    ("examples/optim_zhang_rattled/01-rattled08-cyclobutadiene.xyz", 2.953, 1.127),
    ("examples/optim_zhang_rattled/01-rattled09-cyclobutadiene.xyz", 2.780, 0.978),
    ("examples/optim_zhang_rattled/01-rattled10-cyclobutadiene.xyz", 3.168, 1.352),
]

# We optimize 10 parameters defining Carbon in the Orgone model
const PARAM_NAMES = ["ζ_core", "ζ_diff", "ζ_wiggle", "w_core", "w_diff", "w_wiggle", "χ", "α", "Z_eff", "k_t"]

function make_orgone(x)
    # Allow negative coefficient for Löwdin hole; enforce sum = 1.0
    c_norm = [x[4], x[5], 1.0 - x[4] - x[5]]

    # Inherit physical identity params from defaults
    def_C = default_orgone_model(3).params[:C]

    c_params = PPP.OrgoneAtomParams(
        ζ=[x[1], x[2], x[3]],
        c=c_norm,
        χ=x[7],
        α=x[8],
        Z_eff=x[9],
        n_star=def_C.n_star,
        k_t=x[10],
        IP=def_C.IP,
        pi_electrons=def_C.pi_electrons
    )

    # Start with default and overwrite Carbon
    model = default_orgone_model(3)
    model.params[:C] = c_params
    return model
end

# Initial guess: mapped directly to c_k (x[6] is unused now, but kept for array length)
x0 = [4.5, 0.85, 1.5, 0.60, 0.55, 0.32, 2.5, 0.05, 3.25, 1.0]

# Bounds (allow negative c_k for hole)
lb = [1.0, 0.1, 0.5, -2.0, -2.0, 0.01, 0.0, 0.0, 1.5, 0.1]
ub = [10.0, 2.0, 5.0, 2.0, 2.0, 2.0, 5.0, 0.5, 5.0, 5.0]

println("="^60)
println("PPP Orgone Parameter Optimization (Rattled Cyclobutadiene)")
println("="^60)
println("\nInitial parameters:")
for (name, val) in zip(PARAM_NAMES, x0)
    @printf("  %12s = %10.5f\n", name, val)
end

println("\nOptim over $(length(DATASET)) molecules...")
println("-"^60)

# SAMIN handles the bounded optimization
result = optim_ppp(make_orgone, x0, DATASET; lb=lb, ub=ub, options=Optim.Options(show_trace=true, iterations=2000))

println("\nOptimization Summary:")
println("  Converged ?:      ", Optim.converged(result))
println("  Iterations:     ", Optim.iterations(result))
println("  Final loss:     ", Optim.minimum(result))

println("-"^60)
println("\nFinal parameters:")
x_opt = Optim.minimizer(result)
for (i, name) in enumerate(PARAM_NAMES)
    @printf("  %12s = %10.5f (was %10.5f)\n", name, x_opt[i], x0[i])
end

println("\n" * "="^60)
println("Comparison: Reference data vs Optim Orgone vs Default Orgone")
println("="^60)
model_opt = make_orgone(x_opt)
model_def = default_orgone_model(3)

@printf("%25s     States   %8s  %8s  %8s\n", "Molecule", "Reference", "Opt Org", "Def Org")
println("-"^60)
for (mol_path, S1_ref, T1_ref) in DATASET
    sys, _, scf = run_ppp_calculation(mol_path, model_opt)
    singlet, triplet = run_cis_calculation(sys, scf)
    S1, T1 = singlet.energies[1], triplet.energies[1]

    sys, _, scf = run_ppp_calculation(mol_path, model_def)
    singlet, triplet = run_cis_calculation(sys, scf)
    S1_def, T1_def = singlet.energies[1], triplet.energies[1]

    mol_name = basename(mol_path)[1:end-4]
    @printf("%25s   Singlets %8.3f  %8.3f  %8.3f\n", mol_name, S1_ref, S1, S1_def)
    @printf("%25s   Triplets %8.3f  %8.3f  %8.3f\n", mol_name, T1_ref, T1, T1_def)
end
