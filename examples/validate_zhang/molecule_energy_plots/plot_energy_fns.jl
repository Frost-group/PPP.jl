using Test
using PPP
using Printf
using StaticArrays
using CSV, DataFrames
using Gnuplot

#Compare Singlet + Triplet excitation energies from 3 different methods:
# 1) Zhang results from 2011 paper
# 2) ORCA calculations from PPPTraining
# 3) PPP SCF (from Zhang model) with CIS calculations

struct MoleculeSTData
    name::String
    singlet_E1::Union{Float64, Missing, Nothing}
    triplet_E1::Union{Float64, Missing, Nothing}
    ST_gap::Union{Float64, Missing, Nothing}
end

"""
    construct_orca_st_data
    Given the ORCA DataFrame, constructs a vector of MoleculeSTData structs with the molecule name, singlet first excitation energy, triplet first excitation energy, and ST gap.

"""

function construct_orca_st_data(df::DataFrame)::Vector{MoleculeSTData}
    n_molecules = Int(size(df)[1] / 2)
    orca_st_data = MoleculeSTData[]
    if isa(df[!,3], Vector{Float64})
        energy_column = df[!,3]
    else
        energy_column = [df[!,3][i] == "missing" ? missing : parse(Float64, df[!,3][i]) for i in 1:length(df[!,3])]
    end
    for mol_idx in 1:n_molecules
        name = df[!,1][2*mol_idx]
        #singlet_E1 = df[!,3][2*mol_idx-1]
        #triplet_E1 = df[!,3][2*mol_idx]
        singlet_E1 = energy_column[2*mol_idx-1]
        triplet_E1 = energy_column[2*mol_idx]
        push!(orca_st_data, MoleculeSTData(name, singlet_E1, triplet_E1, (singlet_E1 - triplet_E1)))
        # push!(orca_st_data, MoleculeSTData(name, singlet_E1, triplet_E1, (isa(singlet_E1, Float64) && isa(triplet_E1, Float64) ? (singlet_E1 - triplet_E1) : missing)))
    end
    return orca_st_data
end

"""
    construct_zhang_st_data
    Works for my manually constructed DataFrame of Zhang data, containing (where available) first excited singlet energy, first excited triplet energy, and ST gaps. Constructs a vector of MoleculeSTData structs with the molecule names and the above information.

    Note: some entries may be "missing" strings, which are converted to missing values in the output structs. In some cases only the ST gap is available.
    Note: the data has been collected as it is labelled in the Zhang 2011 paper. I suspect some of the ST gaps are mislabelled and should be called triplet energy rather than ST gaps.
"""

function construct_zhang_st_data(df::DataFrame)::Vector{MoleculeSTData}
    n_molecules = Int(size(df)[1] / 3)
    zhang_st_data = MoleculeSTData[]
    #E1 = copy(df[!,3])
    energy_column = [df[!,3][i] == "missing" ? missing : parse(Float64, df[!,3][i]) for i in 1:length(df[!,3])]
    #df[!,3] = replace.(df[!,3], "missing" => missing)
    #df[!,3] = passmissing(parse).(Float64, df[!,3])
    for mol_idx in 1:n_molecules
        name = df[!,1][3*mol_idx]
        singlet_E1 = energy_column[3*mol_idx-2]
        triplet_E1 = energy_column[3*mol_idx-1]
        ST_gap = energy_column[3*mol_idx]
        push!(zhang_st_data, MoleculeSTData(name, singlet_E1, triplet_E1, ST_gap))
    end
    return zhang_st_data
end

function generate_ST_with_CI(df_orca::DataFrame, model::PPP.AbstractModel)::Vector{MoleculeSTData}
    CI_mol_st_data = MoleculeSTData[]
    n_molecules = Int(size(df_orca)[1] / 2)
    for mol_idx in 1:n_molecules
        molecule_name = df_orca[!,1][2*mol_idx]
        whole_file_path = "examples/validate_zhang/molecule_energy_plots/molecule_ORCA/" * molecule_name * "_pbe1pbe-cc-pvtz.xyz"
        #println("$whole_file_path")
        try
            system, _, scf_result = PPP.run_ppp_calculation(whole_file_path, model);
            singlet, triplet = run_cis_calculation(system, scf_result);
            ST_gap = singlet.energies[1] - triplet.energies[1]
            push!(CI_mol_st_data, MoleculeSTData(molecule_name, singlet.energies[1], triplet.energies[1], ST_gap))
        catch e
            println("Error processing molecule: $molecule_name. Skipping...")
            push!(CI_mol_st_data, MoleculeSTData(molecule_name, nothing, nothing, nothing))
        end
    end
    return CI_mol_st_data
end

function plot_ST_comparison(orca_data::Vector{MoleculeSTData}, ci_data::Vector{MoleculeSTData}, zhang_data::Vector{MoleculeSTData})
    n_molecules = length(orca_data)
    molecule_names = ["$(orca_data[i].name)" for i in 1:n_molecules]
    #convert missing or nothing to NaN for plotting
    orca_ST_data = [(isa(orca_data[i].ST_gap, Float64) ? orca_data[i].ST_gap : NaN) for i in 1:n_molecules]
    # orca_ST_data = pushfirst!(orca_ST_data, NaN)
    # orca_ST_data = push!(orca_ST_data, NaN) #put NaN at start and end so the plot points aren't at the edges of the plot 
    ci_ST_data = [(isa(ci_data[i].ST_gap, Float64) ? ci_data[i].ST_gap : NaN) for i in 1:n_molecules]
    zhang_ST_data = [(isa(zhang_data[i].ST_gap, Float64) ? zhang_data[i].ST_gap : NaN) for i in 1:n_molecules]

    test_vcat = vcat(orca_ST_data, ci_ST_data, zhang_ST_data)
    max_value = maximum(filter(!isnan, test_vcat))

    x = 1:n_molecules

    xtic_labels = "set xtics (" *
    join(["\"$(molecule_names[i])\" " * string(i) for i in 1:n_molecules], ", ") *
    ") rotate by 90 right"

    @gp(
        "set term pdfcairo enhanced size 6in,5in",
        "set output '../zhang_plots/ST_gap_comparison_all3_updated.pdf'",
        "set xrange [0:" * string(n_molecules + 1) * "]",
        "set yrange [0:" * string(max_value + max_value*0.1) * "]",
        #"set xtics (" * join(["'$(molecule_names[i])' " * string(i) for i in 1:n_molecules], ", ") * ") rotate by 90 right",
        "set grid",
        xtic_labels,
        "set xtics rotate by 45 right",
        "set ytics 0, 0.5, 5.0",
        "set ylabel 'Singlet-Triplet Gap (eV)'",
        "set title 'Singlet-Triplet Gap Comparison: ORCA vs PPP.jl vs Zhang et al.'",

        "set linestyle 1 lc rgb '#2b8cbe' ps 1.2",  #blue
        "set linestyle 2 lc rgb '#31a354' ps 1.2",  #green
        "set linestyle 3 lc rgb '#c51b8a' ps 1.2",  #purple


        x, orca_ST_data, "w p ls 1 title 'ORCA'",
        x, ci_ST_data, "w p ls 2 title 'PPP-SCF + CIS'",
        x, zhang_ST_data, "w p ls 3 title 'Zhang 2011 results'",
    )
    @gp "set output"
end

function plot_energies_comparison(orca_data::Vector{MoleculeSTData}, ci_data::Vector{MoleculeSTData}, zhang_data::Vector{MoleculeSTData})
    n_molecules = length(orca_data)
    molecule_names = ["$(orca_data[i].name)" for i in 1:n_molecules]
    #convert missing or nothing to NaN for plotting
    orca_singlet_data = [(isa(orca_data[i].singlet_E1, Float64) ? orca_data[i].singlet_E1 : NaN) for i in 1:n_molecules]
    # println(orca_singlet_data)
    orca_triplet_data = [(isa(orca_data[i].triplet_E1, Float64) ? orca_data[i].triplet_E1 : NaN) for i in 1:n_molecules]
    # println(orca_triplet_data)

    ci_singlet_data = [(isa(ci_data[i].singlet_E1, Float64) ? ci_data[i].singlet_E1 : NaN) for i in 1:n_molecules]
    # println(ci_singlet_data)
    ci_triplet_data = [(isa(ci_data[i].triplet_E1, Float64) ? ci_data[i].triplet_E1 : NaN) for i in 1:n_molecules]
    # println(ci_triplet_data)

    zhang_singlet_data = [(isa(zhang_data[i].singlet_E1, Float64) ? zhang_data[i].singlet_E1 : NaN) for i in 1:n_molecules]
    # println(zhang_singlet_data)
    zhang_triplet_data = [(isa(zhang_data[i].triplet_E1, Float64) ? zhang_data[i].triplet_E1 : NaN) for i in 1:n_molecules]
    # println(zhang_triplet_data)
 
    zhang_ST_data = [(isa(zhang_data[i].ST_gap, Float64) ? zhang_data[i].ST_gap : NaN) for i in 1:n_molecules]
    # println(zhang_ST_data)

    x = 1:n_molecules
    test_vcat = vcat(orca_singlet_data, orca_triplet_data, ci_singlet_data, ci_triplet_data, zhang_singlet_data, zhang_triplet_data)
    max_value = maximum(filter(!isnan, test_vcat))
    # println("Max value: $max_value")
    
    xtic_labels = "set xtics (" *
    join(["\"$(molecule_names[i])\" " * string(i) for i in 1:n_molecules], ", ") *
    ") rotate by 90 right"

    @gp(
        "set term pdfcairo enhanced size 6in,5in",
        "set output '../zhang_plots/energy_comparison_all3_updated.pdf'",
        "set xrange [0:" * string(n_molecules + 1) * "]",
        "set yrange [0:" * string(max_value + max_value*0.1) * "]",
        "set grid",
        xtic_labels,
        "set xtics rotate by 45 right",
        "set ytics 0, 0.5, 5.0",
        "set ylabel 'Excitation Energy (eV)'",
        "set title 'Singlet and Triplet Excitation Energy: ORCA vs PPP.jl vs Zhang et al.'",

        "set linestyle 1 lc rgb '#2b8cbe' ps 1.2",  #blue
        "set linestyle 2 lc rgb '#31a354' ps 1.2",  #green
        "set linestyle 3 lc rgb '#c51b8a' ps 1.2",  #purple

        "set label 'ORCA' at screen 0.24,0.874 font ',15' left tc ls 1",
        "set label 'PPP-SCF + CIS' at screen 0.24,0.82225 font ',15' left tc ls 2",
        "set label 'Zhang 2011' at screen 0.24,0.8515 font ',15' left tc ls 3",

        "set key font \",15\"",
        x, orca_singlet_data, "w p pt 2 ps 1.2 lc rgb '#2b8cbe' title 'singlet'",  #blue
        x, orca_triplet_data, "w p pt 1 ps 1.2 lc rgb '#2b8cbe' title 'triplet'",
        x, ci_singlet_data, "w p pt 2 ps 1.2 lc rgb '#31a354' notitle", #green
        x, ci_triplet_data, "w p pt 1 ps 1.2 lc rgb '#31a354' notitle",
        x, zhang_singlet_data, "w p pt 2 ps 1.2 lc rgb '#c51b8a' notitle", #purple
        x, zhang_triplet_data, "w p pt 1 ps 1.2 lc rgb '#c51b8a' notitle",
        #x, zhang_ST_data, "w p pt 8 lc rgb '#d62728' title 'Zhang ST gap'",
    )
    @gp "set output"
end

df_orca_st = CSV.read("examples/validate_zhang/molecule_energy_plots/singlet_triplet_ORCA.dat", DataFrame; delim=' ', ignorerepeated=true, header=true)

df_zhang_st = CSV.read("examples/validate_zhang/molecule_energy_plots/singlet_triplet_zhang.dat", DataFrame; delim=' ', ignorerepeated=true, header=true)

orca_st_data = construct_orca_st_data(df_orca_st)
zhang_st_data = construct_zhang_st_data(df_zhang_st)
cis_st_data = generate_ST_with_CI(df_orca_st, PPP.Zhang2011Model())
plot_ST_comparison(orca_st_data, cis_st_data, zhang_st_data)
plot_energies_comparison(orca_st_data, cis_st_data, zhang_st_data)

