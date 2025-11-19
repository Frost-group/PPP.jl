using Test
using PPP
using Printf
using StaticArrays
using CSV, DataFrames
using Gnuplot

# Gnuplot multiplot to compare excited state energies: CASPT2 vs Zhang closed form vs CIS

function remove_eqb_offset(energies::Vector{Float64})
    offset = energies[6] # assuming equilibrium bond length at 6th element in data frame column as per Zhang2011Ethene.dat
    return energies .- offset
end

df_zhang_CASPT2 = CSV.read("examples/validate_zhang/Zhang2011Ethene.dat", DataFrame;
              delim=' ',
              ignorerepeated=true,
              comment="#",
              header=false)

df_Energies_ethene = CSV.read("examples/validate_zhang/ValidateZhang_CIS_energy_results.csv", DataFrame; delim=',')

rs = df_zhang_CASPT2[:, 1] #Vector{Float64}
ES0s_Zhang = df_zhang_CASPT2[:, 2]
ES0s_Zhang .*= Ha2eV
ES1s_Zhang = df_zhang_CASPT2[:, 3] # 11B1u
ES1s_Zhang .*= Ha2eV
ET1s_Zhang = df_zhang_CASPT2[:, 4] # 13B1u
ET1s_Zhang .*= Ha2eV

caspt2_excitation_singlet = ES1s_Zhang .- ES0s_Zhang
caspt2_excitation_triplet = ET1s_Zhang .- ES0s_Zhang

# subplots with singlet and triplet

@gp "set multiplot layout 1,2"
# --- singlet plot ---
@gp :- 1 "set xlabel 'C=C bond length (Å)'"
@gp :- "set ylabel 'Excitation energy (eV)'"
@gp :- "set grid"
@gp :- "set yrange [1.5:9.0]"
@gp :- "set key top right spacing 1.5"
@gp :- "set title 'Singlet'"
@gp :- rs caspt2_excitation_singlet "w lp title 'Zhang CASPT2 1^{1}B_{1u} - 1^{1}A_{g}'"
@gp :- rs (df_Energies_ethene.ES1s_PPP .- df_Energies_ethene.ES0s_PPP) "w lp title 'Zhang E_{S_{1}} - E_{S_{0}}'"
@gp :- rs df_Energies_ethene.ESs_CIS "w lp title 'CIS'"

# --- triplet plot ---
@gp :- 2 "unset ylabel"
@gp :- "set format y ''"
@gp :- "set grid"
@gp :- "set yrange [1.5:9.0]"
@gp :- "set key top right spacing 1.5"
@gp :- "set title 'Triplet'"
@gp :- "set xlabel 'C=C bond length (Å)'"
@gp :- rs caspt2_excitation_triplet "w lp title 'Zhang CASPT2 1^{3}B_{1u} - 1^{1}A_{g}'"
@gp :- rs (df_Energies_ethene.ET1s_PPP .- df_Energies_ethene.ES0s_PPP) "w lp title 'Zhang E_{T_{1}} - E_{S_{0}}'"
@gp :- rs df_Energies_ethene.ETs_CIS "w lp title 'CIS'"
Gnuplot.save("../../Chem/zhang_plots/ethene_excited_states_subplots.pdf", term="pdfcairo enhanced size 8in,3in")
@gp :- "unset multiplot"




# @gp(
#     "set term pdfcairo enhanced size 8in,3in",
#     "set output '../../Chem/zhang_plots/ethene_excited_states.pdf'",
#     "set multiplot layout 1,2 rowsfirst",

#     # --- singlet plot ---
#     "set xlabel 'C=C bond length (Å)'",
#     "set ylabel 'Excitation energy (eV)'",
#     "set grid",
#     "set yrange [1.5:9.0]",
#     "set key top right spacing 1.5",
#     "set title 'Singlet'",
#     rs, caspt2_excitation_singlet, "w lp title 'Zhang CASPT2 1^{1}B_{1u} - 1^{1}A_{g}'",
#     rs, (df_Energies_ethene.ES1s_PPP - df_Energies_ethene.ES0s_PPP), "w lp title 'Zhang E_{S_{1}} - E_{S_{0}}'",
#     rs, df_Energies_ethene.ESs_CIS, "w lp title 'CIS'",

#     # --- triplet plot ---
#     "set xlabel 'C=C bond length (Å)'",
#     "unset ylabel",
#     "set format y ''",
#     "set grid",
#     "set yrange [1.5:9.0]",
#     "set key top right spacing 1.5",
#     "set title 'Triplet'",
#     rs, caspt2_excitation_triplet, "w lp title 'Zhang CASPT2 1^{3}B_{1u} - 1^{1}A_{g}'",
#     rs, (df_Energies_ethene.ET1s_PPP - df_Energies_ethene.ES0s_PPP), "w lp title 'Zhang E_{T_{1}} - E_{S_{0}}'",
#     rs, df_Energies_ethene.ETs_CIS, "w lp title 'CIS'",

#     "unset multiplot"
#     )
# @gp "set output"