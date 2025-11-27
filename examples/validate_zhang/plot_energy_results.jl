using Test
using PPP
using Printf
using StaticArrays
using CSV, DataFrames
using Gnuplot

#NOTE: NEED TO SET YOUR OWN DIRECTORY FOR OUTPUT PLOTS

function remove_eqb_offset(energies::Vector{Float64})
    offset = energies[6] # assuming equilibrium bond length at 6th element in data frame column as per Zhang2011Ethene.dat
    return energies .- offset
end

df_zhang_CASPT2 = CSV.read("examples/validate_zhang/Zhang2011Ethene.dat", DataFrame;
              delim=' ',
              ignorerepeated=true,
              comment="#",
              header=false)

df_Energies_ethene = CSV.read("examples/validate_zhang/ValidateZhang_CIS_energy_results_new.csv", DataFrame; delim=',')

rs = df_zhang_CASPT2[:, 1] #Vector{Float64}
ES0s_Zhang = df_zhang_CASPT2[:, 2]
ES0s_Zhang .*= Ha2eV
ES1s_Zhang = df_zhang_CASPT2[:, 3] # 11B1u
ES1s_Zhang .*= Ha2eV
ET1s_Zhang = df_zhang_CASPT2[:, 4] # 13B1u
ET1s_Zhang .*= Ha2eV

caspt2_excitation_singlet = ES1s_Zhang .- ES0s_Zhang
caspt2_excitation_triplet = ET1s_Zhang .- ES0s_Zhang

# first excited state (singlet)

@gp(
    "set term pdfcairo enhanced size 4in,3in",
    "set output '../../chem/zhang_plots/ethene_excited_singlet_NEW.pdf'",
    "set xlabel 'C=C bond length (Å)'",
    "set ylabel 'Excitation energy (eV)'",
    "set yrange [1.5:9.0]",
    "set grid",
    "set key top right spacing 1.5",
    rs, caspt2_excitation_singlet, "w lp title 'Zhang CASPT2 1^{1}B_{1u} - 1^{1}A_{g}'",
    rs, (df_Energies_ethene.ES1s_PPP - df_Energies_ethene.ES0s_PPP), "w lp title 'Zhang E_{S_{1}} - E_{S_{0}}'",
    rs, df_Energies_ethene.ESs_CIS, "w lp title 'CIS'"
    )
@gp "set output"

# triplet state
@gp(
    "set term pdfcairo enhanced size 4in,3in",
    "set output '../../chem/zhang_plots/ethene_excited_triplet_NEW.pdf'",
    "set xlabel 'C=C bond length (Å)'",
    "set ylabel 'Excitation energy (eV)'",
    "set grid",
    "set yrange [1.5:9.0]",
    "set key top right spacing 1.5",
    rs, caspt2_excitation_triplet, "w lp title 'Zhang CASPT2 1^{3}B_{1u} - 1^{1}A_{g}'",
    rs, (df_Energies_ethene.ET1s_PPP - df_Energies_ethene.ES0s_PPP), "w lp title 'Zhang E_{T_{1}} - E_{S_{0}}'",
    rs, df_Energies_ethene.ETs_CIS, "w lp title 'CIS'"
    )
@gp "set output"


# ground state comparison
@gp(
    "set term pdfcairo enhanced size 4in,3in",
    "set output '../../chem/zhang_plots/ethene_ground_state_NEW.pdf'",
    "set xlabel 'C=C bond length (Å)'",
    "set ylabel 'Energy (eV)'",
    "set grid",
    "set key top left spacing 1.5",
    rs, remove_eqb_offset(df_Energies_ethene.PPPs), "w lp title 'PPP SCF'",
    rs, remove_eqb_offset(df_Energies_ethene.Huckels), "w lp title 'Huckel'",
    rs, remove_eqb_offset(df_Energies_ethene.ES0s_PPP), "w lp title 'Zhang E_{S_{0}}'",
    rs, remove_eqb_offset(ES0s_Zhang), "w lp title 'Zhang CASPT2 1^{1}A_{g}'"
    )
@gp "set output"