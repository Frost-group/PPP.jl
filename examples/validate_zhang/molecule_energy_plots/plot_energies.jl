using Test
using PPP
using Printf
using StaticArrays
using CSV, DataFrames
using Gnuplot
include("plot_energy_fns.jl")

df_orca_st = CSV.read("examples/validate_zhang/molecule_energy_plots/singlet_triplet_ORCA.dat", DataFrame; delim=' ', ignorerepeated=true, header=true)

df_zhang_st = CSV.read("examples/validate_zhang/molecule_energy_plots/singlet_triplet_zhang.dat", DataFrame; delim=' ', ignorerepeated=true, header=true)

orca_st_data = construct_orca_st_data(df_orca_st)
zhang_st_data = construct_zhang_st_data(df_zhang_st)

cis_st_data = generate_ST_with_CI(df_orca_st, PPP.Zhang2011Model())