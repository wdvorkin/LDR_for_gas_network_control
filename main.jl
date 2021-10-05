#!/usr/bin/env julia
using CSV, DataFrames, PowerModels, Tables, ArgParse
using JuMP, Ipopt, Mosek, MosekTools
using Statistics, Distributions, Measures, StatsBase
using SparseArrays, LinearAlgebra
using CDDLib, Polyhedra

# cd to directory
cd(dirname(@__FILE__))
PowerModels.silence()

# load functions
include("scr/aux_fun.jl")
include("scr/pwr_fun.jl")
include("scr/gas_fun.jl")

# simulation settings
set = Dict(
# gas and power networks
:gas_case => "case_48", :pwr_case => "pglib_opf_case118_ieee.m",
# control horizon data
:T => 5, :k_t => [1 4 7 10 13],
# wind data
:w̅ => 100, :σ => args["σ²"],
# number of uncertainty samples
:Ns => args["Ns"],
# constraint violation tolerance
:ε̅_p => args["ε̅_p"], :ε̅_g => args["ε̅_g"],
# variability parameters
:α_ϱ => args["α_ϱ"], :α_ϑ => args["α_ϑ"], :α_ψ => args["α_ψ"],
# netowrk topology opt data
:b_valve_edge => [21 30], :M̅ => 10e3, :M̲ => 10e3,
)

## solve CC OPF
net = load_pwr_data(set)
# wind power forecast
ξ = wind_forecast(set)
# solve OPF models
sol_opf = solve_chance_constrained_OPF(net,set,ξ)
# save gas consumption at gas-power coupling points
CSV.write("data/gas_uncertainty.csv",  Tables.table(net[:Λ]*net[:M_p]*sol_opf[:P][set[:T]]), writeheader=false)

## solve CC gas network optimization
net             = load_gas_data(set)
# solve non-convex gas network optimization
sol_non_convex  = non_convex_opt(net)
# linearize Weymouth equation
lin_res         = linearization(net,sol_non_convex[:model])
# solve the linearized model
sol_lin         = linearized_opt(net,lin_res)
# optimize chance-constrained gas network optimization
sol_cc  = chance_constrained_gas_opt(net,lin_res,set)
ofs_res_sto = out_of_sample(sol_cc,net,set,ξ)

## display some and export all results
println("expected cost:                \$$(round(sol_cc[:cost]/1000,digits=1))×10^3")
println("pressure variability measure:   $(round(sol_cc[:ϱ_std_tot]/1000,digits=1))×10^3")
open("$(outdir)/summary.txt","a") do io
   println(io,"simulation settings:")
   println(io,"α_ϱ = $(args["α_ϱ"]), α_ϑ = $(args["α_ϑ"]), α_ψ = $(args["α_ψ"]), σ = $(args["σ²"])")
   println(io,"Ns = $(args["Ns"]), ε̅_p = $(args["ε̅_p"]), ε̅_g = $(args["ε̅_g"])")
   println(io,"results:")
   println(io,"exp_cost = ", sol_cc[:cost])
   println(io,"variability = ", sol_cc[:ϱ_std_tot])
   println(io,"mean_comp_deployment = ", ofs_res_sto[:mean_compr_deployment])
   println(io,"mean_valve_deployment = ", ofs_res_sto[:mean_valve_deployment])
   println(io,"mass_con_violation_exp = ", ofs_res_sto[:mass_violation_exp])
   println(io,"press_con_violation_exp = ", ofs_res_sto[:press_violation_exp])
end
