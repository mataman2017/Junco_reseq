#!/bin/env julia

mappingfile = "pop.txt"

using DataFrames

tm = CSV.read(mappingfile, DataFrame, header=true, delim='\t')

taxonmap = Dict(row[:individual] => row[:species] for row in eachrow(tm))

using PhyloNetworks

genetrees = readMultiTopology("Trees_concat_astral_225_2-20_FILT_miss_auto_neutral_GC.txt");

q, t = countquartetsintrees(genetrees, taxonmap, showprogressbar=true)

using DataFrames

df_sp = writeTableCF(q, t)

CSV.write("tableCF_species.csv", df_sp)

d_sp = readTableCF!(df_sp)
