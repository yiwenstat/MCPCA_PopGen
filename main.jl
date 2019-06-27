## Use MCPCA_PopGen ##

## load packages
using Discretizers, DelimitedFiles, Statistics, StatsBase, Random, LinearAlgebra, Clustering, Distributions, RCall
## include related functions
include("./Discretize.jl")
include("./utils.jl")
include("./MCPCA_sample_disc.jl")
include("./MCPCA_sample_disc_wrapper.jl")
include("./MCPCA_PopGen.jl")

## parameter
## number of MCPCs
kq = 10

## input data
DS = readdlm("DosageGenotype.txt");
n, p = size(DS)
phimat_ds = MCPCA_PopGen(DS, kq, discretize_method="none");
D, V = eigen( phimat_ds' * phimat_ds ./n);
D = sort(D, rev=true);

## use equal width binning, equal frequency binning, and Jenks binning
phimat_int = MCPCA_PopGen(DS, kq, discretize_method="interval");

phimat_freq = MCPCA_PopGen(DS, kq, discretize_method="counts");

phimat_jenks = MCPCA_PopGen(DS, kq, discretize_method="Jenks");
