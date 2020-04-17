# MCPCA_PopGen

## Introduction
Low to medium depth sequencing is cost-effective and allows researchers to increase sample size at the expense of lower accuracy for genotype calling. To incorporate uncertainties and maintain statistical power in downstream analysis, we introduce MCPCA_PopGen to analyze low depth sequencing data. The method uses dosages rather than genotypes to account for uncertainties in genotype calling. It further optimizes the choice of nonlinear transformations of dosages to maximize the Ky-Fan norm of the covariance matrix.

MCPCA_PopGen is an open-source package. The source code of MCPCA is provided by Soheil Feizi using Matlab (available [here](https://github.com/SoheilFeizi/MCPCA)). To make it easier to install and implement, we write the entire package MCPCA_PopGen in Julia language. 

The package includes the following files:
- main.jl: an example about applying MCPCA_PopGen to dosage data.
- MCPCA_PopGen.jl, Discretize.jl : the MCPCA_PopGen method
- MCPCA_sample_disc_wrapper.jl, MCPCA_sample_disc.jl, utils.jl: the MCPCA method in Julia lanugage.
- getJenksBreaks.R, jenksBrks.c: get Jenks breaks; ported from R package BAMMtools.
- DosageGenotype.txt: dosage data.
- ms: use _ms_ simulator.

## Example

### Prepare the genotype dosages

The genotype dosage data are calculated as the posterior mean of the genotype under additive coding. With values 0, 1 and 2 assigned to the number of minor alleles N in each genotype, the dosage is the conditional mean of N given the data. To call the posterior probability of the genotype likelihood, we use the program ANGSD. An example code
```{angsd}
./angsd -bam bam.filelist -GL 1 -out outfile -doMaf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doGeno 8 -doPost 1 -postCutoff 0.95 
```
gives output like this:
```{angsd}
chr1	762273	0.000000	0.000934	0.999066	0.000000	0.999998	0.000002	...
chr1	880238	0.999789	0.000211	0.000000	0.999579	0.000421	0.000000	...
chr1	886788	0.000000	0.007465	0.992535	0.984030	0.015970	0.000000	...
chr1	887560	0.998309	0.001691	0.000000	0.986626	0.013374	0.000000	...
```
- Column 3 - 5: Genotype likelihood for the first individual
- Column 6 - 8: Genotype likelihood for the second individual

The genotype dosages are calculated using equation DS = E( N|Data ) = Pr(1|Data) + 2Pr(2|Data).

### Apply MCPCA to genotype dosages

#### Load packages and include related functions
```{julia}
using Discretizers, DelimitedFiles, Statistics, StatsBase, Random, LinearAlgebra, Clustering, Distributions, RCall
include("./Discretize.jl")
include("./utils.jl")
include("./MCPCA_sample_disc.jl")
include("./MCPCA_sample_disc_wrapper.jl")
include("./MCPCA_PopGen.jl")
```
#### Set parameters
```{julia}
## number of MCPCs
kq = 10
````

#### Input data
```{julia}
DS = readdlm("DosageGenotype.txt");
n, p = size(DS)
phimat_ds = MCPCA_PopGen(DS, kq, discretize_method="none");
D, V = eigen( phimat_ds' * phimat_ds ./n);
D = sort(D, rev=true);
```

#### Use equal width binning, equal frequency binning, and Jenks binning
```{julia}
phimat_int = MCPCA_PopGen(DS, kq, discretize_method="interval");
phimat_freq = MCPCA_PopGen(DS, kq, discretize_method="counts");
phimat_jenks = MCPCA_PopGen(DS, kq, discretize_method="Jenks");
```

The principal components of phimat_int, phimat_freq, or phimat_jenks are the optimized maximally correlated principal components and can be used to visualize the data.

## Reference:
- Miao Zhang, Yiwen Liu, Hua Zhou, Jin Zhou, and Joseph Watkins. (2019). A novel non-linear dimension reduction approach to infer population structure for low-coverage sequencing data.
- Soheil Feizi and David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471
- Rabosky DL, Grundler M, Anderson C, Title P, Shi JJ, Brown JW, Huang H, Larson JG. BAMM tools: an R package for the analysis of evolutionary dynamics on phylogenetic trees. Methods in Ecology and Evolution. 2014 Jul;5(7):701-7.
- Korneliussen, T. S., Albrechtsen, A., and Nielsen, R. (2014). ANGSD: analysis of next generation sequencing data. BMC bioinformatics, 15(1):356.
