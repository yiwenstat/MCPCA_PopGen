# MCPCA_PopGen

## Introduction
Low to medium depth sequencing is cost-effective and allows researchers to increase sample size at the expense of lower accuracy for genotype calling. To incorporate uncertainties and maintain statistical power in downstream analysis, we introduce MCPCA_PopGen to analyze low depth sequencing data. The method uses dosages rather than genotypes to account for uncertainties in genotype calling. It further optimizes the choice of nonlinear transformations of dosages to maximize the Ky-Fan norm of the covariance matrix.

MCPCA_PopGen is an open-source package. The source code of MCPCA is provided by Soheil Feizi using Matlab (available \href{here}{https://github.com/SoheilFeizi/MCPCA}). To make it easier to install and implement, we write the entire package MCPCA_PopGen in Julia language.

## Reference:
Miao Zhang, Yiwen Liu, Hua Zhou, Jin Zhou, and Joseph Watkins. (2019). A novel non-linear dimension reduction approach to infer population structure for low-coverage sequencing data.
Soheil Feizi and David Tse, Maximally Correlated Principal Component Analysis, arXiv:1702.05471
