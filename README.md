# MICRFnew
This repository contains MICRF (Mutation Integration with Conditional Random Field), an integrative network method to discover more risk genes, described in .

Problem:
Under current relative small sample size, only a few risk genes were discovered by {\it de novo} mutations. Biological networks are helpful for improving the power.  

Model:

Implementation:
It is implemented based on R and matlab.


Data integration:
1. de novo mutations for DDD
2. Networks: three protein-protein interaction networks (STRING, HPRD, iRefIndex) and two co-expression networks (Pearson correlation: CORR and top 5 ranked neighbors: CoEXP)
3. Combined networks: CoEXP and PrePPI, called CoPrePPI

#Dependences:

R packages: WGCNA for co-expression modules; igraph for network betweenness centrality

Matlab toolbox: UGM: https://www.cs.ubc.ca/~schmidtm/Software/UGM.html

#Related methods:
TADA
DAWN
MAGI
//HotNet2


#File list:
DDD_denovo_mutations.R: collected DDD de novo mutation lists

coexp.R: build co-expression networks




#Simulations:



