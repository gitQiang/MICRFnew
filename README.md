# MICRFnew
This repository contains MICRF (Mutation Integration with Conditional Random Field), an integrative network method to discover more risk genes, described in .

It is implemented based on R and matlab.


Data integration:
1. de novo mutations for DDD
2. Networks: three protein-protein interaction networks (STRING, HPRD, iRefIndex) and two co-expression networks (Pearson correlation and top 5 ranked neighbors)
3. Combined networks: HPRD  + coexp + GGM (based on glasso R packages)

#Dependences:

R packages: WGCNA
Matlab toolbox: UGM: https://www.cs.ubc.ca/~schmidtm/Software/UGM.html

#Related methods:
TADA
DAWN
MAGI
HotNet2


#File list:
DDD_denovo_mutations.R: collected DDD de novo mutation lists

coexp.R: build co-expression networks




#Simulations:



