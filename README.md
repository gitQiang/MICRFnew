# MICRF
This repository contains MICRF (Mutation Integration with Conditional Random Field), an integrative network method to discover more risk genes, described in .

Motivation:
Under current relative small exome sequencing samples, only a few risk genes were discovered by statistical methods based on {\it de novo} mutations. Several existing methods showed that biological networks are helpful for improving the power of risk discovery. However, existing methods are limited by their selected networks and modeled biological networks as unweighted. To overcome the limitations of existing models and maximum the benefit of network information, we developed a new integration method based on a conditional random field model, MICRF (Mutation Integration with Conditional Random Field).  

Implementation: R and matlab.

Data integration:

1. de novo mutations from 5542 trios of neuron development disorders (NDD)

2. Biological networks: three protein-protein interaction networks (STRING, HPRD+STRING, iRefIndex) and two co-expression networks (Pearson correlation cutoff based: CORR and neighbors in top 5 ranked correlation: CoEXP); Combined networks: CoEXP and PrePPI, called CoPrePPI

#Dependences:

Matlab toolbox: Mark Schmidt, UGM: https://www.cs.ubc.ca/~schmidtm/Software/UGM.html

R packages: WGCNA for building co-expression modules; igraph for network edge betweenness centrality

#Related methods:
TADA: is a Bayesian hierarchical model for finding statistical significance risk genes, which incorporates de novo mutations, inherited rare variants, and variants identified within case-control data. It is downloaded from http://wpicr.wpic.pitt.edu/WPICCompGen/TADA/TADA_homepage.htm. Origin TADA used the fraction of LOF and damaging missense mutations to estimate the mutation rates of LOF and damaging missense based on gene level mutation rates, we changed it directly deal with mutation type specific mutation rate.

DAWN: a network method based on hidden Markov random field model to label risk genes based on neighbor information. It is got under authors' request. Origin DAWN contains several bugs to generate NAN results, we fixed them based on minimum changes rule to eliminate NAN and avoid exiting by any exception. 

MAGI: based on a combinatorial optimization algorithm which simultaneously integrated PPIs and gene expression data to discover modules enriched for de novo mutations.

#Data
De novo mutation lists:

#File list:
DDD_denovo_mutations.R: collected DDD de novo mutation lists

coexp.R: build co-expression networks: CORR and CoEXP. They are built based on BrainSpan Microarray expression data.

#Useage:

#Input file format:


# Contacts:
Yufeng Shen: ys2411@cumc.columbia.edu

Qiang Huang: qh2159@cumc.columbia.edu
