# sigma-spore-phage
Analysis of sporulation-like sigma factors in phages

This repository contains several different analyses all related to the publication **Phage-encoded sigma factors alter bacterial dormancy** (DOI:10.1128/msphere.00297-22).
A second repo containing analysis of sporulation by flow-cytometry is at https://github.com/LennonLab/sigma-spore-phage-flow

* All the analyses were done in Rstudio. I used the ["here"](https://here.r-lib.org/) package throughout the repo to call files using relative paths. To start up Rstusio at the right local directory I created an empty start.R file that can be used if you have .R files associated with Rstudio. 


## Bioinformatic analysis of phage sigma factors

This analysis spans over three interconnected folders. Analysis was run in the following order:  
1. vogdb -  Retrieval and analysis of phage-encoded sigma factors from the  Virus Orthologous Groups (VOG) database. This includes genomic clustering of phages and matching phages to host.

2. TIGR - Retrieving profiles of bacterial-encoded sigma-factor families from TIGRfams, and using these profiles to classify phage sigma factors. 
  
3. phylo-clust -   Phylogenetic analysis of (clustered) VOG genes along with bacterial sigma factors. This folder also includes the merging of TIGR and phylogeny classifications.


## Phage virulence
 * Virulence  - Growth curve based assay comparing phages SP10 WT and g120 knockout

## Expression with inducible sigma factors
 * RNAseq - Expression analysis following induced expression in *Bacillus subtilis* in exponential growth

## Combining bioinformatic and functional responses
 * synthesis

Within each folder a more detailed description is given in the readme file.
