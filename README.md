# sigma-spore-phage
Analysis of sporulation-like sigma factors in phages

This repository contains several different analyses all related to the pulication **TBA**.
A second repo containing analysis of sporulation by flow-cytometry is at https://github.com/LennonLab/sigma-spore-phage-flow

* All the analyses were done in Rstudio. I used the ["here"](https://here.r-lib.org/) package throughout the repo to call files using relative paths. To start up Rstusio at the right local directory I created an empty start.R file that can be used if you have .R files asociated with Rstudio. 


## Bioinformatic analysis of phage sigma factors

This analysis spans over three interconnected folders. Analysis was run in the following order:  
1. TIGR - Retrieving profiles of bacterial-encoded sigma-factor families from TIGRfams.  
2. vogdb -  Retrieval and analysis of phage-encoded sigma factors from the  Virus Orthologous Groups (VOG) database.  
3. phylo -   Phylogenetic analysis of VOG and bacterial sigma factors


## Phage virulence
 * Virulence  - Growth curve based assay comparing phages SP10 WT and g120 knockout

## Expression with inducible sigma factors
 * RNAseq - Expression analysis following induced expression in *Bacillus subtilis* in exponential growth

## Combining bioinformatic and functional responses
 * synthesis

Within each folder a more detailed description is given in the readme file.
