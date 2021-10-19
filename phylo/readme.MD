# Phylogenetic analysis of VOG and bacterial sigma factors

**Some of the code in this section relies on the Windows Subsystem for Linux being called from R.**

**Some analyses in this section were carried out on the Indiana University [Carbonate](https://kb.iu.edu/d/aolp) cluster with batch jobs submitted using the Slurm workload manager.**


* A_get_faa.R
    * collecting the bacterial sigma factors as described by [Burton et al. 2019](https://doi.org/10.1128/mbio.01899-19). This code uses [table S6[(https://journals.asm.org/doi/suppl/10.1128/mBio.01899-19/suppl_file/mbio.01899-19-st006.xlsx) of the Burton et al. as input. Table was converted to cdv file and save as *phylo/data/Burton_S6.csv*
    * Collection of bacterial sigma factors relies on WSL with blast+ (I used v2.9.0). 
    * Bacterial and phage sigma factors ([from VOG](vogdb/data/vog_sigma_clean_Whost.RData))are saved in one multifasta for alignment in next step.

B_align-trim.sh
batch_bootstrap.sh
batch_parallel_bootstrap.sh
batch_raxML-ng.sh
C_model-selection.sh
D1_tree.sh
D2_bootstrap_tree.sh
duplicated-seqs.R
E_tree-plot.R
F_ditances.R
get_taxa.R
sup_bact_hscan_tigr.R
sup_bact_parse-hscan_tigr.R
