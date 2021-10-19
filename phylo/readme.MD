# Phylogenetic analysis of VOG and bacterial sigma factors

**Some of the code in this section relies on the Windows Subsystem for Linux being called from R.**

**Some analyses in this section were carried out on the Indiana University [Carbonate](https://kb.iu.edu/d/aolp) cluster with batch jobs submitted using the Slurm workload manager.**


* A_get_faa.R
    * collecting the bacterial sigma factors as described by [Burton et al. 2019](https://doi.org/10.1128/mbio.01899-19). This code uses [table S6[(https://journals.asm.org/doi/suppl/10.1128/mBio.01899-19/suppl_file/mbio.01899-19-st006.xlsx) of the Burton et al. as input. Table was converted to cdv file and save as *phylo/data/Burton_S6.csv*
    * Collection of bacterial sigma factors relies on WSL with blast+ (I used v2.9.0). 
    * Bacterial and phage sigma factors ([from VOG](vogdb/data/vog_sigma_clean_Whost.RData))are saved in one multifasta for alignment in next step.

* B_align-trim.sh
  This was executed on Carbonate interactive job
    * Multiple sequence alignment with MAFFT (v7.475) 
    * alignment trm with trimAl v1.4.rev22 build[2015-05-21]
    
* C_model-selection.sh
  This was executed on Carbonate interactive job
    * Substitution model test with modeltest-NG. Software installed from Docker image using Singularity (v3.6.4), as described in the code. Not sure about version. When run I got "modeltest x.y.z". Thcontainer naming suggests it is v.0.1.6.
    
* D_tree.sh
  This code relies on helper script [*batch_raxML-ng.sh*](phylo/code/batch_raxML-ng.sh) for batch job submission
    * check MSA with raxml-ng --check
    * phylogeny with  RAxML-NG (v0.9.0).
    
* E_bootstrap_tree.sh
  This code relies on helper script [*batch_bootstrap.sh*](phylo/code/batch_bootstrap.sh) for batch job submission
    * bootstrapping phylogeny from previous step with  RAxML-NG (raxml-ng --bootstrap). Ran 500 bootstraps.
    * check for convergence of bootstraps (raxml-ng --bsconverge). Bootstrapping test converged after 350 trees.
    
batch_bootstrap.sh
batch_parallel_bootstrap.sh
batch_raxML-ng.sh
C_model-selection.sh
D_tree.sh
D2_bootstrap_tree.sh
duplicated-seqs.R
E_tree-plot.R
F_ditances.R
get_taxa.R
sup_bact_hscan_tigr.R
sup_bact_parse-hscan_tigr.R
