# Phylogenetic analysis of VOG and bacterial sigma factors

**Some of the code in this section relies on the Windows Subsystem for Linux being called from R.**

**Some analyses in this section were carried out on the Indiana University [Carbonate](https://kb.iu.edu/d/aolp) cluster with batch jobs submitted using the Slurm workload manager.**

* A1_clust_phage_proteins.R
	* clustering phage sigma factors using cd-hit (v 4.7) as installed on Carbonate. The ourpose of this is to reduce redudnancy of closely related genes in the data.
	results written to phylo-clust/data/faa_table_clustered.csv.
	
* A2_get_bact_faa.R
    * collecting the bacterial sigma factors as described by [Burton et al. 2019](https://doi.org/10.1128/mbio.01899-19). This code uses [table S6](https://journals.asm.org/doi/suppl/10.1128/mBio.01899-19/suppl_file/mbio.01899-19-st006.xlsx) of the Burton et al. as input. Table was converted to csv file and saved as *phylo/data/Burton_S6.csv*
    * Collection of bacterial sigma factors relies on WSL with blast+ (I used v2.9.0). 
    * Bacterial and phage sigma factors ([from VOG](vogdb/data/sigmas_to_align.faa))are saved in one multifasta for alignment in next step.

* B_align-trim.sh
  This was executed on Carbonate interactive job
    * Multiple sequence alignment with MAFFT (v7.475) 
    * alignment trm with trimAl v1.4.rev22 build[2015-05-21]
    
* C_model-selection.sh
  This was executed on Carbonate interactive job
    * Substitution model test with modeltest-NG. Software installed from Docker image using Singularity (v3.6.4), as described in the code. Not sure about version. When run I got "modeltest x.y.z". The container naming suggests it is v.0.1.6.
    
* D_tree.sh
  This code relies on helper script [*batch_raxML-ng.sh*](code/batch_raxML-ng.sh) for batch job submission
    * check MSA with raxml-ng --check
    (swapped suplicate sequence headers using helper script *duplicated-seqs.R* )
    * phylogeny with  RAxML-NG (v0.9.0).
    
* Dsup_duplicated-seqs.R
		Checking if during phylogenetic analysis of trimmed sequences resulted in any duplicate sequences being removed.
		
* E_bootstrap_tree.sh
  This code relies on helper script [*batch_bootstrap.sh*](code/batch_bootstrap.sh) for batch job submission
    * bootstrapping phylogeny from previous step with  RAxML-NG (raxml-ng --bootstrap). Ran 500 bootstraps.
    * check for convergence of bootstraps (raxml-ng --bsconverge). Bootstrapping test converged after 450 trees.
    
* Plot tree
This required first collecting some data on bacterial genes  

    * F1_bact_hscan_tigr.R - classify bacterial sigmas in the same way done for phage sigmas in "vogdb/D_hscan_tigr.R/"
    Relies on the windows subsytem for Linux (wsl) on which hmmer has been installed  
  
    * F2_bact_parse-hscan_tigr.R - parse data from previous script  
    Requires function from "vogdb/code/parse_hmmer_tbl.R"  
  
    * F3_get_bact_taxa.R - retrieve taxonomy data on bacterial strains in tree.
    Relies on the windows subsytem for Linux (wsl) on which [Entrez Direct](https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/) has been installed  
  
Finally, plotting with ggtreeExtra:  
    * F4_tree-plot.R
    The tree in Newick format with Transfer Bootstrap estimates is saved as **phylo/data/PlotPhylogeny_wData.tree**
  
Final tweaks to the plot as shown in the manuscript was done on exported [pptx file](plots/sigma_circle_rooted.pptx).


* G_classify.R
	* combining classification of phage sigma factors from phylogenetic analysis with HMM analysis (TIGR) to finalize classification of phage sigma genes as sporulation-like or other.

