# Retrieval and analysis of phage-encoded sigma factors from the  Virus Orthologous Groups (VOG) database

**Some of the code in this section relies on the Windows Subsystem for Linux being called from R.**

VOG data was downloaded from <http://vogdb.org/download> (v209; 6/JAN/2022) and placed in the *vogdb/data/vogdb_downloads* folder.

* *vog.annotations.tsv*
* *vog.species.list*
* *vog.members.tsv*
* *vog.faa.tar.gz*

Virus-Host database data downloaded from <ftp://ftp.genome.jp/pub/db/virushostdb/> (6/JAN/2022) and saved to *vogdb/dat/vogdb_downloads/virushostdb.tsv*


* A_get_vog_files.R  
    * find VOGs of sigma factors  

* B1_clst_getFastas.R  
    * Retrieving genomic sequences of phages for clustering
    
* B2_clst_CdHit.sh 
    * Clustering genomic sequences of phages with cd-hit
    
* B3_clst_parseCdHit.R  
    * Parse clustering results


* C_collect_faa.R  
    * get fastas of phage-encoded sigmas (proteins)

* D_match_phage_host.R  
    * match VOG phages to host using Virus-Host database. 
    * summary of sigma factor presence by host and viral taxa.