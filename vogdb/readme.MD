# Retrieval and analysis of phage-encoded sigma factors from the  Virus Orthologous Groups (VOG) database

**Some of the code in this section relies on the Windows Subsystem for Linux being called from R.**

VOG data was downloaded from <http://vogdb.org/download> (24/Nov/2020) and placed in the *vogdb/vogdb_downloads* folder.

* *vog.annotations.tsv*
* *vog.species.list*
* *vog.members.tsv*
* *vog.faa.tar.gz*

Virus-Host database data downloaded from <ftp://ftp.genome.jp/pub/db/virushostdb/> (24/Nov/2020) and saved to *vogdb/vogdb_downloads/virushostdb.tsv*


* A_get_vog_files.R  
    * find VOGs of sigma factors  

* B_collect_faa.R  
    * get fastas of phage-encoded sigmas

* C_phage_host.R  
    * match VOG phages to host using Virus-Host database. 
    * summary of sigma factor presence by host and viral taxa.
    * summary of sigma factor number by phage genome in various host and phage taxonomic levels.

* D_hscan_tigr.R  
  *This script relies on the windows subsytem for Linux (wsl) on which hmmer has been installed*
    * Preform HMMscan for each of the VOG sigma factors against database of TIGRfam sigma database (from script A).  


* E_parse_hscan.R  
  *This script requires [helper functions](https://github.com/arendsee/rhmmer/blob/master/R/parse.R) sourced from parse_hmmer_tbl.R.*  
    * Classify VOG sigmas by HMMscan
    * Plot TIGR by host phylum
    * Plot TIGRfams for Firmicutes phages genes

* F_align_to_hmm.R 
    * trying it out 


