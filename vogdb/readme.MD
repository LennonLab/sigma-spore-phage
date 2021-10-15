# Retrieval and analysis of phage-encoded sigma factors from the  Virus Orthologous Groups (VOG) database

**Some of the code in this section relies on the Windows Subsystem for Linux being called from R.**

VOG data was downloaded from <http://vogdb.org/download> (24/Nov/2020) and placed in the *vogdb/vogdb_downloads* folder. 
* vog.annotations.tsv
* vog.species.list
* vog.members.tsv
* vog.faa.tar.gz

Virus-Host database data downloaded from <ftp://ftp.genome.jp/pub/db/virushostdb/> (24/Nov/2020) and saved to *vogdb/vogdb_downloads/virushostdb.tsv*


* A_get_vog_files.R  
** find VOGs of sigma factors  

* B_collect_faa.R  
** get fastas of phage-encoded sigmas

* B2_phage_host.R  
match VOG phages to host using Virus-Host database

* C_hscan_tigr.R  
* D_parse_hscan.R
* E_align_to_hmm.R 
* contingency.R
* parse_hmmer_tbl.R
* plot_hosyPhylu_vFam.R

