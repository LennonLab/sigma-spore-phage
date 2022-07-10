For this analysis I downloaded TIGRFAM files from:  
<https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/>  
(downloaded 21/Nov/2020)

Extracted HMMs (*TIGRFAMs_15.0_HMM.tar.gz*) into folder *TIGR/data/tigr_hmm*  
Extracted info (*TIGRFAMs_15.0_INFO.tar.gz*) into folder *TIGR/data/tigr_info*  

Then ran code in *code/A_parse_tigr_info*.  
The final step of this code is making an HMM database. The code for this uses HMMER installed in Windows Subsystem for Linux being called from R.

Main analysis is in *B_hscan_tigr.R*, comparing each phage protein to the profile database.

Finaly, analysis results are parse in *C_parse_hscan.R*.


