library(tidyverse)
library(here)

# I want to find the sigma factor related tigr-fams
# downloaded TIGRFAM files from: 
# https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/
# (21/Nov/2020)

# Extracted HMMs (TIGRFAMs_15.0_HMM.tar.gz) into folder "TIGR/data/tigr_hmm"
# Extracted info (TIGRFAMs_15.0_INFO.tar.gz) into folder "TIGR/data/tigr_info"

# I will write the first 3 rows of each info file. E.g
# ID  rpmI_bact
# AC  TIGR00001
# DE  ribosomal protein bL35

# info=list.files(here("TIGR/data","tigr_info"), full.names = T)
# # i=info[2]
# 
# #make df to collect results
# tigr <- tibble()
# 
# for (i in info){
#   tmp <- 
#     read_lines(i,n_max = 3)%>%
#     as_tibble()%>%
#     separate(1,into = c("field","value"),sep = "  ")%>%
#     pivot_wider(names_from = 1,values_from=2)
#   
#   tigr <-   bind_rows(tigr,tmp)  
# }
# write_csv(tigr, here("TIGR/data","tigr_info_3lines.csv"))
tigr <- read_csv( here("TIGR/data","tigr_info_3lines.csv"))

# filter out all descripritions woth sigma factor
tigr.sigma <- tigr%>%
  filter(str_detect(DE, pattern = regex("sigma",ignore_case = T)))%>%
  filter(str_detect(DE, pattern = regex("factor",ignore_case = T)))%>%
  # remove anti-sigma facotr
  filter(str_detect(DE, pattern = regex("anti",ignore_case = T), negate = T))

write_csv(tigr.sigma, here("TIGR/data","tigr_info_sigma.csv"))
tigr.sigma <- read_csv( here("TIGR/data","tigr_info_sigma.csv"))

#copy all sigma related HMMS to separate folder
if (! dir.exists(here("TIGR/data","tigr_sigma_hmm"))){
  dir.create(here("TIGR/data","tigr_sigma_hmm"))
}
hmm.all <- list.files(here("TIGR/data","tigr_hmm"), full.names = T)
# i=tigr.sigma$AC[2]

for(i in tigr.sigma$AC){
  
  hmm.match <- 
    str_detect(hmm.all,i)%>%
    which()%>%
    hmm.all[.]
  
  file.copy(hmm.match, here("TIGR/data", "tigr_sigma_hmm"))
  
  
}
 # make HMM database in linux CI for hmmscan
setwd(here("TIGR/data/tigr_sigma_hmm/"))

#make a single file 
system("wsl cat *.HMM > sigma.fam")
#make database
system(("wsl hmmpress sigma.fam"))
#clean up
system(("wsl rm *.HMM"))

