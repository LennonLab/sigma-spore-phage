library(tidyverse)
library(here)
# Downloads come from http://vogdb.org/download (24/Nov/2020)

# import list of vogs
d.vog <- read_tsv(here("vogdb","vogdb_downloads","vog.annotations.tsv"))%>%
  rename("GroupName"="#GroupName")

# find vogs of sigma factors
d.vog.sigma <- d.vog%>%
  filter(str_detect(ConsensusFunctionalDescription, regex("sigma",ignore_case = T)))%>%
  filter(str_detect(ConsensusFunctionalDescription, regex("factor",ignore_case = T)))%>%
  filter(!str_detect(ConsensusFunctionalDescription, regex("anti",ignore_case = T)))


# get the sequences out of tar.gz file
setwd(here("vogdb","vogdb_downloads"))
for(i in d.vog.sigma$GroupName){
  
  wsl <- paste0("wsl tar -zxvf vog.faa.tar.gz ",i,".faa")
  system(wsl)
}

#move files
f <- dir(pattern = ".faa$")
file.copy(f, here("vogdb","data","vog_sigma_faa"))
file.remove(f)
setwd(here())

