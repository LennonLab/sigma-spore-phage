library(tidyverse)
library(here)
# Downloads come from http://fileshare.csb.univie.ac.at/vog/

# for reproducibility I will use the current latest version (Jan/6/2022)
vog_version <- "209"
url_vog <- paste0("http://fileshare.csb.univie.ac.at/vog/vog",vog_version,"/")

# for latest version uncomment theses lines
# latest_vog <- readLines(url("http://fileshare.csb.univie.ac.at/vog/latest/release.txt"))
# url_vog <- paste0("http://fileshare.csb.univie.ac.at/vog/vog",latest_vog,"/")

#make folder for downloads
dir_vog <- here("vogdb/data/vogdb_downloads/")
if (!dir.exists(dir_vog)){
  dir.create(dir_vog)
}

#download the relavent files
download_files <- c("vog.annotations.tsv.gz", 
                    "vog.faa.tar.gz", 
                    "vog.members.tsv.gz", 
                    "vog.species.list" )

for (f_i in download_files){
  if(!file.exists(here(dir_vog, f_i))){
    download.file(url = paste0(url_vog,f_i),
                  destfile = here(dir_vog, f_i))
  }
}


# import list of vogs
d.vog <- read_tsv(here(dir_vog,"vog.annotations.tsv.gz"))%>%
  rename("GroupName"="#GroupName")

# find vogs of sigma factors
d.vog.sigma <- d.vog%>%
  filter(str_detect(ConsensusFunctionalDescription, regex("sigma",ignore_case = T)))%>%
  filter(str_detect(ConsensusFunctionalDescription, regex("factor",ignore_case = T)))%>%
  filter(!str_detect(ConsensusFunctionalDescription, regex("anti",ignore_case = T)))


# get the sequences out of tar.gz file
  #make folder for fastas
  dir_faa <- here("vogdb","data","vog_sigma_faa")
  if (!dir.exists(dir_faa)){
    dir.create(dir_faa)
  }
  
untar(tarfile = here(dir_vog, "vog.faa.tar.gz"), 
      files = paste0(d.vog.sigma$GroupName,".faa"), 
      exdir = dir_faa)
