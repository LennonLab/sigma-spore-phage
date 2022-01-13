library(tidyverse)
library(here)
library(rentrez)

# Reviewer asked we examined TIGRFAMs of bacterial sigma factors
# to make sure they do not include prophages

# sigma TIGFAMS
tigr.sigma <- read_csv(here("TIGR/data","tigr_info_sigma.csv"))
# get IDs 
ac_tigr <- tigr.sigma$AC


# #get seeds
# download.file(url = "https://ftp.ncbi.nlm.nih.gov/hmm/TIGRFAMs/release_15.0/TIGRFAMs_15.0_SEED.tar.gz",
#               destfile = here("TIGR/data/TIGRFAMs_15.0_SEED.tar.gz"))

# extract sigma seeds
untar(here("TIGR/data/TIGRFAMs_15.0_SEED.tar.gz"), 
      files = paste0(ac_tigr,".SEED"),
      exdir = here("TIGR/data/sigma_seeds"))


# Get seed headers
f.seed <- list.files(here("TIGR/data/sigma_seeds"), full.names = T)

d.seeds <- tibble()
# extract headers
for (tigr in ac_tigr){
  d.seeds <- 
    read_lines(here("TIGR/data/sigma_seeds", paste0(tigr,".SEED")), skip = 1) %>% 
      str_remove(" .*$") %>% 
      tibble(header=., tigr = tigr) %>% 
      filter(header != "//") %>% 
      separate(header, into = letters[1:4], sep = "\\|") %>% 
    bind_rows(d.seeds)
}

# get genomic data
# using edirect utility to find refseq identical protein (ipg)
d.seeds$genome <- NA 
for (i in 1:nrow(d.seeds)){
  print(i)
  if (! d.seeds$a[i] %in% c("GB", "SP")) next
  
  # Records removed
  if (d.seeds$b[i] %in% c("AAZ47680.1", "AAZ47880.1", "AAZ46269.1",
                          "AAZ47364.1","AAZ54493.1","AAZ56630.1",
                          "P17211/50-276")) next
  str_replace(d.seeds$b[i],"/.*","")
  
  wsl <- paste("wsl",
               "efetch -db protein  -id",
               str_replace(d.seeds$b[i], regex("\\.\\d$"), "") %>% #version independent 
                 str_replace(.,"/.*",""), # prblem in one of the rows
               '-format ipg | grep -m 1 \"RefSeq\\|INSDC\"')
  
  d.seeds$genome[i] <- system(wsl, intern = T)
  
}

#parse results
ipg <- c("ipg", "source", "nuc.acc","nuc.from","nuc.to", 
         "strand", "prot.acc", "name","organism","strain")
d.seeds <- d.seeds %>% 
  separate(genome, into = ipg, sep = "\t")
write_csv(d.seeds, here("TIGR/data/seed_genomic.csv"))
# d.seeds <- read_csv(here("TIGR/data/seed_genomic.csv"))


# PHASTER
# get prophage loci for genomes with seeds
# using URLAPI of phaster https://phaster.ca/instructions#urlapi

if (! dir.exists(here("TIGR/data/phaster"))){
  dir.create(here("TIGR/data/phaster"))
}

# there are seeds that come from the same genome
genome.acc <- unique(d.seeds$nuc.acc)
for (i in 1:length(genome.acc)){
  print(i)
  if ( is.na(genome.acc[i])) next
  
  url <- paste0("http://phaster.ca/phaster_api?acc=",
                genome.acc[i])
  dest <- here("TIGR/data/phaster",paste0(genome.acc[i],".json"))
  download.file(url = url, 
                destfile = dest)
  Sys.sleep(10) # prevents some errors
}

# parse phaster

d.phaster <- tibble()

for (acc in genome.acc){
  if(is.na(acc)) next
  x <- read_lines(here("TIGR/data/phaster",paste0(acc,".json")))
  
  x <- strsplit(x[1], split = '\\\\n') %>% unlist()
  
  # find data lines
  n_underline <- which(str_detect(x,"------"))
  if(is_empty(n_underline)) next
  if (n_underline >= length(x)-1) next
  
  phaster_cols <- strsplit(x[n_underline-1], split = "\\s+") %>% unlist()
  
  x.data <- read.table(text = x[(n_underline+1):(length(x)-1)]) %>% 
    as_tibble()
  colnames(x.data) <- phaster_cols[-1]
  
  x.data$acc <- acc
  
  d.phaster <- bind_rows(x.data, d.phaster)
}



# now we have the genomic postion of each seed 
# and the prophage predictions for the same region
# we can check if any seed protein is in a prophage region
d.seeds$in_prophage <- NA

for(i in which(!is.na(d.seeds$nuc.acc))){
  d.seeds$in_prophage[i] <- 
    d.phaster %>% 
    filter(acc == d.seeds$nuc.acc[i]) %>% 
    separate(REGION_POSITION, into = c("from", "to"), sep = "-", convert = T) %>% 
    select(from,to) %>% 
    mutate(in_prophage_from = d.seeds$nuc.from[i] >= from,
           in_prophage_to = d.seeds$nuc.to[i] <= to) %>% 
    mutate(in_prophage = in_prophage_from & in_prophage_to) %>% 
    pull(in_prophage) %>% any(.)
}


write_csv(d.seeds, here("TIGR/data/seed_genomic.csv"))
write_csv(d.phaster, here("TIGR/data/seed_phaster.csv"))

# how well did we do?
table(d.seeds$in_prophage, useNA = "a")
# FALSE  TRUE  <NA> 
#   365     7   184 

d.seeds %>% 
  group_by(tigr, in_prophage) %>% 
  summarise(n=n()) %>% 
  pivot_wider(names_from = in_prophage, values_from = n) %>% 
  view()

# which TIGRFAMS have prophage members?
d.seeds %>% 
  filter(in_prophage) %>% 
  pull(tigr) %>% 
  unique()

# All are from "TIGR02937", which is the most general sigma factor TIGRFAM!!

# Focus on sporulation sigmas ---------------------------------------------

# get IDs sporulation sigma
spore_tigr <- tigr.sigma %>% 
  filter(str_detect(ID, "sigG|sigF|sigmaE|sigmaK")) %>% 
  pull(AC)

d.seeds %>% filter(tigr%in% spore_tigr) %>% view()

# none are in prophage.
# only one seed is NA

# The NA sequence is:
# >OMNI|NTL02CB2537/3-229
# FLSYLWDVLGSALFLTGYVTGNTSFPKPLSDEEEKYYLDRLKKGDALAKDILVERNLRLVAHIVKKYSYPGKD
# IDDLISIGTVGLIKAIDSFDNSKGTRLATYAARCIENEILMLIRNNKKTKGEVYLQDPIGVDKEGNEISLMDI
# LSSDEDSIIEIVSTKIEVKKLYGKIDTCLRGREKKVIQMRYGLKDGRPRTQREIAGILNISRSYVSRIEKKAL
# KKLYKELN

# Blastp identifies this as 
# "RNA polymerase sporulation sigma factor SigK [Clostridium botulinum]"
# NC_009495.1 2679478-2680182 (-) WP_012047887.1	

#get phaster
url <- paste0("http://phaster.ca/phaster_api?acc=NC_009495.1")
dest <- here("TIGR/data/phaster/NC_009495.1.json")
download.file(url = url, 
              destfile = dest)

# parse phaster

  x <- read_lines(here("TIGR/data/phaster/NC_009495.1.json"))
  x <- strsplit(x[1], split = '\\\\n') %>% unlist()
  
  # find data lines
  n_underline <- which(str_detect(x,"------"))
  is_empty(n_underline) #FALSE
  n_underline >= length(x)-1 #FALSE
  
  phaster_cols <- strsplit(x[n_underline-1], split = "\\s+") %>% unlist()
  
  x.data <- read.table(text = x[(n_underline+1):(length(x)-1)]) %>% 
    as_tibble()
  colnames(x.data) <- phaster_cols[-1]
  
  x.data$acc <- acc
  view(x.data)
  
  # While this genome has 3 prophage regions
  # sigK (at 2679478-2680182) is not within them
  # being in a position greater than all 3
