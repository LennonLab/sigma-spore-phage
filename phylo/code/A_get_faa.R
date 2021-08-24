library(here)
library(tidyverse)
library(seqinr)
# This code runs blast using a wsl linuxx system on windows

# This code collects the bacterial sigma factors as described by Burton et al. 2019
# > Sigma homologs to SigA, SigH, and SigM were identified in the genomes of
# > 24 model organisms (listed in Table S6) using BLAST version 2.2.31 and an 
# > E value threshold of 1E-2 and compared to all SigN homologs found in the database.






# Get bacterial faa -------------------------------------------------------
# Get bacterial faa files specified by Burton et al. 2019

d.burton <- read_csv(here("phylo","data","Burton_S6.csv"))

#fix names to remove spaces
d.burton <- d.burton%>%
  mutate(organism_name=str_replace_all(organism_name," ","_"))

if (!dir.exists(here("phylo","data","bacteria_faa/"))){
  dir.create(here("phylo","data","bacteria_faa/"))
} 
setwd(here("phylo","data","bacteria_faa/"))

for (i in 1:nrow(d.burton)){

  file.name <- paste0(str_remove(d.burton$ftp_path[i],".*/"),"_protein.faa.gz")
  
  ftp <- paste0(d.burton$ftp_path[i],"/",file.name)
  
  download.file(ftp, destfile = paste0(d.burton$organism_name[i],".faa.gz"))
}

##############
# make blast database of 
setwd(here("phylo","data", "bacteria_faa/"))

# wsl <- "gunzip -c *.gz | makeblastdb -dbtype prot -in -"
system("wsl gunzip -c *.gz | makeblastdb -dbtype prot -in - -out bact24 -title bact24  -parse_seqids")

# blast sigmas against protein db

#make folder for results
if (!dir.exists(here("phylo","data","blast_bact/"))){
  dir.create(here("phylo","data","blast_bact/"))
} 
setwd(here("phylo","data","blast_bact/"))

#list of query genes
query <-  dir(here("phylo","data","Bs168_sig_faa"),pattern = ".faa$",full.names = F) 

# blast Bs sigmas agains 24 bacteria
for(q in query){
  
  #get gene name
  # gene <- str_extract(q,"/sig.")%>%str_remove("/")
  gene <- str_extract(q,"sig.")
  # blast command for linux
  wsl <- paste0("wsl blastp -query ../Bs168_sig_faa/",q," -db ../bacteria_faa/bact24 -evalue 0.01 -outfmt 6  -out ",gene,"_hits.tab" )
  # do blast
  system(wsl)
}



#collect all results 
hits<-  dir(".",pattern = ".hits.tab$",full.names = T) 
#column names
fmt6.names <- c("qseqid" ,"sseqid" ,"pident" ,"length", "mismatch" ,"gapopen" ,'qstart' ,'qend' ,"sstart" ,"send" ,"evalue" ,"bitscore")

d.hits <-  tibble()
#read in results
for (i in hits){
  d.hits <- read_tsv(i, col_names = fmt6.names) %>% 
    mutate(gene = str_extract(q,"sig.")) %>% 
    bind_rows(d.hits, .)
}

# how many hits ?
hits <- 
  d.hits%>%
  filter(!duplicated(sseqid))%>%
  select(sseqid)
#243


# get fasta of hits from balstdb
setwd(here("phylo","data", "bacteria_faa/"))
system("wsl > ../bacterial_sigma.faa")
wsl <- paste("wsl blastdbcmd -entry ", 
            hits$sseqid,
             " -db bact24 >> ../bacterial_sigma.faa")
sapply(wsl, system)

#  Get bacterial features for blast hits ----------------------------------


if (!dir.exists(here("phylo","data","bacteria_features/"))){
  dir.create(here("phylo","data","bacteria_features/"))
} 
setwd(here("phylo","data","bacteria_features/"))

d.features <- tibble()

for (i in 1:nrow(d.burton)){
  
  file.name <- paste0(str_remove(d.burton$ftp_path[i],".*/"),"_feature_table.txt.gz")
  
  ftp <- paste0(d.burton$ftp_path[i],"/",file.name)
  
  new.file.name <- paste0(d.burton$organism_name[i],"._features.tsv.gz")
  download.file(ftp, destfile = new.file.name)
  
  d.features <- read_tsv(new.file.name) %>% 
    filter(`# feature` == "CDS") %>% 
    select(assembly, genomic_accession, product_accession, name, symbol) %>% 
    mutate(sp = d.burton$organism_name[i]) %>% 
    filter(product_accession %in% hits$sseqid) %>% 
    bind_rows(d.features,.)
}


write_csv(d.features, here("phylo/data/bacterial_features.csv"))

#cleanup
setwd("..")
unlink(here("phylo","data","bacteria_features/"), recursive = T)
unlink(here("phylo","data","bacteria_faa"), recursive = T)


# combine bacterial and viral sigmas for alignment ------------------------------------------

# load baterial fastas collected above
bact.fa <- read.fasta(here("phylo","data","bacterial_sigma.faa")) 

# load viral sigmas from vog HMM analysis
load(here("vogdb/data/vog_sigma_clean_Whost.RData"))

# keep only viral sigmas of bacillus phages
phage.fa <- d.faa %>% 
  filter(phylum == "Firmicutes") %>% 
  filter(str_detect(family.etc, regex("bacillus", ignore_case = T)))

# write a combined fasta with minimal header
write.fasta(sequences = c(getSequence(bact.fa), phage.fa$seq), 
            names = c(paste0(names(bact.fa),"-bacteria"),
                      paste0(phage.fa$protein,"-phage")),
            file.out = here("phylo/data/sigmas_to_align.faa"))
