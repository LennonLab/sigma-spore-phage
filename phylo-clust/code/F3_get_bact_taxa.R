library(here)
library(tidyverse)
library(seqinr)

# Get bacterial taxa specified by Burton et al. 2019

d.burton <- read_csv(here("phylo-clust","data","Burton_S6.csv"), trim_ws = T) %>% 
  rename( assembly_accession = `# assembly_accession`)

#use e-utilities installed locally on carbonate
#sh -c "$(curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
eutils="/N/u/danschw/Carbonate/my_tools/edirect"
d <- tibble()

for (tx in d.burton$taxid){
  ranks <- system(paste0(eutils,"/efetch -db taxonomy -id ",tx," -format xml | ",
                         eutils,"/xtract  -pattern LineageEx -sep ',' -element Rank "), intern = T)
  lineage <- system(paste0(eutils,"/efetch -db taxonomy -id ",tx," -format xml | ",
                           eutils,"/xtract  -pattern LineageEx -sep ',' -element ScientificName "), intern = T)
  ranks <- unlist(str_split(ranks, pattern = ","))
  lineage <- unlist(str_split(lineage, pattern = ","))
  
  d <- tibble(tx,ranks,lineage) %>%  
    bind_rows(d,.)
}

d2 <- d %>% 
  filter(ranks != "clade", ranks != "no rank") %>% 
  pivot_wider(names_from = ranks, values_from = lineage)%>% 
  left_join(select(d.burton, assembly_accession, taxid, organism_name ), . ,
            by = c("taxid" = "tx"))
write_csv(d2,(here("phylo-clust","data","bacterial_taxonomy.csv")))
