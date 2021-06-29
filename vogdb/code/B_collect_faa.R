library(tidyverse)
library(here)
library(seqinr)

f <- list.files(here("vogdb","data","vog_sigma_faa"),".faa$",full.names = T)

d.faa <- tibble()

for(i in f){
  faa <- read.fasta(i, whole.header = TRUE)
  
  d.faa <- tibble(header=getName(faa),seq=getSequence(faa))%>%
    mutate(protein=str_extract(header,".*\\.. "))%>%
    mutate(protein=str_remove(protein, " $"))%>%
    mutate(sp=str_extract(header,"\\[.*\\]"))%>%
    mutate(sp=str_remove(sp, "\\["))%>%
    mutate(sp=str_remove(sp, "\\]"))%>%
    mutate(description=str_extract(header," .*\\["))%>%
    mutate(description=str_remove(description, " \\["))%>%
    separate(col = protein,into = c("taxon","protein"),sep = "\\.",extra = "merge")%>%
    mutate(vog=str_extract(i, "VOG\\d*"))%>%
    bind_rows(d.faa)
}


#test for uniqness of protein
anyDuplicated(d.faa$protein)
# yes. how many?
duplicated(d.faa$protein)%>%sum() #11

#lets have a look
dup <- duplicated(d.faa$protein)%>%which()
dup <- d.faa$protein[dup]

d.faa%>%
  filter(protein %in% dup)%>%
  select(protein,vog)%>%
  table()

#all but 1 are duplicates found in the same vog (VOG00048)
#are they completely identical rows?
d.faa%>%
  filter(protein %in% dup)%>%
  distinct()

# remove duplicated row from main table
d.faa <- d.faa%>%distinct()

#still leaves 1 protein (YP_003084147.1) that is a member of 2 vogs (VOG00048, VOG26037)
# Cyanophage PSS2,	 type III RNAP sigma factor
dup <- duplicated(d.faa$protein)%>%which()
dup <- d.faa$protein[dup]
d.faa%>%
  filter(protein %in% dup)%>%
  select(protein,vog)%>%
  table()

# I will keep the row from the smaller vog (VOG26037) has only 2 sequences 
n.rm <- d.faa%>%
  mutate(n=row_number())%>%
  filter(protein %in% dup)%>%
  filter(vog=="VOG00048")%>%
  pull(n)%>%
  as.numeric()

d.faa <- d.faa[-n.rm,]

#test for uniqness of sequence
anyDuplicated(d.faa$seq) # yes. how many?
duplicated(d.faa$seq)%>%sum() #160
# are they from the same phage?
d.faa%>%
  select(sp, seq)%>%
  duplicated()%>%
  sum()
#no

# export fasta

setwd(here("vogdb","data","vog_sigma_clean"))
for(i in 1:nrow(d.faa)){
  file.name <- paste0(d.faa$protein[i],".faa")
  write.fasta(sequences = d.faa$seq[i],
              names =d.faa$header[i],
              file.out = file.name)
  
  #compress
  system(paste("wsl gzip",file.name))
}

# write table of sequences
write_csv(select(d.faa,-seq), here("vogdb","data","vog_sigma_clean.csv"))

# also save table as R data to maintain sequence
setwd( here("vogdb","data"))
save(d.faa,file="vog_sigma_clean.RData")
# load( here("vogdb","data","vog_sigma_clean.RData"))
