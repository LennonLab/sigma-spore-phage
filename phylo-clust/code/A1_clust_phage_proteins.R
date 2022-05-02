library(here)
library(tidyverse)
library(seqinr)
library(foreach)

# load viral sigmas from VOG
load(here("vogdb/data/vog_sigma_clean_Whost.RData"))


# cluster with cd-clust at different levels ----
if(!dir.exists(here("phylo-clust/data/clust_phage"))){
  dir.create(here("phylo-clust/data/clust_phage"))
}

# > make multi fastsa ------------------------
write.fasta(sequences = d.faa$seq,
            names = d.faa$protein, 
            file.out = here("phylo-clust/data/clust_phage/vog_phage.faa"))

cdhit.cmd <- "module load cd-hit; cd-hit"

v_aS = 1 #seq(0.5,1,0.05) # aS does not matter here
v_c = seq(0.65,1,0.01)

d.clusts <- tibble(aS = rep(v_aS, length(v_c)),
                   C = rep(v_c, each= length(v_aS)),
                   n.clusters = NA)
foreach(i = seq(1:nrow(d.clusts))) %do%
  {
    
    cmd <- 
      paste(cdhit.cmd, 
            "-aS", d.clusts$aS[i],
            "-c", d.clusts$C[i],
            "-d 0 -g 1 -sc 1", 
            "-i", here("phylo-clust/data/clust_phage/vog_phage.faa"), 
            "-o", here("phylo-clust/data/clust_phage", 
                       paste0("clust_aS",d.clusts$aS[i],"_c",d.clusts$C[i]))
            
      )
    
    cd.x <- system(cmd, intern = T)
    
    # extract number of clusters
    d.clusts$n.clusters[i] <- 
      cd.x %>% str_detect("clusters") %>% cd.x[.] %>% 
      str_remove(.,".*finished *") %>% parse_number()
    
    # -aS  alignment coverage for the shorter sequence, default 0.0
    #      if set to 0.9, the alignment must covers 90% of the sequence
    # -c   sequence identity threshold, default 0.9
    #      this is the default cd-hit's "global sequence identity" calculated as:
    #      number of identical amino acids in alignment
    #      divided by the full length of the shorter sequence
    # -d	length of description in .clstr file, default 20
    # 	if set to 0, it takes the fasta defline and stops at first space
    # -g	1 or 0, default 0
    #      by cd-hit's default algorithm, a sequence is clustered to the first 
    #      cluster that meet the threshold (fast cluster). If set to 1, the program
    #      will cluster it into the most similar cluster that meet the threshold
    #      (accurate but slow mode)
    #      but either 1 or 0 won't change the representatives of final clusters     
    # -sc	sort clusters by size (number of sequences), default 0, output clusters by decreasing length
    # 	if set to 1, output clusters by decreasing size
    
  }



d.clusts %>% 
  ggplot(aes(x = C, y = n.clusters))+
  geom_line()+
  geom_point(size=2, shape=21, fill = "white")+
  geom_text(aes(label = C), nudge_y = 10, size =3 )

################
#
# Using C = 0.95
#
################

# parse cd-hit-est clstr ----
clstr <- read_lines(here("phylo-clust/data/clust_phage",paste0("clust_aS",1,"_c",0.95,".clstr")))

# initialize
d.clstr <- tibble()
n.clstr = 0

for(i in 1:length(clstr)){
  
  # assign clstr number
  n.clstr <- if_else(str_detect(clstr[i],"^>Cluster",), 
                     parse_number(clstr[i]),n.clstr)
  #row with clstr number has no more data
  if(str_detect(clstr[i],"^>Cluster")) next
  
  n.member <- str_remove(clstr[i], "\t.*")
  acc <- str_remove(clstr[i],".*, >") %>% str_remove("\\.\\.\\..*")
  aa <- str_extract(clstr[i],"\t.*aa") %>% parse_number()
  is.rep <- str_detect(clstr[i], "\\*$")
  perc.id <- if_else(is.rep, 100,
                     str_extract(clstr[i], "at.*%") %>% parse_number())
  
  
  d.clstr <- tibble(protClstr = n.clstr,
                    protein = acc,
                    aa = aa,
                    perc.id_protClstr = perc.id ,
                    protClstr_rep = is.rep) %>% 
    bind_rows(d.clstr, .)
  
}

# add cluster data to fasta table
d.faa <- left_join(d.faa, d.clstr, by = "protein" )

# replace representatives if needed -----------

# replacment function
swap_cluster_rep <- function(d, phage, cluster){
  current_rep <- which(d$protClstr == cluster &
                         d$protClstr_rep)
  new_rep <- which(str_detect(d$sp, phage) &
                     d$protClstr == cluster)
  d$protClstr_rep[current_rep] <- FALSE
  d$protClstr_rep[new_rep] <- TRUE
  return(d)
}

# > SPO1 -------
spo1.clst <- d.faa %>% 
  filter(str_detect(sp, "SPO1")) %>% pull(protClstr)
d.faa %>% 
  filter(protClstr %in% spo1.clst) %>% view
#all double sequence clusters with Camphawk phage being the other member

d.new <- d.faa
for(i in spo1.clst){
  d.new <- swap_cluster_rep(d.new, "SPO1", i)
}
#check
d.new %>% 
  filter(protClstr %in% spo1.clst) %>% view

# > SP-10 -----------
sp10.clst <- d.faa %>% 
  filter(str_detect(sp, "SP-10")) %>% pull(protClstr)
d.faa %>% 
  filter(protClstr %in% sp10.clst) %>% view
# all in single sequence "clusters"

# > Eldridge -------------
eld.clst <- d.faa %>% 
  filter(str_detect(sp, "Eldridge")) %>% pull(protClstr)
d.faa %>% 
  filter(protClstr %in% eld.clst) %>% view
# all in single sequence "clusters"

# > Goe3 ------------
goe3.clst <- d.faa %>% 
  filter(str_detect(sp, "Goe3")) %>% pull(protClstr)
d.faa %>% 
  filter(protClstr %in% goe3.clst) %>% view
# all in single sequence "clusters"

# > T4 (Escherichia virus T4) -------------
t4.clst <- d.faa %>% 
  filter(str_detect(sp, "Escherichia virus T4")) %>% pull(protClstr)
d.faa %>% 
  filter(protClstr %in% t4.clst) %>% view
# T4 gp55 is a member of the largest cluster with 82 members. All are 100% identical
# making T4 cluster rep.

for(i in t4.clst){
  d.new <- swap_cluster_rep(d.new, "Escherichia virus T4", i)
}
#check
d.new %>% 
  filter(protClstr %in% t4.clst) %>% view


# > anthracis phage Fah -------------
Fah.clst <- d.faa %>% 
  filter(str_detect(sp, "Fah")) %>% pull(protClstr)
d.faa %>% 
  filter(protClstr %in% Fah.clst) %>% view
# in cluster with sequences from tw oher phages


for(i in Fah.clst){
  d.new <- swap_cluster_rep(d.new, "Fah", i)
}
#check
d.new %>% 
  filter(protClstr %in% Fah.clst) %>% view

# > anthracis phage Bcp1 -------------
bcp.clst <- d.faa %>% 
  filter(str_detect(sp, "Bcp1")) %>% pull(protClstr)
d.faa %>% 
  filter(protClstr %in% bcp.clst) %>% view

# also needs replacing
for(i in bcp.clst){
  d.new <- swap_cluster_rep(d.new, "Bcp1", i)
}
#check
d.new %>% 
  filter(protClstr %in% bcp.clst) %>% view


# Save results
save(d.new, file = here("phylo-clust/data/faa_table_clustered.RData"))
write_csv(d.new %>% select(-seq), file = here("phylo-clust/data/faa_table_clustered.csv"))
