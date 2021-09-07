library(here)
library(tidyverse)
library(cowplot)
library(scales)
source(here("vogdb/code/parse_hmmer_tbl.R"))

# parse hmmserch results
l.res <- list.files(here("phylo/data/hscan_bactXtigr/"), full.names = T)
# variable to store reults
hits.all <- tibble()


for (i in 1:length(l.res)){
  cur.hits <- read_tblout(l.res[i])
  
  # make row for 0 hits
  if(nrow(cur.hits)==0){
    #add row
    cur.hits[1,] <- NA
    # get protein name
    cur.hits$query_name[1] <- 
      str_remove(l.res[i],".*hscan_vogXtigr/")%>%
      str_remove(".txt")
    #assign no_hit fam
    cur.hits$description[1] <- "no_hit:no_hit"
    
    #assign non-significant values to E-values
    cur.hits <- cur.hits%>%
      mutate(across(where(is.numeric),~1))
  }
  
  hits.all <- bind_rows(hits.all,cur.hits)
  

}

# separate TIGR descriptions
hits.all <- hits.all%>%
  separate(description,into = c("id","description"),sep = ":" )



##############
# get 3 best hits for each protein
# defining best by sequence evalue()

best3 <- hits.all%>%
  group_by(query_name)%>%
  mutate(sequence_rate=-log10(sequence_evalue))%>%
  slice_max(sequence_rate, n=3)%>%
  # define hit order
  arrange(desc(sequence_rate))%>%
  mutate(place=row_number())%>%
  mutate(place=paste0("tigr.hit",place))%>%
  ungroup()

# hist(best3$sequence_rate, breaks = 100)
# all(!is.na(best3$sequence_rate))
# all(!is.na(best3$sequence_rate))

# find the top most specific hit
# if best hit is a general term:
  # "sigma70-ECF"
  # "SigBFG"
# look down the list to find a more specific one.
# if a more specific one does is not present use the general one
# in case the 2 general terms are listed and no 3rd term than use "SigBFG" as it is more specific


best3<- best3%>%
  pivot_wider(query_name, names_from=place,values_from=id)%>%
  add_column(tigr.hit=NA)%>%
  mutate(tigr.hit=as.character(tigr.hit))



g.terms <- c("sigma70-ECF","SigBFG")
# g.terms <- c("sigma70-ECF")

for(i in 1:nrow(best3)){
  if(!best3$tigr.hit1[i] %in% g.terms|
     is.na(best3$tigr.hit2[i])){
    best3$tigr.hit[i] <- best3$tigr.hit1[i]
    next
  }
  if(!best3$tigr.hit2[i] %in% g.terms|
     is.na(best3$tigr.hit3[i])){
    best3$tigr.hit[i] <- best3$tigr.hit2[i]
    next
  }
  best3$tigr.hit[i] <- best3$tigr.hit3[i]

}

# these rules did the job. 


best3%>%
  group_by(query_name)%>%
  group_by(tigr.hit)%>%
  summarise(n=n())%>%
  mutate(tigr.hit=fct_reorder(tigr.hit,n))%>%
  ggplot(aes(tigr.hit,n))+
  geom_col()+
  geom_label(aes(label=paste0("n=",n)), y=150)+
  coord_flip()+
  theme_cowplot()+
  xlab("best specific TIGRFAM match")  


#looking at e-values of best hits
best3.evalue <- hits.all%>%
  group_by(query_name)%>%
  slice_max(sequence_evalue, n=3)%>%
  # define hit order
  arrange(desc(sequence_evalue))%>%
  mutate(place=row_number())%>%
  ungroup()


best3.evalue%>%
  filter(id!="no_hit")%>%
  # pivot_wider(query_name, names_from=place,values_from=sequence_evalue)%>%
  ggplot(aes(x=place, y=-log10(sequence_evalue)))+
  geom_jitter(aes(color=sequence_score),alpha=0.7, width = 0.3)+
  geom_hline(yintercept = c(3), color="red")+
  geom_hline(yintercept = c(5), color="pink")+
  facet_wrap(~id)+
  scale_y_log10()+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")+
  theme_cowplot()


# all have an e-value lower than 1e-3 
# I don't think there is a problem


# Save data ---------------------------------------------------------------
best3 %>% 
  select(query_name,tigr.hit) %>% 
  write_csv(file = here("phylo/data/bacterial_sigma_tigr.csv"))
