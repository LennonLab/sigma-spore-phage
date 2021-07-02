library(here)
library(tidyverse)
library(cowplot)
source(here("vogdb/code/parse_hmmer_tbl.R"))

# parse hmmserch results
l.res <- list.files(here("vogdb/data/hscan_vogXtigr/"), full.names = T)
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


# need to solve the no hit rows which do not have taxon
# #separate query name
#   separate(query_name, 
#            into = c("taxon","protein"),
#            sep = "\\.",
#            extra="merge", remove = FALSE)


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
# if best hit is a genral term:
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
# SigBFG was completely excluded (always had a more specific term to take over)
# sigma70-ECF was retained only for proteins that had no other hit

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
  xlab("best specific TIGRFAM match")+
  ggsave2(filename = here("vogdb/figures/vogXtigr.png"),
          width = 7,height = 7)


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


# all have an e-value lower than 1e-3 exxpet for a singlr sigH hit
# I don't think there is a problem




################################################
# combining vogdb and virus-host db info on proteins

# load in the data made in previous script
load(here("vogdb","data","vog_sigma_clean_Whost.RData"))

#combining of data frames will be by protein number.
# extract it from hscan query_name
# remove leading taxon number from query name if present
#leaving protein id that starts with letter

best3 <- best3%>%
  mutate(protein=str_remove(best3$query_name, regex("^[0-9]*\\.?")))

##QC
# anyDuplicated(best3$protein) # 0
# is.na(best3$protein)%>%sum() # 0
# all(d.faa$protein %in% best3$protein) #TRUE
# all(best3$protein %in% d.faa$protein) #TRUE


d.faa <- 
  best3%>%
  select(-query_name)%>%
  full_join(d.faa, ., by="protein")


############################################
# Do tigr hits corespond to VOG?
############################################
d.faa%>%
  select(vog,tigr.hit,tigr.hit1,tigr.hit2,tigr.hit3 )%>%
  pivot_longer(-vog, names_to="hit", values_to="tigr.fam")%>%
  group_by(vog,hit,tigr.fam)%>%
  summarise(n=n())%>%
  filter(!is.na(tigr.fam))%>%
  ggplot(aes(x=vog, y=tigr.fam))+
  geom_tile(aes(fill=n))+
  facet_wrap(~hit)+
  scale_fill_viridis_b()+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=30,hjust = 1),
        legend.position = "bottom")+
  ggsave2(filename = here("vogdb/figures/vogXtigr_HitPlace.png"),
          width = 10,height =10)

# mostly no clear 1 to 1 correspondence

############################################
# parsing sigma factor type by 
# number of sigma factors in phage
###########################################

# add column of sigma factor content per genome
n.sig <- d.faa%>%
  group_by(taxon, `virus name`)%>%
  summarise(n.sig=n())%>%
  ungroup()

d.faa <- n.sig%>%
  select(-`virus name`)%>%
  right_join(.,d.faa, by="taxon")

###################################
# classes of sigma factors
sig.class <- read_csv(here("vogdb/data/tigr_Bs_sigma.csv"))

d.faa <- sig.class%>%
  select(ID,sig.class, sig.group)%>%
  left_join(d.faa,., by=c("tigr.hit"="ID"))

d.faa$sig.class[d.faa$tigr.hit=="no_hit"] <- "no_hit"
d.faa$sig.group[d.faa$tigr.hit=="no_hit"] <- "no_hit"


###################
# add sigma factor group classification
# (see pfam section)
load(here("pfam/data/phage_sigma_groups.RData"))
#parse out protein details to combine
d.group <- d.group%>%
  separate(query_name,into = c("taxon","protein"),
           sep = "\\.",remove = F,extra = "merge")%>%
  #fix parsing for no_hit rows
  mutate(protein=if_else(sigma.group=="no_hit",query_name,protein))%>%
  mutate(taxon=if_else(sigma.group=="no_hit","--",taxon))

# # QC
# all(d.faa$protein %in% d.group$protein)
# all(d.group$protein %in% d.faa$protein)

d.faa <-
  d.group%>%
  select(protein,sigma.group)%>%
  left_join(d.faa,.,by="protein")


##########################
# Plot by phylum
##########################
n.labs <- d.faa%>%
  group_by(protein,phylum)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(phylum)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))

tigr.spore=c("spore_sigmaK","spore_sigmaE","spore_sigF","spore_sigG")

p.phylum <- 
  d.faa%>%
  # arrange by frequency of sigma factors
  mutate(tigr.hit=fct_infreq(tigr.hit))%>%
  group_by(phylum,tigr.hit)%>%
  summarise(n=n())%>%
    group_by(phylum)%>%
    mutate(perc=n/sum(n))%>%
    ungroup()%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))%>%
  #add color to sporulation sigma factors
  mutate(spore.sigma=if_else(tigr.hit%in%tigr.spore,"spore","other"))%>%
    ggplot( aes(fct_rev(tigr.hit),perc, group = phylum)) + 
  #line for 0
  geom_hline(yintercept = 0, color="grey")+
  geom_col(aes(fill=spore.sigma)) + 
  geom_text(data=n.labs, aes(label=lab), x=2,y=0.5)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage sigma factor genes") +
  xlab("Sigma Factor TIGRFAM")+
  facet_grid(~phylum)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5),
        legend.position = "bottom")+
  coord_flip()+
  panel_border()

p.phylum+
  ggtitle("phage sigma factor type by host phylum")+
  ggsave2(here("vogdb","figures","tigr_sigma_HostPhylum.png"),
          width = 12,height = 8)
###
# same plot different way AND better!!
p.phylum <-
  d.faa%>%
  # arrange by frequency of sigma factors
  mutate(tigr.hit=fct_infreq(tigr.hit))%>%
  group_by(phylum,tigr.hit)%>%
  summarise(n=n())%>%
  ungroup()%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))%>%
  ggplot( aes(fct_rev(tigr.hit),n, group = phylum)) + 
  #line for 0
  geom_hline(yintercept = 0, color="grey")+
  geom_col(aes(fill=phylum)) + 
  ylab("Phage sigma factor genes") +
  xlab("Sigma Factor TIGRFAM")+
  # facet_grid(~phylum)+
  theme_cowplot()+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow = 2))+
  scale_fill_viridis_d()+
  coord_flip()+
  panel_border()

p.phylum+
  ggtitle("phage sigma factor type by host phylum")+
  ggsave2(here("vogdb","figures","tigr_Nsigma_HostPhylum.png"),
          width = 8,height = 8)


#########
# is there a pattern in pfam based group? no
p.phylum <-
  d.faa%>%
  # arrange by frequency of sigma factors
  mutate(sigma.group=fct_infreq(sigma.group))%>%
  group_by(phylum,sigma.group)%>%
  summarise(n=n())%>%
  ungroup()%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))%>%
  ggplot( aes(fct_rev(sigma.group),n, group = phylum)) + 
  #line for 0
  geom_hline(yintercept = 0, color="grey")+
  geom_col(aes(fill=phylum)) + 
  ylab("Phage sigma factor genes") +
  xlab("Sigma Factor group")+
  # facet_grid(~phylum)+
  theme_cowplot()+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow = 2))+
  coord_flip()+
  panel_border()

p.phylum+
  ggtitle("phage sigma factor group by host phylum")+
  ggsave2(here("vogdb","figures","pfamGroup_sigma_HostPhylum.png"),
          width = 10,height = 10)


###################################
# Focus on phages with 3 sigma factors
d.nsig3 <- 
d.faa%>%
  filter(n.sig>0)%>%
  filter(phylum=="Firmicutes")%>%
  # arrange by frequency of sigma factors
  mutate(tigr.hit=fct_infreq(tigr.hit))
  
tigr.spore=c("spore_sigmaK","spore_sigmaE","spore_sigF","spore_sigG")

d.nsig3%>%
  # assign a positional index to each sigma factor
  group_by(taxon)%>%
  arrange(tigr.hit)%>%
  mutate(sig.position=row_number())%>%
  ungroup()%>%
  mutate(`virus name`=str_replace(`virus name`, "phage", "\u03c6"))%>%
  mutate(`virus name`=str_replace(`virus name`, "virus", "\u03c6"))%>%
  mutate(`virus name`=as_factor(`virus name`))%>%
  # Add an asterisk for sporulation factor
  mutate(spore.t=ifelse(tigr.hit%in%tigr.spore,"\u00B7",""))%>%
  #plot
  ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
  geom_tile(aes(fill=tigr.hit), color="black")+
  geom_text(aes(label=spore.t), size=15)+
  theme_cowplot()+
  scale_fill_viridis_d()+
  facet_wrap(~n.sig, scales = "free")+
  theme(axis.text.x = element_blank())+
  ggsave2(here("vogdb/figures/","sigma_TIGR_content_firmicutes.png"),
          width = 13, height = 18)

###############
d.nsig3%>%
  # assign a positional index to each sigma factor
  group_by(taxon)%>%
  arrange(sigma.group)%>%
  mutate(sig.position=row_number())%>%
  ungroup()%>%
  mutate(`virus name`=str_replace(`virus name`, "phage", "\u03c6"))%>%
  mutate(`virus name`=str_replace(`virus name`, "virus", "\u03c6"))%>%
  mutate(`virus name`=as_factor(`virus name`))%>%
  #plot
  ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
  geom_tile(aes(fill=sigma.group), color="black")+
  theme_cowplot()+
  scale_fill_viridis_d()+
  facet_wrap(~n.sig, scales = "free")+
  theme(axis.text.x = element_blank())+
  ggsave2(here("vogdb/figures/","sigma_pfam_content_firmicutes.png"),
          width = 13, height = 18)

