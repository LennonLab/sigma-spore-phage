library(here)
library(tidyverse)
library(cowplot)
library(scales)
library(foreach)
source(here("TIGR/code/parse_hmmer_tbl.R"))


# parse hmmserch results -----------------------------------------

l.res <- list.files(here("TIGR/data/hscan_vogXtigr/"), full.names = T)
# variable to store reults
hits.all <- tibble()


foreach (i = 1:length(l.res)) %do%
  {
    
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
hits.all <-
  hits.all%>%
  separate(description,into = c("id","description"),sep = ":" )




# get 3 best hits for each protein ---------------------
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
  xlab("best specific TIGRFAM match")  
  # ggsave2(filename = here("vogdb/figures/vogXtigr.png"),
  #         width = 7,height = 7)


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


# all have an e-value lower than 1e-3 expect for a few sigH hits
# I don't think there is a problem





# combining vogdb and virus-host db info on proteins ------

# load in the data made in previous script
load(here("vogdb","data","vog_sigma_clean_Whost.RData"))

#combining of data frames will be by protein number.
# extract it from hscan query_name
# remove leading taxon number from query name if present
#leaving protein id that starts with letter

best3 <- best3%>%
  mutate(protein=str_remove(best3$query_name,".faa") %>% 
         str_remove(regex("^[0-9]*\\.?")))

##QC
# anyDuplicated(best3$protein) # 0
# is.na(best3$protein)%>%sum() # 0
# all(d.faa$protein %in% best3$protein) #TRUE
# all(best3$protein %in% d.faa$protein) #TRUE


d.faa <- 
  best3%>%
  select(-query_name)%>%
  full_join(d.faa, ., by="protein")


# crASS phages ------------------------------------------------------------

#assign crASS phages to Bacteroidetes
# based on Koonin and Yutin 2020
# https://doi.org/10.1016/j.tim.2020.01.010
d.faa <- d.faa %>% 
  mutate(phylum = 
           if_else(str_detect(sp, regex("crAssphage")),
                   "Bacteroidetes/Chlorobi group", phylum))

# _________________----------------------------------------------
# parsing sigma factor type by number of sigma factors in phage --------


# add column of sigma factor content per genome
n.sig <- d.faa%>%
  group_by(taxon, `virus name`)%>%
  summarise(n.sig=n())%>%
  ungroup()

d.faa <- n.sig%>%
  select(-`virus name`)%>%
  right_join(.,d.faa, by="taxon")

# classes of sigma factors -------------------------------------

sig.class <- read_csv(here("TIGR/data/tigr_Bs_sigma.csv"))

d.faa <- sig.class%>%
  select(ID,sig.class, sig.group)%>%
  left_join(d.faa,., by=c("tigr.hit"="ID"))

d.faa$sig.class[d.faa$tigr.hit=="no_hit"] <- "no_hit"
d.faa$sig.group[d.faa$tigr.hit=="no_hit"] <- "no_hit"



# _________________----------------------------------------------
# Plot by phylum ---------------------------------------------

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


# summarize by spore function
  d.plot <- 
    d.faa%>%
    #classify by spore function
    mutate(sigma = case_when(tigr.hit %in% tigr.spore ~ tigr.hit,
                             tigr.hit == "no_hit" ~ "no hit",
                             TRUE ~ "other")) %>% 
    mutate(sigma = str_remove(sigma, "spore_")) %>% 
    mutate(sigma = str_replace(sigma, "sigma", "sig")) %>% 
    mutate(sigma = str_replace(sigma, "Sig", "sig")) %>% 
    mutate(sigma = fct_relevel(sigma, "sigF", "sigG", "sigK", "sigE", "other", "no hit") %>% 
             fct_rev()) %>%
    group_by(phylum,sigma)%>%
    summarise(n=n(), .groups = "drop") %>% 
    group_by(phylum)%>%
    mutate(perc=n/sum(n)) %>% 
    
    mutate(phylum=str_remove(phylum,"/.*"))%>%
    mutate(phylum=fct_reorder(phylum,n)) %>% 
    left_join(., n.gene.phylum, by = "phylum")
  
  p.phylum <- d.plot %>% 
    ggplot(aes(x=phyl.lab, y = perc,fill = sigma)) + 
    geom_bar(position="fill", stat="identity", color = "transparent", size = 0,width = 0.5) +
    
    xlab("Host Phylum") +
    ylab(NULL)+
    guides(fill = guide_legend("Sigma\nFactor\nType", reverse = T))+
    
    scale_y_continuous(labels=scales::percent, position = "left") +
    scale_fill_viridis_d(direction = -1) + 
    theme_classic(base_size = 8)+
    panel_border(color = "black")+
    coord_flip()
  
  
    ggsave2(here("TIGR","figures","tigr_Nsigma_HostPhylum.png"),
            p.phylum+theme(legend.position = "none"),
            width = 4,height = 1.7)   
    
    # export  to ppt ----------------------------------------------------------
    
    #export to pptx using officer and rvg
    library (officer)
    library(rvg)
    
    read_pptx() %>%
      add_slide(layout = "Blank", master = "Office Theme" ) %>%
      ph_with(dml(ggobj = p.phylum+theme(legend.position = "none")),
              location = ph_location(type = "body",
                                     left = 0, top = 0, width = 3.2, height = 1.6)) %>%
      print(target = here("TIGR","figures","tigr_Nsigma_HostPhylum.pptx"))
    



# # _________________----------------------------------------------
# # Focus on phages of Firmicutes --------------------------
# d.nsig <- 
# d.faa%>%
#   #classify by spore function
#   mutate(sigma = case_when(tigr.hit %in% tigr.spore ~ tigr.hit,
#                            tigr.hit == "no_hit" ~ "no hit",
#                            TRUE ~ "other")) %>%
#   mutate(sigma = str_remove(sigma, "spore_")) %>% 
#   mutate(sigma = str_replace(sigma, "sigma", "sig")) %>% 
#   mutate(sigma = str_replace(sigma, "Sig", "sig")) %>% 
#   mutate(sigma = fct_relevel(sigma, "sigF", "sigG", "sigK", "sigE", "other", "no hit") %>% 
#            fct_rev()) %>% 
#   # assign a positional index to each sigma factor
#   group_by(taxon)%>%
#   arrange(tigr.hit)%>%
#   mutate(sig.position=row_number())%>%
#   ungroup()%>%
#   mutate(`virus name`=str_remove(`virus name`, ".*phage "))%>%
#   mutate(`virus name`=str_remove(`virus name`, ".*virus"))%>%
#   mutate(`virus name`=as_factor(`virus name`)) %>% 
#       
#   separate(family.etc, into = c("family","genus","species","strain"),sep=";") %>% 
#   # genus is also in viral sp name as first word
#   mutate(genus2=str_extract(sp,regex(".*? ")) %>% trimws()) %>% 
#   mutate(genus = trimws(genus)) %>% 
#   mutate(genus.plot = if_else(is.na(genus), genus2, genus))
#   
#     
# 
# # Plot ---------------------------------------------------------------
#     # extract viral family
#     # d.sp <- 
#     d.nsig <- d.nsig %>% 
#       mutate(viral.family=str_extract(`virus lineage`,
#                                       regex("caudovirales;.*", ignore_case = T)))%>%
#       mutate(viral.family=str_remove(viral.family,
#                                      regex("caudovirales; ", ignore_case = T)))%>%
#       mutate(viral.family=str_remove(viral.family,
#                                      regex(";.*", ignore_case = T)))%>%
#       mutate(viral.family = if_else(is.na(viral.family), "non-Caudovirales*", viral.family)) %>% 
#       mutate(viral.family = fct_relevel(viral.family, "non-Caudovirales*", after = 0)) %>% 
#       filter(phylum == "Firmicutes")
#     
#     
#     lab_text_size=3
#     
#     l.plot <- list() 
# 
#     for (vfam in unique(d.nsig$viral.family)){    
# l.plot[[vfam]] <- 
#       d.nsig %>% 
#       mutate(genus.plot = fct_infreq(genus.plot)) %>% 
#       filter(viral.family ==vfam) %>% 
#       
#       # filter(! genus.plot %in% c("Bacillus", "Priestia")) %>%
#       # filter(! genus.plot %in% c("Staphylococcus", "Enterococcus")) %>%
#       mutate(sigma = fct_relevel(sigma, "sigF", "sigG", "sigK", "sigE", "other", "no hit") %>% 
#                fct_rev()) %>% 
#       
#       # arrange phage by n.sigmas
#       mutate(`virus name` = fct_reorder(`virus name`, n.sig)) %>% 
#       ggplot(aes(x=sig.position, y=fct_rev(`virus name`)))+
#       geom_tile(aes(fill=sigma), color="black")+
#       theme_classic()+
#       panel_border(color = "black", size=0.1)+
#       scale_fill_viridis_d(direction = -1, drop = FALSE)+
#       facet_grid(genus.plot~., scales = "free", space = "free")+
#       geom_text(aes(label = genus.plot), size=lab_text_size,
#                 x = Inf, y = Inf, hjust = 1.1, vjust = 1.05) +
#       theme(strip.background = element_blank(),
#             strip.text = element_blank(),
#             legend.position = "none",
#             axis.title.y = element_blank(),
#             panel.spacing = unit(0, "lines"),
#             plot.title = element_text(hjust = 0.5))+
#       xlab("Sigma Factor Gene")+
#       expand_limits(x=5)+
#   scale_x_continuous(breaks = c(1:3))+
#   ggtitle(vfam)
#     }
# 
# #count viral families for proportional plotting    
# n.fams <- d.nsig %>% group_by(viral.family) %>% summarise(n=n())
# 
# pA <- egg::ggarrange(l.plot$Myoviridae+
#                        theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
#                # ggplot()+ theme_void(),
#                l.plot$Podoviridae+
#                  theme(axis.title.x = element_blank(),axis.text.x = element_blank()),
#                # ggplot()+ theme_void(),
#                l.plot$Siphoviridae,
#           ncol = 1,heights = c(7,6,61))
# 
# 
#   
# ggsave2(here("TIGR/figures/","sigma_TIGR_content_Firmi.png"),
#         plot = plot_grid(l.plot$Herelleviridae, pA , ncol = 2),
#         width = 8, height = 11)


# Save data ---------------------------------------------------------------
  
  save(d.faa, file = here("TIGR/data/vog_sigma_clean_tigr.RData"))

  write_csv(d.faa %>% select(-seq), file = here("TIGR/data/vog_sigma_clean_tigr.csv"))




