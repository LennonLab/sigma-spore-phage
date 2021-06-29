library(tidyverse)
library(here)
library(cowplot)
library(ggh4x)

#load data frame of phage sigma factors generated in (B)
load(file=here("vogdb","data","vog_sigma_clean.RData"))

d.faa <- d.faa%>%
  mutate(taxon=as.numeric(taxon))


# import all viruses used to assemble VOGs
vog.sp <-  read_tsv(here("vogdb","vogdb_downloads","vog.species.list"), trim_ws = T) %>% 
  as_tibble( .name_repair = "universal")

#import virus-host data
# downloaded from ftp://ftp.genome.jp/pub/db/virushostdb/ (24/Nov/2020)
vh.db <- read_tsv(here("vogdb","vogdb_downloads","virushostdb.tsv"))

##############################################
#  virus duplicates in VHDB data
# These reflect multiple hosts
duplicated(vh.db$`virus tax id`)%>%sum() # 3461
vh.db%>%
  filter(str_detect(`host lineage`,regex("bacteria",ignore_case = T)))%>%
  group_by(`virus tax id`, `virus name`)%>%
  summarise(n=n())%>%
  ggplot(aes(x=n))+
  geom_histogram()+
  scale_y_log10()+
  ggtitle("VHDB host number for phages" )


#-----------
d.sp <- left_join(vog.sp, vh.db, by = c("tax.id" = "virus tax id"))


# I will deal with this only in the viruses present in sigma df

d.faa <- left_join(d.faa,vh.db,by=c("taxon"="virus tax id"))

d.faa %>% 
  group_by(taxon, `virus name`, protein)%>%
  summarise(n=n())%>%
  arrange(desc(n))
#TsarBomba has s digma factors and 4 hosts
# lets look at one of the proteins
d.faa%>%
  filter(protein== 'YP_009206960.1' )%>%
  select(`host name`, evidence)

# `host name`                             evidence          
# <chr>                                   <chr>             
# 1 Bacillus cereus                         Literature        
# 2 Bacillus thuringiensis                  Literature        
# 3 Bacillus thuringiensis serovar kurstaki Literature, RefSeq
# 4 Bacillus anthracis str. Sterne          Literature  

# All are in the same genus, and one has Refseq listed inder evidence. 
# The genome paper (PMID: 26472830) lists the Refseq host as its isolation host
# I will use that to select the host kept

# First I will break up the host lineage data so I verify that hosts are od simmilat taxonomy
# the host lineage is not uniform havin 4-9 levels
# str_count(d.faa$`host lineage`,";")%>%range()
# str_count(d.faa$`host lineage`,";")%>%hist()
# the last looks like this: 
# Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales
# Which corresponds to 
# Domain; Phylum; class; order. The next level down would be family

#but some have an intermediate rank between domain and phylum
# specifically FCB group AND Terrabacteria group 
# I will remove those for consistency

d.faa <- d.faa%>%
  mutate(`host lineage`=str_remove(`host lineage`," FCB group;"))%>%
  mutate(`host lineage`=str_remove(`host lineage`," Terrabacteria group;"))

d.faa <- d.faa%>%
  separate(`host lineage`, sep="; ",
           into = c("domain", "phylum", "class", "order", "family.etc"), 
           extra = "merge")


# Back to phage duplicates due to multiple hosts
dup <- d.faa%>%
  group_by(taxon, `virus name`, protein)%>%
  summarise(n.host=n())%>%
  filter(n.host>1)%>%
  ungroup()%>%
  # keep only one protein for each phage
  filter(!duplicated(taxon))
  # group_by(taxon)%>%
  # summarise(nn=n())
#there are 38 phages with duplicate hosts

look <- d.faa%>%
  filter(taxon %in% dup$taxon)%>%
  filter(protein %in% dup$protein)%>%
  select(taxon,`virus name`,`host name`,domain, phylum, class, order, family.etc, evidence)

# look%>%
#   group_by(taxon,order)%>%
#   summarize(n=n())%>%
#   filter(n<2)
# 
# look%>%
#   ggplot(aes(x=order,y=fct_rev(`virus name`)))+
#   geom_tile()+
#   theme(axis.text.x = element_text(angle=40, hjust=1))

# all hosts for each phage are the same down to order
# some phages have refseq in evidence for one of thier hosts
# I think this is the isolation host (it is for TsarBomba)

#remove host duplicate
#add host number to assist and to remember
d.faa <- dup%>%
  select(taxon,n.host)%>%
  left_join(d.faa, ., by="taxon")
#add 1 host for all the rest
d.faa$n.host[is.na(d.faa$n.host)] <- 1

d.faa <- d.faa%>%
  mutate(filt=interaction(taxon,protein))%>%
  filter(!duplicated(filt))

# # save data
# write.csv(here("vogdb","data","vog_sigma_clean_Whost.csv"))
# save(d.faa,file = here("vogdb","data","vog_sigma_clean_Whost.RData"))
# # load(file = here("vogdb","data","vog_sigma_clean_Whost.RData"))

# Which host phyla are present in the sigma factor dataset?

d.faa%>%
  group_by(phylum)%>%
  summarise(n=n())%>%
  mutate(phylum=fct_reorder(phylum,n))%>%
  ggplot(aes(phylum,n))+
  geom_col()+
  coord_flip()+
  ggtitle("phage sigma factor by host phylum")+
  theme_cowplot()

#number of sigma factors per phage
d.faa%>%
  group_by(taxon,phylum)%>%
  summarise(n=n())%>%
  mutate(phylum=fct_reorder(phylum,n))%>%
  ggplot(aes(phylum,n))+
  geom_boxplot(fill="grey")+
  geom_jitter(shape=21, alpha=0.5, width=0.3, color="blue", height = 0.1)+
  coord_flip()+
  scale_y_continuous(breaks = seq(2,12,2))+
  ggtitle("phage sigma factor content by host phylum")+
  theme_cowplot()


#How many unique viruses have sigma factors?
length(unique(d.faa$taxon))
#471

# Phylum independent summary of sigma factors/genome
d.faa%>%
  group_by(taxon)%>%
  summarise(n.sigma=n())%>%
  ungroup()%>%
  group_by(n.sigma)%>%
  summarise(n.genomes=n())%>%
  mutate(perc=100*n.genomes/sum(n.genomes))

# n.sigma n.genomes  perc
# <int>     <int> <dbl>
# 1       1       397 84.3 
# 2       2        50 10.6 
# 3       3        24 5.10

# summary of sigma factors/genome
d.faa%>%
  group_by(taxon, phylum)%>%
  summarise(n.sigma=n())%>%
  filter(n.sigma>1)%>%
  ungroup()%>%
  group_by(n.sigma, phylum)%>%
  summarise(n.genomes=n())%>%
  mutate(perc=100*n.genomes/sum(n.genomes))

##########################
# Plot by phylum
##########################
# summary table for N by phylum
n.labs <- d.faa%>%
  group_by(taxon,phylum)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(phylum)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))

#--------------
p.phylum <- 
  d.faa%>%
  group_by(taxon,phylum)%>%
  summarise(n=n())%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=fct_reorder(phylum,n))%>%
  ggplot( aes(n, group = phylum)) + 
  geom_bar(aes(y = ..prop..), stat="count") + 
  geom_text(data=n.labs, aes(label=lab), x=3,y=0.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_nested("Host Phylum" + phylum~., nest_line = TRUE)+
  theme_classic(base_size = 13)+
  panel_border(color = "black")

p.phylum+
  ggsave2(here("vogdb","figures","sigma_HostPhylum.png"),
          width = 4,height = 8)
#-------------
# p.phylum <- 
# d.faa%>%
#   group_by(taxon,phylum)%>%
#   summarise(n=n())%>%
#   mutate(phylum=str_replace(phylum,"/","\n"))%>%
#   mutate(phylum=fct_reorder(phylum,n))%>%
# ggplot( aes(n, group = phylum)) + 
#   geom_bar(aes(y = ..prop..), stat="count") + 
#   geom_text(data=n.labs, aes(label=lab), x=2,y=1)+
#   scale_y_continuous(labels=scales::percent) +
#   ylab("Phage genomes") +
#   xlab("Sigma factors per genome")+
#   facet_grid(~phylum)+
#   theme_cowplot()+
#   panel_border()
#   
# p.phylum+
#   ggtitle("phage sigma factor content by host phylum")+
#   ggsave2(here("vogdb","figures","sigma_HostPhylum.png"),
#           width = 10,height = 3.5)

# make contingency table 
contin.t <- 
d.faa%>%
  #summarise sigma factors per phage
  group_by(taxon,phylum)%>%
  summarise(n.sig=n())%>%
  ungroup()%>%
  select(-taxon)%>%
  table()

m1 <- chisq.test(contin.t)
m2 <- fisher.test(contin.t,workspace = 2e8)

##########################
# Plot by order
##########################
# summary table for N by order
n.labs <- d.faa%>%
  group_by(taxon,order)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(order)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))%>%
  mutate(order=str_replace(order,"/","\n"))%>%
  mutate(order=fct_reorder(order,n))



d.faa%>%
  filter(!duplicated(protein))%>%
  group_by(taxon,order)%>%
  summarise(n=n())%>%
  mutate(order=str_replace(order,"/","\n"))%>%
  mutate(order=fct_reorder(order,n))%>%
  ggplot( aes(n, group = order)) + 
  geom_bar(aes(y = ..prop..), stat="count") + 
  geom_text(data=n.labs, aes(label=lab), x=2,y=.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_wrap(~order)+
  theme_cowplot()+
  panel_border()+
  ggtitle("phage sigma factor content by host order")+
  ggsave2(here("vogdb","figures","sigma_HostOrder.png"),
          width = 10,height = 6)

##########################
# Plot by genus within firmicutes
##########################

# extracting genus data for  firmicutes
firmi <- d.faa%>%
  filter(str_detect(phylum,regex("firmicutes", ignore_case = T)))%>%
  # filter(!str_detect(order,regex("lacto", ignore_case = T)))%>%
  separate(family.etc, into = c("family","genus","species","strain"),sep=";")

firmi$sp[is.na(firmi$genus)]
#   "Staphylococcus phage 812" "Geobacillus virus E2" 
# genus is also in viral sp name as first word
firmi <- firmi%>%
  mutate(genus2=str_extract(sp,regex(".*? ")))

# summary table for N by order
n.labs <- firmi%>%
  group_by(taxon,genus2)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(genus2)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))

p.genus <- 
firmi%>%
  group_by(taxon,genus2)%>%
  summarise(n=n())%>%
  ggplot( aes(n, group = genus2)) + 
  geom_bar(aes(y = ..prop..), stat="count") + 
  geom_text(data=n.labs, aes(label=lab), x=2,y=.9)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_wrap(~genus2)+
  theme_cowplot()+
  panel_border()

p.genus+
  ggtitle("phage sigma factor content by host genus",
          "Firmicute infecting phages")+
  ggsave2(here("vogdb","figures","sigma_HostGenus_firmicutes.png"),
          width = 7,height = 7)



##########################
# Plot by viral type
##########################

# First I will break up the virus lineage data so
# the host lineage is not uniform havin 4-10 levels
str_count(d.faa$`virus lineage`,";")%>%range()
str_count(d.faa$`virus lineage`,";")%>%hist()
# the last looks like this: 
# Viruses; Duplodnaviria; Heunggongvirae; Uroviricota; Caudoviricetes;
#  Caudovirales; Siphoviridae; unclassified Siphoviridae

# All but 1 belong to the tailed phages (Caudovirales)
# the single exception is Salisaeta icosahedral phage 1
# indeed this is " a tailless bacteriophage" (PMID: 22509017)
  # tmp <- (!str_detect(d.faa$`virus lineage`,"Caudovirales"))%>%which()
  # d.faa$`virus name`[tmp]
# I will only look at the family level below tailed phage (order Caudovirales)
d.faa <- d.faa%>%
  mutate(viral.family=str_extract(`virus lineage`,
                                  regex("caudovirales;.*", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                  regex("caudovirales; ", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex(";.*", ignore_case = T)))

# d.faa$`virus name`[is.na(d.faa$viral.family)]
# One virus left without family assignment
# "Salisaeta icosahedral phage 1" 
# d.faa$`virus lineage`[is.na(d.faa$viral.family)]
# I will list is as unclassified(SSIP-1)
d.faa$viral.family[is.na(d.faa$viral.family)] <- "unclassified (SSIP-1)"

# make contingency table 
contin.t <- 
  d.faa%>%
  #summarise sigma factors per phage
  group_by(taxon,viral.family)%>%
  summarise(n.sig=n())%>%
  ungroup()%>%
  select(-taxon)%>%
  table()

m1 <- chisq.test(contin.t)
m2 <- fisher.test(contin.t,workspace = 6e8)



#########################
# save data
d.faa%>%
  select(-seq)%>%
write_csv(here("vogdb","data","vog_sigma_clean_Whost.csv"))
save(d.faa,file = here("vogdb","data","vog_sigma_clean_Whost.RData"))
# load(file = here("vogdb","data","vog_sigma_clean_Whost.RData"))

# summary table for N by viral.family
n.labs <- d.faa%>%
  group_by(viral.family, taxon)%>%
  summarise(n=n())%>%
  ungroup()%>%
  group_by(viral.family)%>%
  summarise(n=n())%>%
  mutate(lab=paste0("n=",n))



p.vir <-
  d.faa%>%
  group_by(taxon,viral.family, phylum)%>%
  summarise(n.sig=n())%>%
  ungroup()%>%
  group_by(viral.family, n.sig,phylum)%>%
    summarise(n.genomes=n())%>%
  ungroup()%>%
  group_by(viral.family)%>%
    mutate(perc=n.genomes/sum(n.genomes))%>%
  mutate(phylum=str_replace(phylum,"/","\n"))%>%
  mutate(phylum=str_replace(phylum,"Proteobacteria","Proteo\nbacteria"))%>%
  mutate(viral.family=fct_reorder(viral.family,n.genomes))%>%
  ggplot( aes(n.sig,perc)) + 
  geom_col(aes(fill=phylum)) + 
  geom_text(data=n.labs, aes(label=lab), x=2,y=1)+
  scale_y_continuous(labels=scales::percent) +
  ylab("Phage genomes") +
  xlab("Sigma factors per genome")+
  facet_wrap(~viral.family)+
  theme_cowplot()+
  panel_border()+
  theme(legend.position = "bottom")+
  guides(fill=guide_legend(nrow=2))
  


p.vir+
  ggtitle("phage sigma factor content by viral family")+
  ggsave2(here("vogdb","figures","sigma_ViralFamily.png"),
          width = 10,height = 6)

# Who are the siphoviruses with 3 sigma factors?
d.faa%>%
  filter(viral.family=="Siphoviridae")%>%
  group_by(`virus name`)%>%
  summarise(n.sigma=n())%>%
  filter(n.sigma>1)%>%
  arrange(desc(n.sigma))


# `virus name`                   n.sigma
# <chr>                            <int>
# 1 Bacillus phage Slash                 3
# 2 Bacillus phage Stahl                 3
# 3 Bacillus phage Staley                3
# 4 Bacillus phage Stills                3
# 5 Bacillus phage BceA1                 2
# 6 Brevibacillus phage Jenst            2
# 7 Cyanophage PSS2                      2
# 8 Staphylococcus phage SpaA1           2
# 9 Synechococcus phage S-CBS2           2

# They are all Bacillus phages(n=4)
# with 2 sigmas there are other hosts

#combine plots
bottom_row <- plot_grid(p.genus,p.vir,
                        labels = LETTERS[2:3])
plot_grid(p.phylum,bottom_row,
          rel_heights= c(1, 1.5),
          nrow = 2, labels = c('A', ''))+
  ggsave2(here("vogdb","figures","sigma_taxonomy.png"),
            width = 12,height = 10)
