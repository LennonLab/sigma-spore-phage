library(here)
library(tidyverse)
library(cowplot)
library(scales)
library(seqinr)

#will use hmmalign installed in wsl

if (! dir.exists(here("vogdb/data/hmm_align"))){
  dir.create(here("vogdb/data/hmm_align"))
}

# load and select phage sequences ----------------------------------------------
# load in the data made in previous script
load(here("vogdb","data","vog_sigma_clean_tigr.RData"))

seq_select <- d.faa %>% 
  filter(sig.class == "sporulation") %>% 
  filter(str_detect(sp, regex("Eldridge", ignore_case = T))|
           str_detect(sp, regex("SP-10", ignore_case = T))|
           str_detect(sp, regex("Goe3", ignore_case = T)))


write.fasta(sequences = seq_select$seq,
            names = seq_select$protein, 
            file.out = here("vogdb/data/hmm_align","viral_sporulation_sigmas.faa"),
            open = "w")  



# > get sporulation sigma profile -------------------------------------------

tigr.sigma <- read_csv( here("TIGR/data","tigr_info_sigma.csv"))

# accession of HMM profile
sigf.AC <- tigr.sigma %>% 
  filter(str_detect(ID, regex("sigf", ignore_case = T))) %>% 
  pull(AC)

#path to HMM profile
sigf.hmm <- 
  str_subset(list.files(here("TIGR/data","tigr_hmm"), full.names = T),
                       sigf.AC) 
# copy hmm profile to wd
file.copy(sigf.hmm, here("vogdb/data/hmm_align"))

sigf.hmm <- 
  str_subset(list.files(here("TIGR/data","tigr_hmm")),
             sigf.AC) 

# > align to hmm ----------------------------------------------  
setwd(here("vogdb/data/hmm_align"))

# full alignment

wsl <- paste(" wsl cat ", 
             "viral_sporulation_sigmas.faa",
             "bact_sigF.faa |", # Bacillus and Clostr. sigF
             "hmmalign --amino --informat FASTA",
             "-o SporeSigma_X_sigfProfile_full.aln",
             sigf.hmm, "-")
system(wsl)


wsl <- paste(" wsl cat ", 
             "viral_sporulation_sigmas.faa",
             "bact_sigF.faa |", # Bacillus and Clostr. sigF
             "hmmalign --amino --informat FASTA --trim",
             "-o SporeSigma_X_sigfProfile_trim.aln",
             sigf.hmm, "-")
system(wsl)


# pid
pid_app <- "~/hmmer-3.3.1/easel/miniapps/esl-alipid"

wsl <- paste(" wsl",
             pid_app,
             "--amino",
             "SporeSigma_X_sigfProfile_trim.aln > pid_sigF.tsv")

system(wsl)


# > parse results -----------------------------------------------------------


pid_sigF <- read.table(here("vogdb/data/hmm_align", "pid_sigF.tsv")) %>% 
  as_tibble()
colnames(pid_sigF) <- c("seqname1","seqname2","perc_id","nid","denomid","perc_match","nmatch","denommatch")

# names of proteins
protein.names <- read_csv(here("vogdb/data/hmm_align","seq_names.csv"))

# replace names
pid_sigF <- left_join(pid_sigF, protein.names, by = c("seqname1" = "protein")) %>% 
  mutate(seqname1 = name) %>% 
  rename(phage1 = phage) %>% 
  select(-name) %>% 
  left_join(., protein.names, by = c("seqname2" = "protein")) %>% 
  mutate(seqname2 = name) %>% 
  rename(phage2 = phage) %>% 
  select(-name)


#organize data: phage vs host

# keep only phage-host comparisons
pid_sigF <- pid_sigF %>% 
  filter(phage1 != phage2)



# Repeat with sigG --------------------------------------------------------


# > get sporulation sigma profile -------------------------------------------

# accession of HMM profile
sigg.AC <- tigr.sigma %>% 
  filter(str_detect(ID, regex("sigg", ignore_case = T))) %>% 
  pull(AC)

#path to HMM profile
sigg.hmm <- 
  str_subset(list.files(here("TIGR/data","tigr_hmm"), full.names = T),
             sigg.AC) 
# copy hmm profile to wd
file.copy(sigg.hmm, here("vogdb/data/hmm_align"))

sigg.hmm <- 
  str_subset(list.files(here("TIGR/data","tigr_hmm")),
             sigg.AC) 

# > align to hmm ----------------------------------------------  
setwd(here("vogdb/data/hmm_align"))

# full alignment

wsl <- paste(" wsl cat ", 
             "viral_sporulation_sigmas.faa",
             "bact_sigG.faa |", # Bacillus and Clostr. sigG
             "hmmalign --amino --informat FASTA",
             "-o SporeSigma_X_sigGProfile_full.aln",
             sigg.hmm, "-")
system(wsl)


wsl <- paste(" wsl cat ", 
             "viral_sporulation_sigmas.faa",
             "bact_sigG.faa |", # Bacillus and Clostr. sigG
             "hmmalign --amino --informat FASTA --trim",
             "-o SporeSigma_X_sigGProfile_trim.aln",
             sigg.hmm, "-")
system(wsl)


# pid
pid_app <- "~/hmmer-3.3.1/easel/miniapps/esl-alipid"

wsl <- paste(" wsl",
             pid_app,
             "--amino",
             "SporeSigma_X_sigGProfile_trim.aln > pid_sigG.tsv")

system(wsl)


# > parse results -----------------------------------------------------------


pid_sigG <- read.table(here("vogdb/data/hmm_align", "pid_sigG.tsv")) %>% 
  as_tibble()
colnames(pid_sigG) <- c("seqname1","seqname2","perc_id","nid","denomid","perc_match","nmatch","denommatch")

# names of proteins
protein.names <- read_csv(here("vogdb/data/hmm_align","seq_names.csv"))

# replace names
pid_sigG <- left_join(pid_sigG, protein.names, by = c("seqname1" = "protein")) %>% 
  mutate(seqname1 = name) %>% 
  rename(phage1 = phage) %>% 
  select(-name) %>% 
  left_join(., protein.names, by = c("seqname2" = "protein")) %>% 
  mutate(seqname2 = name) %>% 
  rename(phage2 = phage) %>% 
  select(-name)


#organize data: phage vs host

# keep only phage-host comparisons
pid_sigG <- pid_sigG %>% 
  filter(phage1 != phage2)


# plot --------------------------------------------------------------------


pid <- bind_rows(pid_sigF %>% mutate(host_gene = "sigF"),
                 pid_sigG %>% mutate(host_gene = "sigG"))


p <- pid%>% 
  mutate(seqname1 = str_remove(seqname1, "_.*")) %>% 
  ggplot(aes(seqname1, seqname2))+
  geom_tile(aes(fill = perc_id))+
  geom_text(aes(label = perc_id))+
  scale_fill_viridis_c()+
  facet_wrap(~host_gene, scales = "free_x")+
  ylab("Phage")+
  xlab("Host")+
  theme_classic()+
  panel_border(color = "black")

save_plot(here("vogdb/figures","percID_align2profile.png"),p)

