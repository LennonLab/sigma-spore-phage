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


# get sigBFG sigma profile -------------------------------------------

tigr.sigma <- read_csv( here("TIGR/data","tigr_info_sigma.csv"))

# accession of HMM profile
sig_bfg.AC <- tigr.sigma %>% 
  filter(str_detect(ID, regex("sigbfg", ignore_case = T))) %>% 
  pull(AC)

#path to HMM profile
sig_bfg.hmm <- 
  str_subset(list.files(here("TIGR/data","tigr_hmm"), full.names = T),
             sig_bfg.AC) 
# copy hmm profile to wd
file.copy(from = sig_bfg.hmm, to = here("vogdb/data/hmm_align"), overwrite = T)

sig_bfg.hmm <- 
  str_subset(list.files(here("TIGR/data","tigr_hmm")),
             sig_bfg.AC) 

# align to hmm ----------------------------------------------  
setwd(here("vogdb/data/hmm_align"))

# full alignment

wsl <- paste(" wsl cat ", 
             "viral_sporulation_sigmas.faa",
             "bsu_sigFG.faa|", # Bacillus subtilis. sigF + sigG
             "hmmalign --amino --informat FASTA",
             "-o SporeSigma_X_sigBFGProfile_full.aln",
             sig_bfg.hmm, "-")
system(wsl)


wsl <- paste(" wsl cat ", 
             "viral_sporulation_sigmas.faa",
             "bsu_sigFG.faa|", # Bacillus subtilis. sigF + sigG
             "hmmalign --amino --informat FASTA --trim",
             "-o SporeSigma_X_sigBFGProfile_trim.aln",
             sig_bfg.hmm, "-")
system(wsl)

# names of proteins
protein.names <- read_csv(here("vogdb/data/hmm_align","seq_names.csv"))
x <- readLines("SporeSigma_X_sigBFGProfile_full.aln")
x1 <- x
for(i in 1:nrow(protein.names)){
  x1 <- gsub(protein.names$protein[i], protein.names$name[i],x1)
}
writeLines(x1,"NAMES_SporeSigma_X_sigBFGProfile_full.aln" )

# Perccent ID ----

# pid
pid_app <- "~/hmmer-3.3.1/easel/miniapps/esl-alipid"

wsl <- paste(" wsl",
             pid_app,
             "--amino",
             "SporeSigma_X_sigBFGProfile_trim.aln > pid_sigBFG.tsv")

system(wsl)



# parse results -----------------------------------------------------------


pid_sigBFG <- read.table(here("vogdb/data/hmm_align", "pid_sigBFG.tsv")) %>% 
  as_tibble()
colnames(pid_sigBFG) <- c("seqname1","seqname2","perc_id","nid","denomid","perc_match","nmatch","denommatch")

# names of proteins
protein.names <- read_csv(here("vogdb/data/hmm_align","seq_names.csv"))

# replace names
pid_sigBFG <- left_join(pid_sigBFG, protein.names, by = c("seqname1" = "protein")) %>% 
  mutate(seqname1 = name) %>% 
  rename(phage1 = phage) %>% 
  select(-name) %>% 
  left_join(., protein.names, by = c("seqname2" = "protein")) %>% 
  mutate(seqname2 = name) %>% 
  rename(phage2 = phage) %>% 
  select(-name)


#organize data: phage vs host

# remove phage-phage comparisons
pid_sigBFG <- pid_sigBFG %>% 
  filter(! (phage2 & phage1))

# Plot --------------
p <- pid_sigBFG%>% 
  # mutate(perc_id = paste0(perc_id,"%")) %>% 
  ggplot(aes(seqname1, seqname2))+
  geom_tile(aes(fill = perc_id), color = "white")+
  geom_text(aes(label = paste0(perc_id,"%")))+
  scale_fill_viridis_c()+
  ylab(NULL)+
  xlab(NULL)+
  theme_classic()+
  panel_border(color = "black")

save_plot(here("vogdb/figures","percID_align2profile.png"),p,
          base_height = 2,base_width = 4)


