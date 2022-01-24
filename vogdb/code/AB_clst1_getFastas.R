library(here)
library(tidyverse)

# Reviewer comment:
# for all these phages the sequence heterogeneity is substantial,
# so the data might be clearer if there was subdivision into 
# groups of sequence-related phage genomes.


# Based on response from Wrighton Lab (CSU) members will cluster with the "ClusterGenomes" tool 
# on the Cyverse discovery environement (https://de.cyverse.org/apps/agave/ClusterGenomes-1.1.3u3)
# Described here: https://dx.doi.org/10.17504/protocols.io.gwebxbe
# HOWEVER 
# that method did not give back results :(
# As an alternative I will preform clustering using cd-hit-est 
# (https://github.com/weizhongli/cdhit/wiki)
# This has been done in the MVP (Microbe Versus Phage) database
# http://mvp.medgenius.info/help
# altough they used the parameter recomended below, I could not match the phages in MVP to VOG
# So instead I am going ro preform clustering with CD-HIT on Carbonate.

## Parameters
# According to the "Minimum Information about an Uncultivated Virus Genome (MIUViG)" 
# (https://doi.org/10.1038/nbt.4306)
    # We propose to formalize the use of species-rank virus groups and to name these 
    # “virus operational taxonomic units” (vOTUs) to avoid confusion because species 
    # groups have been variously named “viral population,” “viral cluster” or “contig 
    # cluster” in the literature4,7,60. We suggest standard thresholds of 95% average
    # nucleotide identity over 85% alignment fraction (relative to the shorter sequence)


# The input needed for clustering is the genomic sequences of all if the phages in the analysis.
# Those are all the phages used to assemble the VOG data base.


# get VOG phage Accessions ------------------------------------------------


# import all viruses used to assemble VOGs (from A_get_vog_files.R).
vog.sp <-  read_tsv(here("vogdb/data","vogdb_downloads","vog.species.list"), trim_ws = T) %>% 
  as_tibble( .name_repair = "universal") %>% 
  # keep only phages
  filter(phage.nonphage == "phage")

# use virus host database to find genome accesion
# download.file(url = "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv",
#               destfile = here("vogdb/data", "virushostdb.tsv"))
# downloaded on (6/JAN/2022)
vh.db <- read_tsv(here("vogdb/data","virushostdb.tsv"))

# All in v-h db?
sum(!vog.sp$tax.id %in% vh.db$`virus tax id`) #0
# yes!

#keep only relevant phages
vh.db <- vh.db %>% 
  filter(vh.db$`virus tax id` %in% vog.sp$tax.id)
  
# Duplicates in v-h db
vh.db.dups <- vh.db %>% 
  group_by(`virus tax id`) %>% 
  summarise(n = n()) %>% 
  filter(n>1)

# do any of the duplicates have different accessions
# or can I just get rid of duplicates?

vh.db %>% 
  filter(`virus tax id` %in% vh.db.dups$`virus tax id`) %>% 
  group_by(`virus tax id`, `refseq id`) %>% 
  summarise(n = n()) %>% 
  ungroup() %>% 
  group_by(`virus tax id`) %>% 
  summarise(n = n()) %>% 
  filter(n>1)
  
# There is a single phage (taxid = 31754)
# Lactococcus phage BK5-T
# in one row it has only "NC_002796" and in another it has both NC_002796, NC_056724.
# I'll get rid of the row with two accessions
rmv <- which(
  (vh.db$`virus tax id`==31754) & 
  (str_detect(vh.db$`refseq id`,","))
)

vh.db <- vh.db[-rmv,]

#remove other duplicates
vh.db <- vh.db %>% 
  filter(!duplicated(`virus tax id`))


# add ref acc to vog
vog.sp <- 
  vh.db %>% 
  select(tax.id = `virus tax id`, acc = `refseq id`) %>% 
  left_join(vog.sp, . , by = "tax.id")

# clean data
# there are a few phages that have multiple refseq_ids (separated by comma)
vog.sp$acc %>% str_count(",") %>% table(.)
#    0    1    2    4 
# 4509    5    8    1 
#14 in total 
multi_ref <- 
  vog.sp %>% filter(str_count(acc,",") > 0)

# most are segmented phages.
    #segmentsed phages
    # Pseudomonas phage phi2954
    # Pseudomonas phage phiYY
    # Pseudomonas phage phi12
    # Pseudomonas virus phi6
    # Pseudomonas phage phi13
    # Pseudomonas phage phiNN
    # Campylobacter virus IBB35 (contigs)
    # Pseudomonas phage phi8 (contigs?)
    # Salmonella phage SP069 (contigs?)


# this is a small number so I will just filter them out.
vog.sp <- vog.sp %>%
  filter(!str_count(acc,",") > 0)

#save data
write_csv(vog.sp, here("vogdb/data", "vog-members.csv"))

# Get fastas for clustering -----------------------------------------------

# import virus-host database genome sequences
# download.file(url = "ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.genomic.fna.gz",
#               destfile = here("vogdb/data", "virushostdb.genomic.fna.gz"))

# index the genome data
d.vhdb_fna <- tibble(
  # number of the header line
  header_line = readLines(here("vogdb/data", "virushostdb.genomic.fna.gz")) %>% 
    grep(">",x = . ,value = F),
  # the header
  header = readLines(here("vogdb/data", "virushostdb.genomic.fna.gz")) %>% 
    grep(">",x = . ,value = T)
)

# Parse header 
d.vhdb_fna <- d.vhdb_fna %>% 
  mutate(header = str_remove(header,">")) %>% 
  separate(header, into = c("locus","gb_acc","sp"), sep = " \\[|\\] ", 
           remove = F, extra = "merge")

# end line for each fasta
last_line <- system (paste("wc -l", here("vogdb/data", "virushostdb.genomic.fna.gz")),
                     intern = T) %>% parse_number()
d.vhdb_fna$end_line <- c(d.vhdb_fna$header_line[-1]-1, last_line)


#are all the phages in there?
sum(!vog.sp$acc %in% d.vhdb_fna$gb_acc) #0
# yes!

# filter to keep only relevant fasta indices
d.vhdb_fna <- d.vhdb_fna %>% 
  filter(gb_acc %in% vog.sp$acc)

# unzip for accelrating read_lines
setwd(here("vogdb/data"))
system(paste("gunzip",here("vogdb/data", "virushostdb.genomic.fna.gz")))
setwd(here())

# write VOG fastas to new file

fa_dir <- here("vogdb/data", "vog_genomic")
if(!dir.exists(fa_dir)){
  dir.create(fa_dir)
}

# using parallel for loop to speed it up a bit.
library(foreach)
library(doParallel)
cl <- makePSOCKcluster(4)
registerDoParallel(cl)
foreach(i = 1:nrow(d.vhdb_fna),.packages = c("here","readr")) %dopar% {
  x <- read_lines(here("vogdb/data", "virushostdb.genomic.fna"),
             skip = d.vhdb_fna$header_line[i]-1,
             n_max = d.vhdb_fna$end_line[i] - d.vhdb_fna$header_line[i]+1
  ) 
  write_lines(x,here("vogdb/data", "vog_genomic", paste0(d.vhdb_fna$locus[i],".fna")),)
}
stopCluster(cl)

# merge all VOG fastas to a single multi fasta (using linux)
setwd(here())
sys.cmd <- "cat vogdb/data/vog_genomic/*.fna > vogdb/data/vog_all.fna"
system(sys.cmd, intern = T)


# cluster genomes ---------------------------------------------------------

# cd-hit was next run AB_clst2_CdHit.sh

# clean up ----------------------------------------------------------------
unlink(here("vogdb/data", "virushostdb.genomic.fna.gz"))
unlink(here("vogdb/data", "virushostdb.genomic.fna"))
unlink(here("vogdb/data/vog_genomic/"), recursive = T )
