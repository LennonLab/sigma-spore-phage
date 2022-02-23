library(here)
library(tidyverse)

# Read in data --------------------------------

# parse cd-hit-est clusters
clstr <- read_lines(here("vogdb/data/cd-clusts/vog_cd-clust_MIUViG.clstr"))
  
  # initialize
  d.clstr <- tibble()
  n.clstr = 0

for(i in 1:length(clstr)){
  
  # assign cluster number
  n.clstr <- if_else(str_detect(clstr[i],"^>Cluster",), 
                    parse_number(clstr[i]),n.clstr)
  #row with cluster number has no more data
  if(str_detect(clstr[i],"^>Cluster")) next
  
  n.member <- str_remove(clstr[i], "\t.*")
  acc <- str_remove(clstr[i],".*, >") %>% str_remove("\\.\\.\\..*")
  nt <- str_extract(clstr[i],"\t.*nt") %>% parse_number()
  is.rep <- str_detect(clstr[i], "\\*$")
  perc.id <- if_else(is.rep, 100,
                     str_extract(clstr[i], "/.*%") %>% parse_number())
  
  
  d.clstr <- tibble(n_clstr = n.clstr,
                    acc = acc,
                    nt = nt,
                    perc_id_to_rep = perc.id,
                    is_rep = is.rep) %>% 
    bind_rows(d.clstr, .)
  
}

  
# read vog members data
d.vog <- read_csv(here("vogdb/data","vog-members.csv"))

# (d.vog$acc %in% d.clstr$acc) %>% table() # all TRUE
# (d.clstr$acc %in% d.vog$acc) %>% table() # all TRUE

# add cluster data
d.vog <- left_join(d.vog, d.clstr, by = "acc")

write_csv(d.vog, here("vogdb/data", "vog-members_cd-clstr.csv"))
