# make contingency table 
contin.t1 <- 
  d.sp%>%
  # name adjust
  mutate(phylum =
           case_when(str_detect(phylum, "Cyanobacteria") ~ "Cyanobacteria",
                     str_detect(phylum, "Bacteroidetes") ~ "Bacteroidetes",
                     TRUE ~ phylum)) %>% 
  # remove singeltons
  group_by(phylum) %>% filter( n() > 1 ) %>% ungroup() %>% 
  mutate(has.sigma = n.sigma>0) %>% 
  select(phylum, has.sigma)%>%
  table()

m1 <- chisq.test(contin.t1, simulate.p.value = TRUE, B = 1e6)
m2 <- fisher.test(contin.t1, simulate.p.value = TRUE, B = 1e6)

write.csv(contin.t1, here("vogdb/data/","Hphyla_sigma_contingency.csv"))

# viral family -------------------------------------------------------------



# make contingency table 
contin.t2 <- 
  d.sp%>%
  #extract viral family
  mutate(viral.family=str_extract(`virus lineage`,
                                  regex("caudovirales;.*", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex("caudovirales; ", ignore_case = T)))%>%
  mutate(viral.family=str_remove(viral.family,
                                 regex(";.*", ignore_case = T))) %>%
  mutate(viral.family = if_else(is.na(viral.family), "non-Caudovirales", viral.family)) %>% 
  # remove singeltons
  group_by(viral.family) %>% filter( n() > 1 ) %>% ungroup() %>%

  mutate(has.sigma = n.sigma>0) %>% 
  select(viral.family, has.sigma)%>%
  table()

m1 <- chisq.test(contin.t2, simulate.p.value = TRUE, B = 1e6)
m2 <- fisher.test(contin.t2, simulate.p.value = TRUE, B = 1e6)

write.csv(contin.t2, here("vogdb/data/","VIRfam_sigma_contingency.csv"))
