library(here)
library(tidyverse)
library(cowplot)

# Collecting DE data for all loci in all treatments
# reading in the FDR corrected p-values for indicating significantly DEx'ed genes.
f <- 
list.files(here("RNAseq","data","DEseq"), pattern = "txt",recursive = T, full.names = T)
d <- read_delim(f[1], delim = "\t")
d <-   select(d,id)
d2 <- d
for (file in f){

  de <- read_delim(file, delim = "\t")
  d <- full_join(d, select(de, id, contains("IHW")))
  d2<- full_join(d2, select(de, id, contains("Fold_Change")))
}

d <- d%>%
      rename_at(vars(-id), ~str_remove(.,"_.*"))
p.val <- d
d2 <- d2%>%
      rename_at(vars(-id), ~str_remove(.,"_.*"))
#fold change
fc <- d2
#clean up
rm(de,d,d2)


#------------------
# locus tags
# list matching delta6 locus_tags with 168:  

d6.spor.genes <- read_csv(here("RNAseq/data/annotations/delta6_Spor_Annot.csv"))

#removing RNA genes
d6.spor.cds <- 
  d6.spor.genes %>%
  filter(str_detect(locus_tag.168,"BSU"))%>%
  filter(!str_detect(locus_tag.168,"RNA"))

#------------------
# The Hypergeometric Distribution   
# Using the example given by R stats package:
# This distribution decribes the probability of sampling, without replacment, 
# of X white balls from a (finite) urn containing M white and K black balls. 
# Sample size is N.

# The question we pose here is:
  # Are DExed genes when inducing the host sigF and sigG 
  # enriched in DExed genes seen when phage genes were induced?

sig.host.v <- c("sigF","sigG")
sig.phage.v <- c("ELDg168","ELDg169","Goe3","SP10","sigF","sigG")

# record p values
p.val.hostBG.cds <- tibble( host.gene=character(),phage.gene=character(),hg=numeric(), fisher=numeric(), chi=numeric())

l.plot <- list()
i=1

for (cur.sig.host in sig.host.v){ 
  for (cur.sig.phage in sig.phage.v){
    
    if(cur.sig.host == cur.sig.phage) next
    
    #prepare data to make contingency table
    urn <- 
      p.val%>%
      # #filter only cds
      # filter(str_detect(id, "A8O17"))%>%
      right_join(., d6.spor.cds, by=c("id"="locus_tag.d6"))%>%
      # select the data for the current genes
      select(sig.host=cur.sig.host,sig.phage=cur.sig.phage)%>%
      # # remove genes that were not assigned a p-value for DE
      # filter(!is.na(sig.host))%>%
      # filter(!is.na(sig.phage))%>%
      # change genes that were not assigned a p-value for DE to 1
      mutate(sig.host=if_else(is.na(sig.host),1,sig.host))%>%
      mutate(sig.phage=if_else(is.na(sig.phage),1,sig.phage))%>%
      # define logical vector of DE based on p-value
      mutate(host.sigma=sig.host<0.05, phage.sigma=sig.phage<0.05)%>%
      # change logical vecotrs to meaningful strings
      mutate(host.sigma=ifelse(host.sigma, "h.DExed", "h.unchaged"),
             phage.sigma=ifelse(phage.sigma, "p.DExed", "p.unchaged"))%>%
      # select only DE and sporulation for summary
      select(host.sigma,phage.sigma)#%>%
    # rename(!!cur.sig.host := host.sigma)%>%
    # rename(!!cur.sig.phage := phage.sigma)
    
    
    #make contingency table 
    urn <- table(urn)#to control order   
    
    # DExed genes observed in both sigmas
    x <- urn["h.DExed","p.DExed"]
    
    # DExed gene observed in host sigma
    m<- sum(urn["h.DExed",])
    
    # non-DExed genes observed in host sigma
    n <- sum(urn["h.unchaged",])
    
    #number of DE genes in phage
    k <- sum(urn[,"p.DExed"])
    
    z <- 0:min(k,m)
    
    chi <- capture.output(summary(urn))
    fisher <- signif(fisher.test(urn)$p.value,4)
    hg <- signif(phyper(x-1, m, n, k, lower.tail = F),5)
    
    #save p-values
  p.val.hostBG.cds <- 
      tibble( host.gene=cur.sig.host,
              phage.gene=cur.sig.phage,
              host.DEG=m,
              phage.DEG=k,
              shared.DEG=x,
              total=m+n,
              hg=hg,
              fisher=fisher,
              chi=str_extract(chi[4], "value = .*")%>%parse_number())%>%
      bind_rows(p.val.hostBG.cds,.)
  
    
  }
}


# parmeters to plot PDF
p.val.hostBG.cds <- 
  p.val.hostBG.cds %>% 
    mutate(M = host.DEG,
           N = total-host.DEG,
           K = phage.DEG,
           X = pmin(K,M),
           who = paste(host.gene, phage.gene, sep = "_")) 

# calculate PDF per interaction
max.x <- max(p.val.hostBG.cds$X)
d.hyp <- tibble(x= seq(max.x))

for(i in 1:nrow(p.val.hostBG.cds)){
  
d.hyp[, p.val.hostBG.cds$who[i]] <- 
  with(p.val.hostBG.cds,
  dhyper(x = seq(max.x), m = M[i], n = N[i], k = K[i]))
}

# plot
p <- d.hyp %>% 
  pivot_longer(-x) %>% 
  separate(name, into = c("host.gene", "phage.gene"), sep = "_", remove = F) %>% 
  filter(value > 1e-10) %>% 
  ggplot(aes(x, value))+
  geom_area(fill = "grey70", color = "grey30")+
  geom_vline(data = p.val.hostBG.cds, aes(xintercept = shared.DEG),
             color = "red", size = 1)+
  facet_wrap(host.gene ~ phage.gene, scales = "free", nrow = 2)+
  theme_classic()+
  panel_border(color = "black")+
  scale_y_continuous(expand = c(0, 0))+
  xlab("shared DExed genes")+
  ylab("PDF")

ggsave(here("RNAseq/plots/enrichment1.png"),plot = p, width = 8, height = 6)

#---------------------
#   Are sporulation genes enriched in the sample of deferentially expressed genes?
# the FDR corrected p-values for indicating significantly DExd genes.