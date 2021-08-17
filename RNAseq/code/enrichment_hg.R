library(here)
library(tidyverse)
library(cowplot)
library(gtools)


# Import data -------------------------------------------------------------

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



# * arrange data -------------------------------------------------


# locus tags
# list matching delta6 locus_tags with 168:  

d6.spor.genes <- read_csv(here("RNAseq/data/annotations/delta6_Spor_Annot.csv"))

#removing RNA genes
d6.spor.cds <- 
  d6.spor.genes %>%
  filter(str_detect(locus_tag.168,"BSU"))%>%
  filter(!str_detect(locus_tag.168,"RNA"))

# Hypergeometric enrichment ------------------
# The Hypergeometric Distribution   
# Using the example given by R stats package:
# This distribution decribes the probability of sampling, without replacment, 
# of X white balls from a (finite) urn containing M white and K black balls. 
# Sample size is N.

# * overlap -----------------------------------------------------------------


# The question we pose here is:
  # Are DExed genes when inducing the host sigF and sigG 
  # enriched in DExed genes seen when phage genes were induced?
  # The match we consider accounts for the direction of change being the same.

sig.host.v <- c("sigF","sigG")
sig.phage.v <- c("ELDg168","ELDg169","Goe3","SP10","sigF","sigG")

# record p values
p.val.hostBG.cds <- tibble( )
#host.gene=character(),phage.gene=character(),hg=numeric(), fisher=numeric(), chi=numeric())


for (cur.sig.host in sig.host.v){ 
  for (cur.sig.phage in sig.phage.v){
    
    if(cur.sig.host == cur.sig.phage) next
    
    #prepare data to make contingency table
    d.urn <- 
      p.val%>%
      # #filter only cds
      right_join(., d6.spor.cds, by=c("id"="locus_tag.d6"))%>%
      # select the data for the current genes
      select(id, p.sig.host=all_of(cur.sig.host),p.sig.phage=all_of(cur.sig.phage))
      
      #add data on direction of change
    d.urn <- fc%>%
      # #filter only cds
      right_join(., d6.spor.cds, by=c("id"="locus_tag.d6"))%>%
      # select the data for the current genes
      select(id, fc.sig.host=cur.sig.host, fc.sig.phage=cur.sig.phage) %>%
      left_join(d.urn, ., by = "id")
    
      # mark significance of change
      d.urn <- d.urn %>% 
        # change genes that were not assigned a p-value for DE to 1
        mutate(p.sig.host=if_else(is.na(p.sig.host),1,p.sig.host))%>%
        mutate(p.sig.phage=if_else(is.na(p.sig.phage),1,p.sig.phage))%>%
        # define logical vector of DE based on p-value
        mutate(host.DExed=p.sig.host<0.05, phage.DExed=p.sig.phage<0.05)
      
      # check direction of change is the same
      d.urn <- d.urn %>% 
      mutate( host.change = case_when(!host.DExed ~ "h.unchaged",
                                      fc.sig.host > 2 ~ "h.up",
                                      fc.sig.host < 0.5 ~ "h.down",
                                      TRUE ~ "h.unchaged")) %>% 
      mutate( phage.change = case_when(!phage.DExed ~ "p.unchaged",
                                      fc.sig.phage > 1 ~ "p.up",
                                      fc.sig.phage < 0.5 ~ "p.down",
                                      TRUE ~ "p.unchaged"))

 
    #make contingency table 
    urn <-   d.urn %>% 
        select(host.change, phage.change) %>% 
        table()

    
    # DExed genes observed in both sigmas
    x <- urn["h.down", "p.down"] + urn["h.up", "p.up"]
    
    # DExed gene observed in host sigma
    m<- sum(urn["h.up",]) +  sum(urn["h.down",])
    
    # non-DExed genes observed in host sigma
    n <- sum(urn["h.unchaged",])
    
    #number of DE genes in phage
    k <- sum(urn[,"p.up"]) +  sum(urn[, "p.down"])
    
    z <- 0:min(k,m)
    
    # chi <- capture.output(summary(urn))
    hg.left <- signif(phyper(x-1, m, n, k, lower.tail =T),5) 
    hg.right <- signif(phyper(x-1, m, n, k, lower.tail =F),5) 
    
    #save p-values
  p.val.hostBG.cds <- 
      tibble( host.gene=cur.sig.host,
              phage.gene=cur.sig.phage,
              host.DEG=m,
              phage.DEG=k,
              shared.DEG=x,
              total=m+n,
              hg.left=hg.left,
              hg.right=hg.right)%>%
      bind_rows(p.val.hostBG.cds,.)
  
    
  }
}


# parameters to plot PDF
p.val.hostBG.cds <- 
  p.val.hostBG.cds %>% 
  mutate(pnl=case_when(phage.gene %in% c("sigF","sigG") ~ "host",
                       TRUE ~ "phage") ) %>% 
    mutate(M = host.DEG,
           N = total-host.DEG,
           K = phage.DEG,
           X = pmin(K,M),
           phage.gene = paste0(pnl,": ", phage.gene),
           host.gene = paste0("overlap.set:", host.gene),
           who = paste(host.gene, phage.gene, sep = "_"))


#adjust p values (joint for left and right)
p.val.hostBG.cds <- p.val.hostBG.cds %>% 
  select(who, hg.left, hg.right) %>% 
  pivot_longer(-who) %>% 
  mutate(adj.p=p.adjust(value,  method = "BH")) %>% 
  pivot_wider(-value, names_from = "name", values_from = "adj.p") %>% 
  rename(adj.hg.left = hg.left, adj.hg.right = hg.right) %>% 
  left_join(p.val.hostBG.cds, .)

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
  filter(value > 1e-8) %>% 
  ggplot(aes(x, value))+
  geom_area(fill = "grey70", color = "grey30")+
  geom_vline(data = p.val.hostBG.cds, aes(xintercept = shared.DEG),
             color = "red", size = 1)+
  geom_text(data = p.val.hostBG.cds, aes(x = shared.DEG, label =  stars.pval(adj.hg.right)), 
            size = 10, y = Inf, vjust = 1, hjust = 1, color = "red")+
  geom_text(data = p.val.hostBG.cds, aes(x = shared.DEG, label =  stars.pval(adj.hg.left)), 
            size = 10, y = Inf, vjust = 1, hjust = 0, color = "blue")+
  facet_wrap(phage.gene ~ host.gene, scales = "free", 
             nrow = 2, dir = 'v', labeller = label_parsed)+
  theme_classic()+
  panel_border(color = "black")+
  labs(caption = paste ("BH adj. P-value:",attr(stars.pval(1),"legend")))+
  theme(plot.caption = element_text(colour = "grey40"))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("shared differentially expressed genes")+
  ylab("hypergeometric PDF")

p
ggsave(here("RNAseq/plots/enrichment_hostBG.png"),plot = p, width = 8, height = 4)

# export analysis test results
write_csv(p.val.hostBG.cds,
          here("RNAseq/data/overlap_enrichment.csv"))
# * sporul.gene enrich ---------------------
#   Are sporulation genes enriched in the sample of deferentially expressed genes?
# the FDR corrected p-values for indicating significantly DExd genes.

# sporulation genes by subtiwiki
d6.spor.cds <- d6.spor.cds %>% 
  mutate(sw.spore  = str_detect(replace_na(category2," "), regex("sporulation", ignore_case = T)))

# record p values
p.val.spore <- tibble()


for(cur.gene in colnames(fc)[-1]){
  # pDR110 has noe DEXed genes. so skippin it to prevent error
  if(cur.gene=="pDR110") next
  #prepare data to make contingency table
  d.urn <- 
    p.val%>%
    # add data on sporulation genes
    right_join(., d6.spor.cds, by=c("id"="locus_tag.d6"))%>%
    # select the data for the current gene
    select(id, gene,pBH=cur.gene, sw.spore) %>% 
    # remove genes for which sporulation status is unknown
    filter(!is.na(sw.spore))
  
  #add data on direction of change
  d.urn <- fc%>%
    # #filter only cds
    right_join(., d6.spor.cds, by=c("id"="locus_tag.d6"))%>%
    # select the data for the current genes
    select(id, fc=cur.gene) %>%
    left_join(d.urn, ., by = "id")
  
  d.urn <- d.urn %>% 
    # change genes that were not assigned a p-value for DE to 1 (no change)
    mutate(pBH=if_else(is.na(pBH),1,pBH))%>%
    # define logical vector of DE based on p-value
    mutate(upregulated= ((pBH<0.05) & (fc > 2)) )%>%
    # change logical vecotrs to meaningful strings
    mutate(upregulated=ifelse(upregulated, "upregulated", "sameORdown"),
           sw.spore=ifelse(sw.spore, "spor.gene", "other"))
  
  #make contingency table 
  urn <- d.urn %>% 
    select(sw.spore, upregulated) %>% 
    table()

  # sporulation genes observed (upregulated)
  x <- urn["spor.gene","upregulated"]
  
  # total sporulation genes
  m<- rowSums(urn)["spor.gene"]
  
  # non-sporulation genes
  n <- rowSums(urn)[["other"]]
  
  #number of upregulated genes
  k <- colSums(urn)[["upregulated"]]
  
  z <- 0:min(k,m)
  
  chi <- capture.output(summary(urn))
  fisher <- signif(fisher.test(urn)$p.value,4)
  hg.left <- signif(phyper(x-1, m, n, k, lower.tail =T),5) 
  hg.right <- signif(phyper(x-1, m, n, k, lower.tail =F),5) 
  
  #save p-values
  p.val.spore <- 
    tibble( gene=cur.gene,
            spor.genes=m,
            DEG=k,
            spor.DEG=x,
            total.gene=n+m,
            hg.left=hg.left,
            hg.right=hg.right,
            fisher=fisher,
            chi=str_extract(chi[4], "value = .*")%>%parse_number())%>%
    bind_rows(p.val.spore,.)
  

}



# parameters to plot PDF
p.val.spore <- 
  p.val.spore %>% 
  mutate(M = spor.genes,
         N = total.gene-spor.genes,
         K = DEG,
         X = pmin(K,M),
         gene =  gene)  
  #adjust p values (joint for left and right)
p.val.spore <- p.val.spore %>% 
  select(gene, hg.left, hg.right) %>% 
    pivot_longer(-gene) %>% 
  mutate(adj.p=p.adjust(value,  method = "BH")) %>% 
    pivot_wider(-value, names_from = "name", values_from = "adj.p") %>% 
  rename(adj.hg.left = hg.left, adj.hg.right = hg.right) %>% 
    left_join(p.val.spore, .)


# calculate PDF per interaction
max.x <- max(p.val.spore$X)
d.hyp <- tibble(x= seq(max.x))

for(i in 1:nrow(p.val.spore)){
  
  d.hyp[, p.val.spore$gene[i]] <- 
    with(p.val.spore,
         dhyper(x = seq(max.x), m = M[i], n = N[i], k = K[i]))
}
# adjustment for panel order
p.val.spore <- 
  p.val.spore %>% 
    mutate(gene = fct_relevel(gene, "SP10", after = 3))
# plot
p <- d.hyp %>% 
  pivot_longer(-x, names_to = "gene") %>% 
  filter(value > 1e-8) %>%
  mutate(gene = fct_relevel(gene, "P[IPTG]-SP10", after = 3)) %>% 
  ggplot(aes(x, value))+
  geom_area(fill = "grey70", color = "grey30")+
  geom_vline(data = p.val.spore, aes(xintercept = spor.DEG),
             color = "red", size = 1)+
  geom_text(data = p.val.spore, aes(x = spor.DEG, label =  stars.pval(adj.hg.right)), 
            size = 10, y = Inf, vjust = 1, hjust = 1, color = "red")+
  geom_text(data = p.val.spore, aes(x = spor.DEG, label =  stars.pval(adj.hg.left)), 
            size = 10, y = Inf, vjust = 1, hjust = 0, color = "blue")+
  facet_wrap(~ gene, scales = "free", nrow = 3, dir = 'h', labeller = label_parsed)+
  theme_classic()+
  panel_border(color = "black")+
  labs(caption = paste ("BH adj. P-value:",attr(stars.pval(1),"legend")))+
  theme(plot.caption = element_text(colour = "grey40"))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("upregulated sporulation genes")+
  ylab("hypergeometric PDF")

p
ggsave(here("RNAseq/plots/enrichment_sporulation.png"),plot = p, width = 4, height = 4)

# export analysis test results
write_csv(p.val.spore,
          here("RNAseq/data/sporulation_gene_enrchment.csv"))

# * SigB enrichment ----------------------------------------------------------




#   Are sigB regulated genes enriched in the sample of deferentially expressed genes?
# the FDR corrected p-values for indicating significantly DExd genes.

# Regulons by subtiwiki
sw.regulons <- 
  read_csv(here("RNAseq/data/annotations/SW_regulations.csv"))

b.regulon <- sw.regulons %>% 
  filter(str_detect(regulon, "SigB")) %>% 
  pull(`locus tag`)

d6.spor.cds <- d6.spor.cds %>% 
  mutate(sigB =  locus_tag.168 %in% b.regulon)

# record p values
p.val.sigB <- tibble()


for(cur.gene in colnames(fc)[-1]){
  # pDR110 has noe DEXed genes. so skippin it to prevent error
  if(cur.gene=="pDR110") next
  #prepare data to make contingency table
  d.urn <- 
    p.val%>%
    # add data on sporulation genes
    right_join(., d6.spor.cds, by=c("id"="locus_tag.d6"))%>%
    # select the data for the current gene
    select(id, gene,pBH=cur.gene, sigB) %>% 
    # remove genes for which sporulation status is unknown
    filter(!is.na(sigB))
  
  #add data on direction of change
  d.urn <- fc%>%
    # #filter only cds
    right_join(., d6.spor.cds, by=c("id"="locus_tag.d6"))%>%
    # select the data for the current genes
    select(id, fc=cur.gene) %>%
    left_join(d.urn, ., by = "id")
  
  d.urn <- d.urn %>% 
    # change genes that were not assigned a p-value for DE to 1 (no change)
    mutate(pBH=if_else(is.na(pBH),1,pBH))%>%
    # define logical vector of DE based on p-value
    mutate(upregulated= ((pBH<0.05) & (fc > 2)) )%>%
    # change logical vecotrs to meaningful strings
    mutate(upregulated=ifelse(upregulated, "upregulated", "sameORdown"),
           sigB=ifelse(sigB, "sigB.reg", "other"))
  
  #make contingency table 
  urn <- d.urn %>% 
    select(sigB, upregulated) %>% 
    table()
  
  # sporulation genes observed (upregulated)
  x <- urn["sigB.reg","upregulated"]
  
  # total sporulation genes
  m<- rowSums(urn)["sigB.reg"]
  
  # non-sporulation genes
  n <- rowSums(urn)[["other"]]
  
  #number of upregulated genes
  k <- colSums(urn)[["upregulated"]]
  
  z <- 0:min(k,m)
  
  chi <- capture.output(summary(urn))
  fisher <- signif(fisher.test(urn)$p.value,4)
  hg.left <- signif(phyper(x-1, m, n, k, lower.tail =T),5) 
  hg.right <- signif(phyper(x-1, m, n, k, lower.tail =F),5) 
  
  #save p-values
  p.val.sigB <- 
    tibble( gene=cur.gene,
            sigB.reg=m,
            DEG=k,
            sigB.DEG=x,
            total.gene=n+m,
            hg.left=hg.left,
            hg.right=hg.right,
            fisher=fisher,
            chi=str_extract(chi[4], "value = .*")%>%parse_number())%>%
    bind_rows(p.val.sigB,.)
  
  
}



# parameters to plot PDF
p.val.sigB <- 
  p.val.sigB %>% 
  mutate(M = sigB.reg,
         N = total.gene-sigB.reg,
         K = DEG,
         X = pmin(K,M),
         gene = paste0("P[IPTG]-", gene))  
#adjust p values (joint for left and right)
p.val.sigB <- p.val.sigB %>% 
  select(gene, hg.left, hg.right) %>% 
  pivot_longer(-gene) %>% 
  mutate(adj.p=p.adjust(value,  method = "BH")) %>% 
  pivot_wider(-value, names_from = "name", values_from = "adj.p") %>% 
  rename(adj.hg.left = hg.left, adj.hg.right = hg.right) %>% 
  left_join(p.val.sigB, .)


# calculate PDF per interaction
max.x <- max(p.val.sigB$X)
d.hyp <- tibble(x= seq(max.x))

for(i in 1:nrow(p.val.sigB)){
  
  d.hyp[, p.val.sigB$gene[i]] <- 
    with(p.val.sigB,
         dhyper(x = seq(max.x), m = M[i], n = N[i], k = K[i]))
}
# adjustment for panel order
p.val.sigB <- 
  p.val.sigB %>% 
  mutate(gene = fct_relevel(gene, "P[IPTG]-SP10", after = 3))
# plot
p <- d.hyp %>% 
  pivot_longer(-x, names_to = "gene") %>% 
  filter(value > 1e-8) %>%
  mutate(gene = fct_relevel(gene, "P[IPTG]-SP10", after = 3)) %>% 
  ggplot(aes(x, value))+
  geom_area(fill = "grey70", color = "grey30")+
  geom_vline(data = p.val.sigB, aes(xintercept = sigB.DEG),
             color = "red", size = 1)+
  geom_text(data = p.val.sigB, aes(x = sigB.DEG, label =  stars.pval(adj.hg.right)), 
            size = 10, y = Inf, vjust = 1, hjust = 1, color = "red")+
  geom_text(data = p.val.sigB, aes(x = sigB.DEG, label =  stars.pval(adj.hg.left)), 
            size = 10, y = Inf, vjust = 1, hjust = 0, color = "blue")+
  facet_wrap(~ gene, scales = "free", nrow = 3, dir = 'h', labeller = label_parsed)+
  theme_classic()+
  panel_border(color = "black")+
  labs(caption = paste ("BH adj. P-value:",attr(stars.pval(1),"legend")))+
  theme(plot.caption = element_text(colour = "grey40"))+
  scale_y_continuous(expand = c(0, 0))+
  xlab("upregulated sigB genes")+
  ylab("hypergeometric PDF")

p
ggsave(here("RNAseq/plots/enrichment_sigB.png"),plot = p, width = 4, height = 4)

# END ---------------------