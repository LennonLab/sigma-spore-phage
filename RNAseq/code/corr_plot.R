library(here)
library(tidyverse)
library(cowplot)
library(gtools)
library(scales)


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
sig.phage.v <- c("ELDg168","ELDg169","Goe3","SP10")#,"sigF","sigG")

max.fc.log2 <- max(fc[-1], na.rm = T) %>% log2() %>% round()
min.fc.log2 <- min(fc[-1], na.rm = T) %>% log2() %>% round()
# brx <- c(-15, -10,-5,0,5,10,15)
brx <- c( -10,0,10)

# rsave plots
l.plots <- list()

for (cur.sig.host in sig.host.v){ 
  for (cur.sig.phage in sig.phage.v){
    
    if(cur.sig.host == cur.sig.phage) next
    # get loci DExed in both strains
    d.sig <- p.val%>%
      #filter only cds
      filter(str_detect(id, "A8O17"))%>%
      # select the data for the current genes
      select(id,sig.host=all_of(cur.sig.host),sig.phage=all_of(cur.sig.phage))%>%
      # change genes that were not assigned a p-value for DE to 1 (no change)
      mutate(sig.host=if_else(is.na(sig.host),1,sig.host))%>%
      mutate(sig.phage=if_else(is.na(sig.phage),1,sig.phage))%>%
      # define logical vector of DE based on p-value
      mutate(sig.host = sig.host <0.05)%>%
      mutate(sig.phage = sig.phage <0.05) %>% 
      mutate(signifcance  =  case_when(sig.host & sig.phage ~ "both",
                                    sig.host ~ "host only",
                                    sig.phage ~ "phage only",
                                    TRUE ~ "neither")) %>% 
      mutate(signifcance = fct_relevel(signifcance, "host only","phage only","both", "neither"))
      
      d.plot <- 
      # based on list of loci get fold-change data for this pair
        fc %>% 
        select(id,fc.host=all_of(cur.sig.host), fc.phage=all_of(cur.sig.phage)) %>% 
        filter(!is.na(fc.host)) %>% 
        filter(!is.na(fc.phage)) %>% 
        left_join(d.sig, ., by="id")
      
      
      # plot
      plot.name <- paste(cur.sig.phage,cur.sig.host, sep=" vs. ")
      l.plots[[plot.name]] <- 
      d.plot %>% 
        mutate(pnl = plot.name) %>% 
        ggplot(aes(x=log2(fc.host), y=log2(fc.phage)))+
        geom_hline(yintercept = 0, color="grey20")+
        geom_vline(xintercept = 0, color="grey20")+
        geom_abline(slope = 1, intercept = 0, linetype=2, color="lightsteelblue")+
        geom_point(aes(color = signifcance),alpha = 0.5)+
        geom_smooth(method = 'lm', formula = 'y ~ x', color = "black")+
        scale_y_continuous(breaks = brx,
                           limits = c(min.fc.log2,max.fc.log2))+
        scale_x_continuous(breaks = brx,
                           limits = c(min.fc.log2,max.fc.log2))+
        xlab(bquote(.(cur.sig.host)~(log[2]~FC))) +
        ylab(bquote(.(cur.sig.phage)~(log[2]~FC))) +
        facet_wrap(~pnl)+
        theme_classic(base_size = 12)+
        panel_border(color = "black")+
        scale_color_viridis_d()+
        theme(legend.position="right")
  }
}


# main figure -------------------------------------------------------------


# From: https://wilkelab.org/cowplot/articles/shared_legends.html
# arrange the two plots in a single row
to.plot <- names(l.plots)[str_detect(names(l.plots), "169")]
prow <- plot_grid(
  l.plots[[to.plot[1]]] + theme(legend.position="none"),
  l.plots[[to.plot[2]]] + theme(legend.position="none"),
  align = 'vh',
  # labels = c("A", "B"),
  hjust = -1,
  nrow = 2
)

# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  l.plots[[1]] + theme(legend.box.margin = margin(12, 0, 0, 0),
                       legend.position = "bottom",
                       legend.text = element_text(size=14),
                       legend.title = element_text(size = 14))+
    labs(colour="differential\nexpression\nsignificance") +
    guides(color = guide_legend(nrow = 2,
                                override.aes = list(size = 4, alpha = 1)))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
p <- plot_grid(prow, legend, rel_heights =  c(2, .5), ncol = 1)

ggsave(here("RNAseq/plots/correlations.png"),plot = p, width = 3, height = 6)

save(p, file = here("RNAseq/plots/correlations-plot.Rdata"))

# supl. figure -------------------------------------------------------------


# From: https://wilkelab.org/cowplot/articles/shared_legends.html
# arrange the two plots in a single row
to.plot <- names(l.plots)[!str_detect(names(l.plots), "169")]
prow <- plot_grid(
  l.plots[[to.plot[1]]] + theme(legend.position="none"),
  l.plots[[to.plot[2]]] + theme(legend.position="none"),
  l.plots[[to.plot[3]]] + theme(legend.position="none"),
  l.plots[[to.plot[4]]] + theme(legend.position="none"),
  l.plots[[to.plot[5]]] + theme(legend.position="none"),
  l.plots[[to.plot[6]]] + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 3, byrow = F
)



# extract the legend from one of the plots
# use legnd from above
  # legend <- get_legend(
  #   # create some space to the left of the legend
  #   l.plots[[1]] + theme(legend.box.margin = margin( 12, 0, 0))
  # )

# add the legend to the row we made earlier.
p <- plot_grid(prow, legend, rel_heights = c(4, .4), ncol = 1)


ggsave(here("RNAseq/plots/correlations-supl.png"),plot = p, width = 6, height = 8)

