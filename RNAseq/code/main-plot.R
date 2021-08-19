library(here)
library(tidyverse)
library(cowplot)

load(here("RNAseq/plots/volcano-plot.Rdata"))
pA <- p

load(here("RNAseq/plots/correlations-plot.Rdata"))
pB <- p
# bottom_row <- plot_grid(NULL,pB,NULL,
#                         rel_widths = c(0.1,1,.1), nrow = 1)

p <- plot_grid(pA, NULL, pB, labels = c("a","","b"), nrow = 1,
          rel_widths = c(1, .02 ,0.5))

ggsave(here("RNAseq/plots/main.png"),plot = p, width = 8, height = 6)


