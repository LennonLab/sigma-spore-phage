library(here)
library(tidyverse)
library(cowplot)

load(here("RNAseq/plots/volcano-plot.Rdata"))
pA <- p

load(here("RNAseq/plots/correlations-plot.Rdata"))
pB <- p
bottom_row <- plot_grid(NULL,pB,NULL,
                        rel_widths = c(0.1,1,.1), nrow = 1)

p <- plot_grid(pA, bottom_row,NULL, labels = letters[1:2], ncol = 1,
          rel_heights = c(1, 0.66, 0.05))

ggsave(here("RNAseq/plots/main.png"),plot = p, width = 6, height = 6)
