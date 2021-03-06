Virulence index
================

``` r
library(here)
```

    ## here() starts at C:/Users/danschw/GitHub/sigma-spore-phage

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.5     v purrr   0.3.4
    ## v tibble  3.1.6     v dplyr   1.0.8
    ## v tidyr   1.2.0     v stringr 1.4.0
    ## v readr   2.1.2     v forcats 0.5.1

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(lubridate)
```

    ## 
    ## Attaching package: 'lubridate'

    ## The following objects are masked from 'package:base':
    ## 
    ##     date, intersect, setdiff, union

``` r
library(broom)
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'

    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
library(scales)
```

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

``` r
library(growthcurver)
library(caTools) #needed for trapezoids
source(here("virulence/code/virulence_functions.R"))
```

``` r
init.host <- 1e8 #CFU/ml
host.vol <- 0.1 #ml
init.titer <- 1e9 #PFU/ml
phage.vol <- 0.1 #ml
```

# The virulence index

Based on paper *Storms, Zachary J., et al. “The Virulence Index: A
Metric for Quantitative Analysis of Phage Virulence.” PHAGE 2019*

I want to commpare SP10 WT to the
![\\Delta g120](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20g120 "\Delta g120")
on both the WT host
(![\\Delta 6](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%206 "\Delta 6"))
anf its
![\\Delta 6\\Delta sigF](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%206%5CDelta%20sigF "\Delta 6\Delta sigF")
derivative.

I ran the experiments twice. In the first run I made a mistake in the
columns to which the phages were added. This resulted in having complete
data from that experiment only for SP10
![\\Delta g120](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5CDelta%20g120 "\Delta g120").
Below I analyze the data from both experiments.

## load and organize data

colonies of first experiment labeled \#1 and \#2. From second experiment
labeled \#3 and \#4

``` r
# read file experiment 1
input <- 
  read_synergy_txt(here("virulence/data/20190823_SP10_VIndx.txt")) %>% 
  # In these files there are reads every 2 min. 
  # I will average on window of 8 reads
  smooth_synergy_txt(raw_input = . , window_size = 8)

#transform to long format
d <- gather(input, key="well", value = "OD600", colnames(input)[-c(1:2)])

#add metadata
meta <- read_csv(here("virulence/data/20190823_meta.csv"))
```

    ## Rows: 96 Columns: 7
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr (4): well, row, host, phage
    ## dbl (3): col, colony, dilution
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
d <-merge(d, meta, by = "well")

# read file experiment 1
input <- 
  read_synergy_txt(here("virulence/data/20190826_SP10_VIndx.txt")) %>% 
  # In these files there are reads every 2 min. 
  # I will average on window of 8 reads
  smooth_synergy_txt(raw_input = . , window_size = 8)

#add metadata
meta <- read_csv(here("virulence/data/20190826_meta.csv"))
```

    ## Rows: 96 Columns: 7
    ## -- Column specification --------------------------------------------------------
    ## Delimiter: ","
    ## chr (4): well, row, host, phage
    ## dbl (3): col, colony, dilution
    ## 
    ## i Use `spec()` to retrieve the full column specification for this data.
    ## i Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#transform to long format
d <- 
  gather(input, key="well", value = "OD600", colnames(input)[-c(1:2)]) %>% 
  merge(., meta, by = "well") %>% 
  bind_rows(d, .)
```

## Examine and clean data

``` r
d%>%
  # filter(host != "blank")%>%
  ggplot(aes(x=Time, y=OD600))+
  geom_line(aes(group= interaction(well,colony),
                color=as.character(dilution)), size=1)+
  facet_grid(host~phage)+
  theme_cowplot()+
  panel_border()+
  scale_color_viridis_d()+
  theme(legend.position = "bottom")
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Overall, there is a good dilution dependence in the timing of bacterial
growth reduction (i.e. lysis). There appear to be a contaminated blank
control, and perhaps some wells that did not receive phage in the
highest phage dilutions. This is expected.

``` r
d%>%
  filter(host == "blank")%>%
  ggplot(aes(x=Time, y=OD600))+
  geom_line( size=1)+
  facet_wrap(well~colony, strip.position = "r" )+
  theme_classic()+
  panel_border()+
  scale_color_viridis_d()+
  theme(legend.position = "bottom")
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

single blank contamination. Not bad. Moving on.

## Calculating the viruence index

Folowing Storms et. al.  
*1. establish limit of integration* *“It is important to stress how the
establishment of the limit of integration plays a significant role in
the assessment of virulence. This limit should be set as the onset of
stationary phase in the phage-free control. This provides a consistent
reference for integration that can be easily identified for any
phage–host system and restricts measurements to the period of cell
growth—a necessary condition for productive infection for many phages.47
Moreover, it ensures that the range of the virulence measurements is
well distributed, as discussed hereunder. In general, we recommend
establishing the limit of integration as the time at which the slope of
OD630 over time reaches ≤0.03 h.”*

### Determining the integration limit

Doing this can be done in 2 ways:  
\#\#\#\# A.using ‘growthcurver’ to find no-phage carying capacity (k)

``` r
k.noPHI <- 
  d %>%
  filter(host != "blank")%>%
  filter(phage=="noPHI") %>% 
  mutate(well_col = paste(well, colony, sep="_")) %>% 
  select(well_col, Time, OD600) %>% 
  pivot_wider(names_from = well_col, values_from = OD600) %>% 
  #focus on time before OD decline
  filter(Time < 7) %>%
  # derive growth parameters with "growthcurver"   
  SummarizeGrowthByPlate(., bg_correct = "none")%>%
  select(well_col = sample, k)

# add data on host colony
k.noPHI <-
  d %>%
  mutate(well_col = paste(well, colony, sep="_")) %>% 
  select(well_col, host, colony) %>% 
  distinct() %>% 
  left_join(k.noPHI, ., by = "well_col")

# for each of the non-infected wells we now find the time it reached carrying capacity

# add column to store time
k.noPHI$Time.k <- NA

for(i in seq(k.noPHI$well_col)){
  w <- k.noPHI$well_col[i]
  k <- k.noPHI$k[i]
  tmp <- 
    d%>%
    mutate(well_col = paste(well, colony, sep="_")) %>% 
    filter(well_col==w) %>% 
    arrange(Time)
  
  k.noPHI$Time.k[i] <-
    tmp$Time[which(tmp$OD600>k)[1]]
  
}
#remove loop vars
rm(tmp, i,w,k)

# use the median time found
int.limit <- median(k.noPHI$Time.k, na.rm = TRUE)

k.noPHI %>% 
  ggplot(aes(host, Time.k)) +
  geom_boxplot(fill="grey80")+
  geom_jitter(width = 0.1, shape=21, size=2, fill="white")+
  geom_hline(yintercept = int.limit)+
  labs(caption = paste("median time to K = ", int.limit, "h"))+
  theme_classic()+
  panel_border(color = "black")
```

    ## Warning: Removed 3 rows containing non-finite values (stat_boxplot).

    ## Warning: Removed 3 rows containing missing values (geom_point).

![](virulence_SP10_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

#### B. by visual inspection of plots

First I will look at the change in OD over time. Storms et al. recommend
*establishing the limit of integration as the time at which the slope of
OD630 over time reaches ≤0.03 h*

``` r
time.diff <- diff(d$Time[1:2])

d %>% 
  filter(phage=="noPHI")%>%
  mutate(well_col = paste(well, colony, sep="_")) %>% 
  split(.$well_col)%>% 
  map_df("OD600")%>%
  map_df(diff)%>%
  map_df(function(x) x/time.diff)-> tmp

matplot(x= unique(d$Time)[-60],tmp, type="l", ylab = "OD/Time diff")
abline(h=0.03, col="red", lwd=3)
abline(v=c(1:floor(max(d$Time))), col="grey")
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

This method show that the integrating limit should be at \~ 5 hours.
Taking both together I will set integration limit at 6 hours. How does
that look?

``` r
int.limit <- 6

d%>%
  filter(phage=="noPHI")%>%
  filter(colony %in% c(1:2))%>%
  filter(Time<int.limit)%>%
  ggplot(aes(x=Time, y=OD600))+
  geom_line(data=filter(d,phage=="noPHI") %>% 
              filter(colony %in% c(1:2)), color="grey", size=1)+
  geom_line(aes(color=host), size=1)+
  facet_wrap(~well)+
  theme_classic()
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
d%>%
  filter(phage=="noPHI")%>%
  filter(colony %in% c(3:4))%>%
  filter(Time<int.limit)%>%
  ggplot(aes(x=Time, y=OD600))+
  geom_line(data=filter(d,phage=="noPHI") %>% 
              filter(colony %in% c(3:4)), color="grey", size=1)+
  geom_line(aes(color=host), size=1)+
  facet_wrap(~well)+
  theme_classic()
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-7-2.png)<!-- -->

Looks OK to me.

## Choosing the dilution range

Need to make sure to use only dilutions in which all wells got phage.

``` r
max.dilut <- 1/(init.titer*phage.vol)
#-# verify visualy
d%>%
  mutate(well_col = paste(well, colony, sep="_")) %>% 
  filter(host != "blank")%>%
  filter(phage!="noPHI")%>%
  ggplot(aes(x=Time, y=OD600))+
  geom_line(aes(color=well_col), size=1)+
  facet_grid(dilution~interaction(host, phage))+
  theme_classic()+
  panel_border(color = "black")+
  theme(legend.position = "none")
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

I think 1e-5 is the limit, maybe 1e-6. Closer look:

``` r
d%>%
  filter(host != "blank")%>%
  filter(phage!="noPHI")%>%
  filter(dilution==1e-7 |dilution==1e-6 |dilution==1e-5)%>%
  ggplot(aes(x=Time, y=OD600))+
  geom_line(aes(color=well), size=1)+
  facet_grid(dilution~interaction(host, phage))+
  theme_cowplot()+
    theme(legend.position = "none")
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

``` r
max.dilut <-1e-6
```

I think 1e-6 is pretty safe to say all dilutions recieved phage.

## Integrate the area under curve

Storms et al.: *"areas underneath the optical density versus time curves
were calculated using the trapezoid rule for each well, from the time of
infection to the time corresponding to the onset of stationary phase in
the phage-free control*

``` r
meta <- 
  d %>% 
  select(-OD600, -Time, -temp) %>% 
  distinct() %>% 
  mutate(well_col = paste(well, colony, sep="_"))

meta$auc <- NA
for(i in seq(meta$well_col)){
  if(meta$host[i]=="blank") next

  tmp <- 
    d%>%
    mutate(well_col = paste(well, colony, sep="_")) %>% 
    filter(well_col==meta$well_col[i])%>%
    filter(Time<=int.limit)
  meta$auc[i] <- trapz(tmp$Time, tmp$OD600)
    
}
rm(tmp)

# summarize no phage control areas
sum.noPHI <- 
  meta%>%
    filter(host != "blank")%>%
    filter(phage=="noPHI")%>%
    group_by(host, colony)%>%
    summarise( A0=mean(auc), sd=sd(auc))
```

    ## `summarise()` has grouped output by 'host'. You can override using the
    ## `.groups` argument.

``` r
# sum.noPHI %>%
  # ggplot(aes(colony, A0))+
  # geom_pointrange(aes(ymin=A0-sd, ymax=A0+sd))+
  # facet_wrap(~host)

vindex <- merge(meta, sum.noPHI)
vindex$Vi <- 1-(vindex$auc/vindex$A0)
vindex$moi <- (init.titer * phage.vol * vindex$dilution) /
                              (init.host * host.vol)
vindex$log.moi <- log10(vindex$moi)

vindex%>%
  filter(phage != "noPHI") %>% 
  filter(dilution>=max.dilut)%>%
  group_by(host,colony, phage, log.moi)%>%
  summarise( Virulence=mean(Vi), mn=min(Vi), mx=max(Vi), n=n())%>%
  mutate(colony = as.character(colony)) %>% 
    ggplot(aes(log.moi, Virulence, color=colony))+
      geom_line(aes(group=colony), size=1)+
      geom_pointrange(aes(ymin=mn, ymax=mx ), shape=21, size=.5, fill = "white")+
      facet_wrap(phage~host, nrow=2,labeller = "label_both")+
      theme_classic()+
      panel_border(color="black")+
      scale_colour_viridis_d()
```

    ## `summarise()` has grouped output by 'host', 'colony', 'phage'. You can override
    ## using the `.groups` argument.

![](virulence_SP10_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

## summarize to single value index

Storms et al.: *The virulence index is defined as the area under the
virulence curve (AP) divided by the theoretical maximum area under the
virulence curve (Amax)* **Vp = Ap / Amax**

``` r
sum.phi <- vindex%>%
  filter(phage!="noPHI")%>%
  group_by(host,colony, phage)%>%
  summarise(  n=n())
```

    ## `summarise()` has grouped output by 'host', 'colony'. You can override using
    ## the `.groups` argument.

``` r
sum.phi$Ap <- NA
sum.phi$Amax <- NA
for(i in seq(nrow(sum.phi))){
 
  tmp <- 
    vindex%>%
    filter(phage!="noPHI")%>%
    filter(host==sum.phi$host[i])%>%
    filter(colony==sum.phi$colony[i])%>%
    filter(phage==sum.phi$phage[i])%>%
    filter(dilution>=max.dilut)%>%
    arrange(log.moi)
  sum.phi$Ap[i] <- trapz(tmp$log.moi, tmp$Vi)
  sum.phi$Amax[i] <-  trapz(tmp$log.moi,rep(1, nrow(tmp)))
  
}

rm(tmp)
sum.phi$Vp <- sum.phi$Ap/sum.phi$Amax

sum.phi%>%
  # filter(n==8)%>%
  ggplot(aes(x=interaction(phage,host), y=Vp))+
  geom_jitter(aes(shape=as.factor(colony)),
              height=0, width=0.1,size=3, stroke=2, fill="white")+
  scale_shape_manual(values=c(15,16,22,21))+
  ylim(0,1)+
  theme_cowplot()
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-11-1.png)<!-- --> \#\#
Stats

``` r
summary(aov(Vp~phage+host+colony, sum.phi  ))
```

    ##             Df   Sum Sq   Mean Sq F value Pr(>F)
    ## phage        1 0.001733 0.0017326   1.563  0.235
    ## host         1 0.002637 0.0026373   2.379  0.149
    ## colony       1 0.000377 0.0003766   0.340  0.571
    ## Residuals   12 0.013301 0.0011084

*No Difference *

``` r
t.test(Vp~phage, sum.phi)
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  Vp by phage
    ## t = 1.2194, df = 12.886, p-value = 0.2446
    ## alternative hypothesis: true difference in means between group del120 and group wt is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.01609505  0.05772004
    ## sample estimates:
    ## mean in group del120     mean in group wt 
    ##            0.2619106            0.2410981

*No Difference *

# plot

``` r
library(ggsignif)
p.final <- sum.phi%>%
  group_by(phage)%>%
  summarise(n=n(), sd=sd(Vp), Vp=mean(Vp), se=sd/sqrt(n))%>%
  # filter(host=="wt")%>%
  ggplot(aes(x=phage, y=Vp))+
  geom_crossbar(aes(ymin=Vp-se,ymax=Vp+se), width=0.3)+
  geom_jitter(data=sum.phi,aes(x=phage, y=Vp), shape = 21,
              height=0, width=0.15,size=3, fill=alpha("grey", 0.5))+
    geom_signif(comparisons = list(c("del120", "wt")), 
                y_position = 0.4, tip_length = 0.3,
                annotations ="NS")+
  scale_shape_manual(values = c(21,23))+
  ylim(0,0.5)+
  ylab("phage virulence (Vp)")+
  theme_classic()+
  panel_border(color = "black")
  # ggsave(here("virulence/plots/virulence_SP10.png"),
  #        width = 3, height = 3)
p.final
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# plot of MOI dependent lysis


d.plot <- d %>%
  filter(host != "blank")%>%
  filter(dilution>=max.dilut)%>%
  mutate(dilution = if_else(phage == "noPHI", 0, dilution)) %>% 
  mutate(moi = (init.titer * phage.vol * dilution) /
                              (init.host * host.vol)) %>% 
  group_by(Time, phage, moi) %>% 
  summarise(n = n(), m = mean(OD600), v = sd(OD600)/sqrt(n), .groups = "drop") %>% 
   mutate(moi = as_factor(moi) %>% fct_inseq()) 

# duplicate no phage data to each of the panels
d.plot <- bind_rows(
  filter(d.plot, phage=="noPHI") %>% mutate(phage = "wt"),
  filter(d.plot, phage=="noPHI") %>% mutate(phage = "del120"),
  filter(d.plot, phage!="noPHI")
)

p.lysis <- d.plot %>% 
  filter(Time <= 10) %>% 
  ggplot(aes(x=Time, y=m, color=moi))+
  geom_vline(xintercept = int.limit, color = "grey")+
  geom_linerange(aes(ymin = m-v, ymax = m+v), alpha = 0.5, show.legend = F)+
  geom_line(size=1)+
  #highlight no phage
   geom_line(data = filter(d.plot, moi=="0" & Time <= 10), size=1.5)+
  facet_wrap(~ phage)+
  theme_classic()+
  panel_border(color = "black")+
  scale_color_viridis_d()+
  guides(color = guide_legend("MOI", ncol =2))+
  ylab("Bacterial Density (OD600)")+
  xlab("Time Post Infection (h)")+
  scale_x_continuous(breaks = seq(0,10,2))+
  theme(panel.spacing = unit(2, "lines"),
        plot.margin = unit(c(0.5,1,0.5,1), "cm"))

p.lysis 
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->
Showing virulence index calculation

``` r
# dunction to adjust facet labels
f_labeller <- function(string){return(paste0("MOI = ", string))}

# duplicate no phage data to each of the panels
mois <- levels(d.plot$moi)
d.plot2 <- 
  filter(d.plot, moi != "0") %>% 
  mutate(moi_panel = moi)  
  
for (i in mois){
  
  if(i == "0") next
  
  d.plot2 <- 
    bind_rows(
      d.plot2,
      filter(d.plot, moi == "0") %>% mutate( moi_panel = i),
    )
}

#panel order
d.plot2 <- d.plot2%>% 
  mutate(moi_panel = as_factor(moi_panel) %>% fct_inseq()) %>% 
  # only plotting part of the curved used for index
  filter(Time <= int.limit)


#separate plotting data sets
d.plot2.wtInf <- d.plot2 %>% 
  filter(moi != "0") %>% 
  filter(phage == "wt") 

d.plot2.wtCtrl <- d.plot2 %>% 
  filter(moi == "0") %>% 
  filter(phage == "wt") 

d.plot2.mutInf <- d.plot2 %>% 
  filter(moi != "0") %>% 
  filter(phage == "del120") 

d.plot2.mutCtrl <- d.plot2 %>% 
  filter(moi == "0") %>% 
  filter(phage == "del120") 

p.wt <- d.plot2.wtCtrl %>% 
  ggplot(aes(x=Time, y=m))+
  geom_area(fill = "grey70")+
  geom_area(data = d.plot2.wtInf, fill = "white" )+
  geom_line(linetype=2)+
  geom_line(data = d.plot2.wtInf)+
  facet_wrap(~ moi_panel, labeller = labeller(moi_panel = f_labeller))+
  theme_classic()+
  panel_border(color = "black")+
   theme(panel.spacing = unit(2, "lines"))+
  scale_color_viridis_d()+
  guides(color = guide_legend("multiplicity\nof\ninfection"))+
  ylab("Bacterial Density (OD600)")+
  xlab("Time Post Infection (h)")+
  scale_x_continuous(breaks = seq(0,10,2))+
  ggtitle("SP10 WT")

p.mut <- d.plot2.mutCtrl %>% 
  ggplot(aes(x=Time, y=m))+
  geom_area(fill = "grey70")+
  geom_area(data = d.plot2.mutInf, fill = "white" )+
  geom_line(linetype=2)+
  geom_line(data = d.plot2.mutInf)+
  facet_wrap(~ moi_panel, labeller = labeller(moi_panel = f_labeller))+
  theme_classic()+
  panel_border(color = "black")+
   theme(panel.spacing = unit(2, "lines"))+
  scale_color_viridis_d()+
  guides(color = guide_legend("multiplicity\nof\ninfection"))+
  ylab("Bacterial Density (OD600)")+
  xlab("Time Post Infection (h)")+
  scale_x_continuous(breaks = seq(0,10,2))+
  ggtitle("SP10 del120")

p.auc <- plot_grid(p.mut, NULL,p.wt, rel_widths = c(1,0.1,1), nrow = 1)
p.auc
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
p.vi <- vindex%>%
  filter(phage != "noPHI") %>% 
  filter(dilution>=max.dilut)%>%
  group_by(phage, moi)%>%
  summarise( Virulence=mean(Vi), mn=min(Vi), mx=max(Vi), n=n())%>%
    ggplot(aes(moi, Virulence))+
  geom_area(fill = "grey70")+
      geom_line(size=1)+
      geom_pointrange(aes(ymin=mn, ymax=mx ), shape=21, size=.5, fill = "white")+
      facet_wrap(~phage, nrow=1,labeller = "label_both")+
      theme_classic()+
      panel_border(color="black")+
      scale_colour_viridis_d()+
  scale_x_log10()+
  ylab("local virulence\n(area between curves)")+
  ylim(NA,1)
```

    ## `summarise()` has grouped output by 'phage'. You can override using the
    ## `.groups` argument.

``` r
p.vi
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

\#\#\#combine plots

``` r
p.bottom <- plot_grid(p.vi, NULL, p.final, rel_widths = c(1,0.05,1), nrow = 1,
                      labels = c("", "d"))
p.all <- plot_grid(p.lysis,NULL,p.auc,NULL,p.bottom, ncol = 1,
          rel_heights = c(1,0.05,2,0.05,1), labels = c("a","","b","","c"))

#add white background
p.all <-ggdraw(p.all) + 
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here("virulence/plots", "virulence_steps.png"), p.all,
      width = 8, height = 8)

    # #export to pptx using officer and rvg
    # library (officer)
    # library(rvg)
    # 
    # read_pptx() %>%
    #   add_slide(layout = "Blank", master = "Office Theme" ) %>%
    #   ph_with(dml(ggobj = p.all),
    #           location = ph_location(type = "body",
    #                                  left = 0, top = 0,
    #                                  width = 8, height = 8)) %>%
    #   print(target = here("virulence/plots", "virulence_steps.pptx"))
p.all
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-18-1.png)<!-- --> \#
main figure

``` r
p.bottom <- plot_grid(NULL, p.final,NULL, rel_widths = c(0.2,1,1.5), nrow = 1)
p.all <-  plot_grid(p.lysis,p.bottom, ncol = 1, labels = c("a", "b"))

#add white background
p.all <-ggdraw(p.all) + 
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave(here("virulence/plots", "virulence_main2.png"), 
      p.all,
      width = 6, height = 4)

    #export to pptx using officer and rvg
    library (officer)
    library(rvg)

    read_pptx() %>%
      add_slide(layout = "Blank", master = "Office Theme" ) %>%
      ph_with(dml(ggobj = p.all),
              location = ph_location(type = "body",
                                     left = 0, top = 0,
                                     width = 6, height = 4)) %>%
      print(target = here("virulence/plots", "virulence_main.pptx"))
p.all
```

![](virulence_SP10_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->
