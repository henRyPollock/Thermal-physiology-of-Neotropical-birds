### exploring thermal phys data from Henry Pollock
## 03.02.21 WJC
## 03.14.22 WJC - UPDATED for new dataset from cameron and final version of figures
# Update August 2022 - WJC - 
#     Collapse stratum to 3 categories matching environmental data
#     Recalculate warm_tol & therm_sm
# Update October 2022 - WJC
#     Remove RMR and BMR from figures

# load libraries
library(ggstatsplot)
library(ggplot2)
library(tidyr)
library(reshape2)
library(gridExtra)

dat<-read.csv("/Users/justincooper/OneDrive - George Mason University/Thermo-physiology/thermal phys data_4.20.21_WJC.csv")
str(dat)

# collapse stratum to understory, midstory, and canopy
unique(dat$stratum)
dat[dat$stratum=="generalist","stratum"]<-"midstory"
dat[
  dat$stratum=="terrestrial"|
    dat$stratum=="understory"|
    dat$stratum=="near-ground",
  "stratum"
]<-"understory"

# Recalculate warm_tol for stratum max temp
dat[dat$stratum=="","warm_tol"]<-dat[dat$stratum=="","utl"]-30
dat[dat$stratum=="understory","warm_tol"]<-dat[dat$stratum=="understory","utl"]-28.2
dat[dat$stratum=="midstory","warm_tol"]<-dat[dat$stratum=="midstory","utl"]-29.3
dat[dat$stratum=="canopy","warm_tol"]<-dat[dat$stratum=="canopy","utl"]-29.5

# Recalculate therm_sm for stratum max temp
dat[dat$stratum=="","therm_sm"]<-dat[dat$stratum=="","uct"]-30
dat[dat$stratum=="understory","therm_sm"]<-dat[dat$stratum=="understory","uct"]-28.2
dat[dat$stratum=="midstory","therm_sm"]<-dat[dat$stratum=="midstory","uct"]-29.3
dat[dat$stratum=="canopy","therm_sm"]<-dat[dat$stratum=="canopy","uct"]-29.5


# group edge into second growth
dat[dat$habitat_final=="edge","habitat_final"]<-"second-growth"
# re-organize factor levels
dat$habitat_final<-factor(dat$habitat_final,levels = c("open","second-growth", "forest_edge","forest_interior"),
                          labels=c("Open","Second-Growth", "Forest Edge","Forest Interior"))
dat$stratum<-factor(dat$stratum,levels=c("understory","midstory","canopy"),
                    labels = c("Understory","Midstory","Canopy"))


# setup for ggstatsplot with habitat type
temp<-melt(dat,id.vars = "habitat_final",measure.vars = colnames(dat[,c(12:17,19:21)]))
set.seed(123)

# remove body mass, BMR, & RMR
# p.1<-ggstatsplot::ggbetweenstats(
#   data = dplyr::filter(
#     temp,
#     variable %in% c("mb")),
#   x = habitat_final,
#   y = value,
#   ggsignif.args = list(textsize = 2, tip_length = 0.01),
#   p.adjust.method = "bonferroni", 
#   type = "non-parametric",
#   palette = "default_jama",
#   package = "ggsci",
#   ylab=c("Body Mass (g)"),
#   xlab=" ",
#   results.subtitle = FALSE,
#   title = "A)",
#   centrality.plotting = FALSE,
#   pairwise.comparisons = FALSE,
#   ggplot.component = list(theme(axis.title = element_text(size = 12),
#                                 axis.text = element_text(size=10),
#                                 plot.title = element_text(size = 16)))
# )
# p.1
# p.2<-ggstatsplot::ggbetweenstats(
#   data = dplyr::filter(
#     temp,
#     variable %in% c("bmr")),
#   x = habitat_final,
#   y = value,
#   ggsignif.args = list(textsize = 2, tip_length = 0.01),
#   p.adjust.method = "bonferroni", 
#   type = "non-parametric",
#   palette = "default_jama",
#   package = "ggsci",
#   ylab=c("BMR"),
#   xlab=" ",
#   results.subtitle = FALSE,
#   title = "A)",
#   centrality.plotting = FALSE,
#   pairwise.comparisons = FALSE,
#   ggplot.component = list(theme(axis.title = element_text(size = 12),
#                                 axis.text = element_text(size=10),
#                                 plot.title = element_text(size = 16)))
# )
# p.2
# p.3<-ggstatsplot::ggbetweenstats(
#   data = dplyr::filter(
#     temp,
#     variable %in% c("rmr")),
#   x = habitat_final,
#   y = value,
#   ggsignif.args = list(textsize = 2, tip_length = 0.01),
#   p.adjust.method = "bonferroni", 
#   type = "non-parametric",
#   palette = "default_jama",
#   package = "ggsci",
#   ylab=c("RMR"),
#   xlab=" ",
#   results.subtitle = FALSE,
#   title = "B)",
#   centrality.plotting = FALSE,
#   pairwise.comparisons = FALSE,
#   ggplot.component = list(theme(axis.title = element_text(size = 12),
#                                 axis.text = element_text(size=10),
#                                 plot.title = element_text(size = 16)))
# )
# p.3

# start here for lettering
p.4<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("uct")),
  x = habitat_final,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("UCT (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "A)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.4
p.5<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("lct")),
  x = habitat_final,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("LCT (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "B)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.5
p.6<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("tnz")),
  x = habitat_final,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("TNZ (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "C)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10), 
                                plot.title = element_text(size = 16)))
)
p.6
p.7<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("utl")),
  x = habitat_final,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("HTL (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "D)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.7
p.8<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("therm_sm")),
  x = habitat_final,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("TSM (ºC)"),
  xlab="Habitat Type",
  results.subtitle = FALSE,
  title = "E)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.8
p.9<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("warm_tol")),
  x = habitat_final,
  y = value,
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("WT (ºC)"),
  xlab="Habitat Type",
  results.subtitle = FALSE,
  title = "F)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.9

# Generate pdf figure with all plots
p.f.1<-combine_plots(
  list(p.4,p.5, p.6, p.7, p.8, p.9),
  plotgrid.args = list(nrow = 4),
  annotation.args = list(
  )
)
# ggsave("~/Downloads/Therm-Habitat-updated.jpeg",p.f.1,width =10, height = 15)
ggsave("~/Downloads/Therm-Habitat-updated.pdf",p.f.1,width =10, height = 15)





# setup for ggstatsplot with stratum
temp<-melt(dat,id.vars = "stratum",measure.vars = colnames(dat[,c(12:17,19:21)]))
set.seed(123)

# remove body mass, BMR, & RMR
# p.1<-ggstatsplot::ggbetweenstats(
#   data = dplyr::filter(
#     temp,
#     variable %in% c("mb")),
#   x = stratum,
#   y = value,
#   ggsignif.args = list(textsize = 2, tip_length = 0.01),
#   p.adjust.method = "bonferroni", 
#   type = "non-parametric",
#   palette = "default_jama",
#   package = "ggsci",
#   ylab=c("Body Mass (g)"),
#   xlab=" ",
#   results.subtitle = FALSE,
#   title = "A)",
#   centrality.plotting = FALSE,
#   pairwise.comparisons = FALSE,
#   ggplot.component = list(theme(axis.title = element_text(size = 12),
#                                 axis.text = element_text(size=10),
#                                 plot.title = element_text(size = 16)))
# )
# p.1
# p.2<-ggstatsplot::ggbetweenstats(
#   data = dplyr::filter(
#     temp,
#     variable %in% c("bmr")),
#   x = stratum,
#   y = value,
#   ggsignif.args = list(textsize = 2, tip_length = 0.01),
#   p.adjust.method = "bonferroni", 
#   type = "non-parametric",
#   palette = "default_jama",
#   package = "ggsci",
#   ylab=c("BMR"),
#   xlab=" ",
#   results.subtitle = FALSE,
#   title = "A)",
#   centrality.plotting = FALSE,
#   pairwise.comparisons = FALSE,
#   ggplot.component = list(theme(axis.title = element_text(size = 12),
#                                 axis.text = element_text(size=10),
#                                 plot.title = element_text(size = 16)))
# )
# p.2
# p.3<-ggstatsplot::ggbetweenstats(
#   data = dplyr::filter(
#     temp,
#     variable %in% c("rmr")),
#   x = stratum,
#   y = value,
#   ggsignif.args = list(textsize = 2, tip_length = 0.01),
#   p.adjust.method = "bonferroni", 
#   type = "non-parametric",
#   palette = "default_jama",
#   package = "ggsci",
#   ylab=c("RMR"),
#   xlab=" ",
#   results.subtitle = FALSE,
#   title = "B)",
#   centrality.plotting = FALSE,
#   pairwise.comparisons = FALSE,
#   ggplot.component = list(theme(axis.title = element_text(size = 12),
#                                 axis.text = element_text(size=10),
#                                 plot.title = element_text(size = 16)))
# )
# p.3

# start lettering here
p.4<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("uct")),
  x = stratum,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("UCT (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "A)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.4
p.5<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("lct")),
  x = stratum,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("LCT (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "B)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.5
p.6<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("tnz")),
  x = stratum,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("TNZ (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "C)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.6
p.7<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("utl")),
  x = stratum,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("HTL (ºC)"),
  xlab=" ",
  results.subtitle = FALSE,
  title = "D)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.7
p.8<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("therm_sm")),
  x = stratum,
  y = value,
  ggsignif.args = list(textsize = 2, tip_length = 0.01),
  p.adjust.method = "bonferroni", 
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("TSM (ºC)"),
  xlab="Vertical Stratum",
  results.subtitle = FALSE,
  title = "E)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.8
p.9<-ggstatsplot::ggbetweenstats(
  data = dplyr::filter(
    temp,
    variable %in% c("warm_tol")),
  x = stratum,
  y = value,
  type = "non-parametric",
  palette = "default_jama",
  package = "ggsci",
  ylab=c("WT (ºC)"),
  xlab="Vertical Stratum",
  results.subtitle = FALSE,
  title = "F)",
  centrality.plotting = FALSE,
  pairwise.comparisons = FALSE,
  ggplot.component = list(theme(axis.title = element_text(size = 12),
                                axis.text = element_text(size=10),
                                plot.title = element_text(size = 16)))
)
p.9

# Generate pdf figure with all plots
p.f.1<-combine_plots(
  list(p.4,p.5, p.6, p.7, p.8, p.9),
  plotgrid.args = list(nrow = 4),
  annotation.args = list(
  )
)
# ggsave("~/Downloads/Therm-Stratum-updated.jpeg",p.f.1,width =11, height = 15)
ggsave("~/Downloads/Therm-Stratum-updated.pdf",p.f.1,width =11, height = 15)

