### PGLS models for therm phys project
# 04.20.21 - WJC
# Using tree script from Cameron L. Rutt
# Update August 2022 - WJC - 
#     Collapse stratum to 3 categories matching environmental data
#     Recalculate warm_tol & therm_sm

# Prepare for PGLS using a Maximum Clade Credibility tree built with TreeAnnotator (BEAST)

library(dplyr)
library(ape)
library(phytools)

setwd("C:/Users/hellr/Downloads")

birds<-read.csv("/Users/hellr/Downloads/thermal phys data_4.20.21_WJC.csv")

# collapse stratum to understory, midstory, and canopy
unique(birds$stratum)
birds[birds$stratum=="generalist","stratum"]<-"midstory"

##added in to collapse midstory and canopy species for analyses asked for on first revision
birds[birds$stratum=="canopy","stratum"]<-"midstory"

birds[
  birds$stratum=="terrestrial"|
    birds$stratum=="understory"|
    birds$stratum=="near-ground",
  "stratum"
  ]<-"understory"


#added in to collapse to two categories - forest vs open
birds[birds$habitat_final=="forest_edge","habitat_final"]<-"forest_interior"
birds[birds$habitat_final=="edge","habitat_final"]<-"open"


# Recalculate warm_tol for stratum max temp
birds[birds$stratum=="","warm_tol"]<-birds[birds$stratum=="","utl"]-30
birds[birds$stratum=="understory","warm_tol"]<-birds[birds$stratum=="understory","utl"]-28.2
birds[birds$stratum=="midstory","warm_tol"]<-birds[birds$stratum=="midstory","utl"]-29.3
birds[birds$stratum=="canopy","warm_tol"]<-birds[birds$stratum=="canopy","utl"]-29.5

# Recalculate therm_sm for stratum max temp
birds[birds$stratum=="","therm_sm"]<-birds[birds$stratum=="","uct"]-30
birds[birds$stratum=="understory","therm_sm"]<-birds[birds$stratum=="understory","uct"]-28.2
birds[birds$stratum=="midstory","therm_sm"]<-birds[birds$stratum=="midstory","uct"]-29.3
birds[birds$stratum=="canopy","therm_sm"]<-birds[birds$stratum=="canopy","uct"]-29.5


# load phylo tree
tree <- read.nexus("MCC_AllBirdsHackett1.tre")

birds$Scientific.name <- gsub(" ", "_", birds$Scientific.name)
birds <- rename(birds, Species = Scientific.name)

# See what birds are "missing" from the consensus tree 
# Remember that Glyph, Megarynchus, Pitangus, and Eucometis were all previously misspelled
setdiff(birds$Species, tree$tip.label) # soooo, we need to change the taxonomy on 13 species

birds$Species[birds$Species == "Cyanoloxia_cyanoides"] <- "Cyanocompsa_cyanoides"
birds$Species[birds$Species == "Momotus_subrufescens"] <- "Momotus_momota"
birds$Species[birds$Species == "Ceratopipra_mentalis"] <- "Pipra_mentalis"
birds$Species[birds$Species == "Gymnopithys_bicolor"] <- "Gymnopithys_leucaspis"
birds$Species[birds$Species == "Poliocrania_exsul"] <- "Myrmeciza_exsul"
birds$Species[birds$Species == "Schiffornis_stenorhyncha"] <- "Schiffornis_turdina"
birds$Species[birds$Species == "Juliamyia_julie"] <- "Damophila_julie"
birds$Species[birds$Species == "Pheugopedius_fasciatoventris"] <- "Thryothorus_fasciatoventris"
birds$Species[birds$Species == "Cantorchilus_leucotis"] <- "Thryothorus_leucotis"
birds$Species[birds$Species == "Troglodytes_musculus"] <- "Troglodytes_aedon"
birds$Species[birds$Species == "Thryophilus_rufalbus"] <- "Thryothorus_rufalbus"
birds$Species[birds$Species == "Pachysylvia_aurantiifrons"] <- "Hylophilus_aurantiifrons"
birds$Species[birds$Species == "Pachysylvia_decurtata"] <- "Hylophilus_decurtatus"

setdiff(birds$Species, tree$tip.label) # character(0) -- Okay, roger that 


#create dataframe of species with n>3
birdstnz<-birds[birds$n_tnz>2,]
birdsutl<-birds[birds$n_utl>2,]


birdstnz <- filter(birdstnz, rowSums(is.na(birdstnz)) != ncol(birdstnz))
birdsutl <- filter(birdsutl, rowSums(is.na(birdsutl)) != ncol(birdsutl))

# pruning the ~10,000 species tree down to our list of 90 species for this comparative analysis
phylo <- drop.tip(tree, setdiff(tree$tip.label, birdstnz$Species))
#plotting the tree
plot(ladderize(phylo), cex=0.6)
edgelabels(phylo$edge.length,cex=0.3)

library(caper)

# Make a comparative data object for PGLS
data <- comparative.data(tree, birds, names.col="Species")
phylo <- drop.tip(data$phy, setdiff(data$phy$tip.label, birds$Species))
plot(ladderize(phylo), cex=0.7)

# Run the PGLS analyses with maximum likelihood branch length transformations for lambda

##### start simple: how each metric relates to phylogeny ######
mb.phy<- pgls(mb~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(mb.phy)

utl.phy<- pgls(utl~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(utl.phy)

bmr.phy<- pgls(bmr~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(bmr.phy)

rmr.phy<- pgls(rmr~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(rmr.phy)

uct.phy<- pgls(uct~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(uct.phy)

lct.phy<- pgls(lct~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(lct.phy)

tnz.phy<- pgls(tnz~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(tnz.phy)

uct_active.phy<- pgls(uct_active~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(uct_active.phy)

therm_sm.phy<- pgls(therm_sm~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(therm_sm.phy)

warm_tol.phy<- pgls(warm_tol~1, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(warm_tol.phy)

##### then with just mb for each metric #####
bmr.phy<- pgls(bmr~mb, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(bmr.phy)

rmr.phy<- pgls(rmr~mb, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(rmr.phy)

hab.phy<- pgls(mb~habitat_final, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(hab.phy)

str.phy<- pgls(mb~stratum, comparative.data(tree,birds,"Species"),lambda = "ML")
summary(str.phy)

####test for differences in traits between temperate and tropical species
##### UTL - upper thermal limit ######

#create new data frame utl which includes species for which UTL was measured
utl<-birds[complete.cases(birds$utl),]

#pare down from bigger tree to include only species for which UTL was measured
utltree<-drop.tip(tree, setdiff(tree$tip.label, utl$Species))

#check to make sure all species are included
utltree$tip.label

#run the pgls function using maximum likelihood to compare UTL across groups
utl.slim <- pgls(utl~southern_limit, comparative.data(utltree,utl,"Species"),lambda = "ML")
summary(utl.slim)
coef(utl.slim)
plot(utl$utl~utl$southern_limit, ylab = "UTL", xlab = "Southern Limit")
abline(a=coef(utl.slim)[1],b=coef(utl.slim)[2])


utl.habi <- pgls(utl~habitat_final, comparative.data(utltree,utl,"Species"),lambda = "ML")
summary(utl.habi)

utlphylo<-phylANOVA(y=utl$utl,x=utl$habitat_final, tree = utltree,nsim=1000)



utl_strat<-utl[utl$stratum!="",]
keep <- utl_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>1) %>%
  pull(stratum)
utl_strat<-utl_strat[utl_strat$stratum %in% keep,]
utl_strat<-utl_strat[,c("Species","utl","stratum")]
utl_strat$stratum<-factor(utl_strat$stratum,levels = unique(utl_strat$stratum)[order(unique(utl_strat$stratum))])
utl.strat <- pgls(utl~stratum, comparative.data(utltree,utl_strat,"Species"),lambda = "ML")
summary(utl.strat)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, utl_strat$Species))
phylANOVA(y=utl_strat$utl,x=utl_strat$strat, tree = temp.tree,nsim=1000)


utl_diet<-utl[utl$diet!="",]
keep <- utl_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
utl_diet<-utl_diet[utl_diet$diet %in% keep,]
utl_diet$diet<-factor(utl_diet$diet,levels = unique(utl_diet$diet)[order(unique(utl_diet$diet))])
utl_diet<-utl_diet[,c("Species","utl","diet")]
utl.diet <- pgls(utl~diet, comparative.data(utltree,utl_diet,"Species"),lambda = "ML")
summary(utl.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, utl_diet$Species))
phylANOVA(y=utl_diet$utl,x=utl_diet$diet, tree = temp.tree,nsim=1000)



##### rmr - upper thermal limit ######

#create new data frame rmr which includes species for which rmr was measured
rmr<-birdstnz[complete.cases(birdstnz$rmr),]
rmr<-rmr[complete.cases(rmr$southern_limit),]

#pare down from bigger tree to include only species for which rmr was measured
rmrtree<-drop.tip(tree, setdiff(tree$tip.label, rmr$Species))

#check to make sure all species are included
rmrtree$tip.label

#run the pgls function using maximum likelihood to compare rmr across groups
rmr.slim <- pgls(rmr~southern_limit+mb, comparative.data(rmrtree,rmr,"Species"),lambda = "ML")
summary(rmr.slim)
coef(rmr.slim)
plot(rmr$rmr~rmr$southern_limit, ylab = "RMR", xlab = "Southern Limit")
abline(a=coef(rmr.slim)[1]+mean(rmr$mb*coef(rmr.slim)[3],na.rm = TRUE),b=coef(rmr.slim)[2])


library(evomap)
pGLS_ci<-gls.ci(rmr$rmr,rmr$southern_limit,vcv(rmrtree))
lines(pGLS_ci$CI.plot$X,pGLS_ci$CI.plot$Lower2.5,lty=2)
lines(pGLS_ci$CI.plot$X,pGLS_ci$CI.plot$Upper2.5,lty=2)

pGLS_pi<-gls.pi(rmr$rmr,rmr$southern_limit,vcv(rmrtree),1)
lines(pGLS_pi$PI.plot$X,pGLS_pi$PI.plot$Lower2.5,lty=2)
lines(pGLS_pi$PI.plot$X,pGLS_pi$PI.plot$Upper2.5,lty=2)

rmr.habi <- pgls(rmr~habitat_final+mb, comparative.data(rmrtree,rmr,"Species"),lambda = "ML")
summary(rmr.habi)
phylANOVA(y=rmr$rmr,x=rmr$habitat_final, tree = rmrtree,nsim=10000)


rmr_strat<-rmr[rmr$stratum!="",]
keep <- rmr_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>2) %>%
  pull(stratum)
rmr_strat<-rmr_strat[rmr_strat$stratum %in% keep,]
rmr_strat$stratum<-factor(rmr_strat$stratum,levels = unique(rmr_strat$stratum)[order(unique(rmr_strat$stratum))])
rmr_strat<-rmr_strat[,c("Species","rmr","stratum","mb")]
rmr.strat <- pgls(rmr~stratum+mb, comparative.data(rmrtree,rmr_strat,"Species"),lambda = "ML")
summary(rmr.strat)
plot(pgls.profile(rmr.strat))
# phyANOVA - new test
library(phytools)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, rmr_strat$Species))
phylANOVA(y=rmr_strat$rmr,x=rmr_strat$strat, tree = temp.tree,nsim=1000)


rmr_diet<-rmr[rmr$diet!="",]
keep <- rmr_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
rmr_diet<-rmr_diet[rmr_diet$diet %in% keep,]
rmr_diet$diet<-factor(rmr_diet$diet,levels = unique(rmr_diet$diet)[order(unique(rmr_diet$diet))])
rmr_diet<-rmr_diet[,c("Species","rmr","diet","mb")]
rmr.diet <- pgls(rmr~diet+mb, comparative.data(rmrtree,rmr_diet,"Species"),lambda = "ML")
summary(rmr.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, rmr_diet$Species))
phylANOVA(y=rmr_diet$rmr,x=rmr_diet$diet, tree = temp.tree,nsim=1000)



##### bmr - upper thermal limit ######

#create new data frame bmr which includes species for which bmr was measured
bmr<-birds[complete.cases(birds$bmr),]
bmr<-bmr[complete.cases(bmr$southern_limit),]

#pare down from bigger tree to include only species for which bmr was measured
bmrtree<-drop.tip(tree, setdiff(tree$tip.label, bmr$Species))

#check to make sure all species are included
bmrtree$tip.label

#run the pgls function using maximum likelihood to compare bmr across groups
bmr.slim <- pgls(bmr~southern_limit+mb, comparative.data(bmrtree,bmr,"Species"),lambda = "ML")
summary(bmr.slim)
coef(bmr.slim)
plot(bmr$bmr~bmr$southern_limit, ylab = "BMR", xlab = "Southern Limit")
abline(a=coef(bmr.slim)[1]+mean(bmr$mb*coef(bmr.slim)[3],na.rm = TRUE),b=coef(bmr.slim)[2])

library(evomap)
pGLS_ci<-gls.ci(bmr$bmr,bmr$southern_limit,vcv(bmrtree))
lines(pGLS_ci$CI.plot$X,pGLS_ci$CI.plot$Lower2.5,lty=2)
lines(pGLS_ci$CI.plot$X,pGLS_ci$CI.plot$Upper2.5,lty=2)


bmr.habi <- pgls(bmr~habitat_final+mb, comparative.data(bmrtree,bmr,"Species"),lambda = "ML")
summary(bmr.habi)
phylANOVA(y=bmr$bmr,x=bmr$habitat_final, tree = bmrtree,nsim=1000)

bmr_strat<-bmr[bmr$stratum!="",]
keep <- bmr_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>2) %>%
  pull(stratum)
bmr_strat<-bmr_strat[bmr_strat$stratum %in% keep,]
bmr_strat$stratum<-factor(bmr_strat$stratum,levels = unique(bmr_strat$stratum)[order(unique(bmr_strat$stratum))])
bmr_strat<-bmr_strat[,c("Species","bmr","stratum","mb")]
bmr.strat <- pgls(bmr~stratum+mb, comparative.data(bmrtree,bmr_strat,"Species"),lambda = "ML")
summary(bmr.strat)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, bmr_strat$Species))
phylANOVA(y=bmr_strat$bmr,x=bmr_strat$strat, tree = temp.tree,nsim=1000)

bmr_diet<-bmr[bmr$diet!="",]
keep <- bmr_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
bmr_diet<-bmr_diet[bmr_diet$diet %in% keep,]
bmr_diet$diet<-factor(bmr_diet$diet,levels = unique(bmr_diet$diet)[order(unique(bmr_diet$diet))])
bmr_diet<-bmr_diet[,c("Species","bmr","diet","mb")]
bmr.diet <- pgls(bmr~diet+mb, comparative.data(bmrtree,bmr_diet,"Species"),lambda = "ML")
summary(bmr.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, bmr_diet$Species))
phylANOVA(y=bmr_diet$bmr,x=bmr_diet$diet, tree = temp.tree,nsim=1000)



##### uct - upper thermal limit ######

#create new data frame uct which includes species for which uct was measured
uct<-birdstnz[complete.cases(birdstnz$uct),]

#pare down from bigger tree to include only species for which uct was measured
ucttree<-drop.tip(tree, setdiff(tree$tip.label, uct$Species))

#check to make sure all species are included
ucttree$tip.label

#run the pgls function using maximum likelihood to compare uct across groups
uct.slim <- pgls(uct~southern_limit, comparative.data(ucttree,uct,"Species"),lambda = "ML")
summary(uct.slim)
coef(uct.slim)
plot(uct$uct~uct$southern_limit, ylab = "UCT", xlab = "Southern Limit")
abline(a=coef(uct.slim)[1],b=coef(uct.slim)[2])

uct.habi <- pgls(uct~habitat_final, comparative.data(ucttree,uct,"Species"),lambda = "ML")
summary(uct.habi)
phylANOVA(y=uct$uct,x=uct$habitat_final, tree = ucttree,nsim=1000)

uct_strat<-uct[uct$stratum!="",]
keep <- uct_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>2) %>%
  pull(stratum)
uct_strat<-uct_strat[uct_strat$stratum %in% keep,]
uct_strat$stratum<-factor(uct_strat$stratum,levels = unique(uct_strat$stratum)[order(unique(uct_strat$stratum))])
uct_strat<-uct_strat[,c("Species","uct","stratum")]
uct.strat <- pgls(uct~stratum, comparative.data(ucttree,uct_strat,"Species"),lambda = "ML")
summary(uct.strat)
plot(uct_strat$stratum,uct_strat$uct)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, uct_strat$Species))
phylANOVA(y=uct_strat$uct,x=uct_strat$strat, tree = temp.tree,nsim=1000)



uct_diet<-uct[uct$diet!="",]
keep <- uct_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
uct_diet<-uct_diet[uct_diet$diet %in% keep,]
uct_diet$diet<-factor(uct_diet$diet,levels = unique(uct_diet$diet)[order(unique(uct_diet$diet))])
uct_diet<-uct_diet[,c("Species","uct","diet")]
uct.diet <- pgls(uct~diet, comparative.data(ucttree,uct_diet,"Species"),lambda = "ML")
summary(uct.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, uct_diet$Species))
phylANOVA(y=uct_diet$uct,x=uct_diet$diet, tree = temp.tree,nsim=1000)



##### lct - upper thermal limit ######

#create new data frame lct which includes species for which lct was measured
lct<-birdstnz[complete.cases(birdstnz$lct),]

#pare down from bigger tree to include only species for which lct was measured
lcttree<-drop.tip(tree, setdiff(tree$tip.label, lct$Species))

#check to make sure all species are included
lcttree$tip.label

#run the pgls function using maximum likelihood to compare lct across groups
lct.slim <- pgls(lct~southern_limit, comparative.data(lcttree,lct,"Species"),lambda = "ML")
summary(lct.slim)
coef(lct.slim)
plot(lct$lct~lct$southern_limit, ylab = "LCT", xlab = "Southern Limit")
abline(a=coef(lct.slim)[1],b=coef(lct.slim)[2])

lct.habi <- pgls(lct~habitat_final, comparative.data(lcttree,lct,"Species"),lambda = "ML")
summary(lct.habi)
phylANOVA(y=lct$lct,x=lct$habitat_final, tree = lcttree,nsim=1000)

lct_strat<-lct[lct$stratum!="",]
keep <- lct_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>1) %>%
  pull(stratum)
lct_strat<-lct_strat[lct_strat$stratum %in% keep,]
lct_strat$stratum<-factor(lct_strat$stratum,levels = unique(lct_strat$stratum)[order(unique(lct_strat$stratum))])
lct_strat<-lct_strat[,c("Species","lct","stratum")]
lct.strat <- pgls(lct~stratum, comparative.data(lcttree,lct_strat,"Species"),lambda = "ML")
summary(lct.strat)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, lct_strat$Species))
phylANOVA(y=lct_strat$lct,x=lct_strat$strat, tree = temp.tree,nsim=1000)


lct_diet<-lct[lct$diet!="",]
keep <- lct_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
lct_diet<-lct_diet[lct_diet$diet %in% keep,]
lct_diet$diet<-factor(lct_diet$diet,levels = unique(lct_diet$diet)[order(unique(lct_diet$diet))])
lct_diet<-lct_diet[,c("Species","lct","diet")]
lct.diet <- pgls(lct~diet, comparative.data(lcttree,lct_diet,"Species"),lambda = "ML")
summary(lct.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, lct_diet$Species))
phylANOVA(y=lct_diet$lct,x=lct_diet$diet, tree = temp.tree,nsim=1000)



##### tnz - thermal neutral zone ######

#create new data frame tnz which includes species for which tnz was measured
tnz<-birds[complete.cases(birds$tnz),]

#pare down from bigger tree to include only species for which tnz was measured
tnztree<-drop.tip(tree, setdiff(tree$tip.label, tnz$Species))

#check to make sure all species are included
tnztree$tip.label

#run the pgls function using maximum likelihood to compare tnz across groups
tnz.slim <- pgls(tnz~southern_limit, comparative.data(tnztree,tnz,"Species"),lambda = "ML")
summary(tnz.slim)
coef(tnz.slim)
plot(tnz$tnz~tnz$southern_limit, ylab = "TNZ", xlab = "Southern Limit")
abline(a=coef(tnz.slim)[1],b=coef(tnz.slim)[2])

tnz.habi <- pgls(tnz~habitat_final, comparative.data(tnztree,tnz,"Species"),lambda = "ML")
summary(tnz.habi)
phylANOVA(y=tnz$tnz,x=tnz$habitat_final, tree = tnztree,nsim=1000)

tnz_strat<-tnz[tnz$stratum!="",]
keep <- tnz_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>1) %>%
  pull(stratum)
tnz_strat<-tnz_strat[tnz_strat$stratum %in% keep,]
tnz_strat$stratum<-factor(tnz_strat$stratum,levels = unique(tnz_strat$stratum)[order(unique(tnz_strat$stratum))])
tnz_strat<-tnz_strat[,c("Species","tnz","stratum")]
tnz.strat <- pgls(tnz~stratum, comparative.data(tnztree,tnz_strat,"Species"),lambda = "ML")
summary(tnz.strat)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, tnz_strat$Species))
phylANOVA(y=tnz_strat$tnz,x=tnz_strat$strat, tree = temp.tree,nsim=1000)


tnz_diet<-tnz[tnz$diet!="",]
keep <- tnz_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
tnz_diet<-tnz_diet[tnz_diet$diet %in% keep,]
tnz_diet$diet<-factor(tnz_diet$diet,levels = unique(tnz_diet$diet)[order(unique(tnz_diet$diet))])
tnz_diet<-tnz_diet[,c("Species","tnz","diet")]
tnz.diet <- pgls(tnz~diet, comparative.data(tnztree,tnz_diet,"Species"),lambda = "ML")
summary(tnz.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, tnz_diet$Species))
phylANOVA(y=tnz_diet$tnz,x=tnz_diet$diet, tree = temp.tree,nsim=1000)



##### uct_active - upper thermal limit ######

#create new data frame uct_active which includes species for which uct_active was measured
uct_active<-birds[complete.cases(birds$uct_active),]

#pare down from bigger tree to include only species for which uct_active was measured
uct_activetree<-drop.tip(tree, setdiff(tree$tip.label, uct_active$Species))

#check to make sure all species are included
uct_activetree$tip.label

#run the pgls function using maximum likelihood to compare uct_active across groups
uct_active.slim <- pgls(uct_active~southern_limit, comparative.data(uct_activetree,uct_active,"Species"),lambda = "ML")
summary(uct_active.slim)
coef(uct_active.slim)
plot(uct_active$uct_active~uct_active$southern_limit, ylab = "UCT Active", xlab = "Southern Limit")
abline(a=coef(uct_active.slim)[1],b=coef(uct_active.slim)[2])

uct_active.habi <- pgls(uct_active~habitat_final, comparative.data(uct_activetree,uct_active,"Species"),lambda = "ML")
summary(uct_active.habi)
phylANOVA(y=uct_active$uct_active,x=uct_active$habitat_final, tree = uct_activetree,nsim=1000)

uct_active_strat<-uct_active[uct_active$stratum!="",]
keep <- uct_active_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>1) %>%
  pull(stratum)
uct_active_strat<-uct_active_strat[uct_active_strat$stratum %in% keep,]
uct_active_strat$stratum<-factor(uct_active_strat$stratum,levels = unique(uct_active_strat$stratum)[order(unique(uct_active_strat$stratum))])
uct_active_strat<-uct_active_strat[,c("Species","uct_active","stratum")]
uct_active.strat <- pgls(uct_active~stratum, comparative.data(uct_activetree,uct_active_strat,"Species"),lambda = "ML")
summary(uct_active.strat)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, uct_active_strat$Species))
phylANOVA(y=uct_active_strat$uct_active,x=uct_active_strat$strat, tree = temp.tree,nsim=1000)


uct_active_diet<-uct_active[uct_active$diet!="",]
keep <- uct_active_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
uct_active_diet<-uct_active_diet[uct_active_diet$diet %in% keep,]
uct_active_diet$diet<-factor(uct_active_diet$diet,levels = unique(uct_active_diet$diet)[order(unique(uct_active_diet$diet))])
uct_active_diet<-uct_active_diet[,c("Species","uct_active","diet")]
uct_active.diet <- pgls(uct_active~diet, comparative.data(uct_activetree,uct_active_diet,"Species"),lambda = "ML")
summary(uct_active.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, uct_active_diet$Species))
phylANOVA(y=uct_active_diet$uct_active,x=uct_active_diet$diet, tree = temp.tree,nsim=1000)



##### therm_sm - thermal safety margin ######

#create new data frame therm_sm which includes species for which therm_sm was measured
therm_sm<-birds[complete.cases(birds$therm_sm),]

#pare down from bigger tree to include only species for which therm_sm was measured
therm_smtree<-drop.tip(tree, setdiff(tree$tip.label, therm_sm$Species))

#check to make sure all species are included
therm_smtree$tip.label

#run the pgls function using maximum likelihood to compare therm_sm across groups
therm_sm.slim <- pgls(therm_sm~southern_limit, comparative.data(therm_smtree,therm_sm,"Species"),lambda = "ML")
summary(therm_sm.slim)
coef(therm_sm.slim)
plot(therm_sm$therm_sm~therm_sm$southern_limit, ylab = "Thermal Safety Margin", xlab = "Southern Limit")
abline(a=coef(therm_sm.slim)[1],b=coef(therm_sm.slim)[2])


therm_sm.habi <- pgls(therm_sm~habitat_final, comparative.data(therm_smtree,therm_sm,"Species"),lambda = "ML")
summary(therm_sm.habi)
phylANOVA(y=therm_sm$therm_sm,x=therm_sm$habitat_final, tree = therm_smtree,nsim=1000)

therm_sm_strat<-therm_sm[therm_sm$stratum!="",]
keep <- therm_sm_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>1) %>%
  pull(stratum)
therm_sm_strat<-therm_sm_strat[therm_sm_strat$stratum %in% keep,]
therm_sm_strat$stratum<-factor(therm_sm_strat$stratum,levels = unique(therm_sm_strat$stratum)[order(unique(therm_sm_strat$stratum))])
therm_sm_strat<-therm_sm_strat[,c("Species","therm_sm","stratum")]
therm_sm.strat <- pgls(therm_sm~stratum, comparative.data(therm_smtree,therm_sm_strat,"Species"),lambda = "ML")
summary(therm_sm.strat)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, therm_sm_strat$Species))
phylANOVA(y=therm_sm_strat$therm_sm,x=therm_sm_strat$strat, tree = temp.tree,nsim=1000)


therm_sm_diet<-therm_sm[therm_sm$diet!="",]
keep <- therm_sm_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
therm_sm_diet<-therm_sm_diet[therm_sm_diet$diet %in% keep,]
therm_sm_diet$diet<-factor(therm_sm_diet$diet,levels = unique(therm_sm_diet$diet)[order(unique(therm_sm_diet$diet))])
therm_sm_diet<-therm_sm_diet[,c("Species","therm_sm","diet")]
therm_sm.diet <- pgls(therm_sm~diet, comparative.data(therm_smtree,therm_sm_diet,"Species"),lambda = "ML")
summary(therm_sm.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, therm_sm_diet$Species))
phylANOVA(y=therm_sm_diet$therm_sm,x=therm_sm_diet$diet, tree = temp.tree,nsim=1000)




##### warm_tol - warming tolerance ######

#create new data frame warm_tol which includes species for which warm_tol was measured
warm_tol<-birds[complete.cases(birds$warm_tol),]

#pare down from bigger tree to include only species for which warm_tol was measured
warm_toltree<-drop.tip(tree, setdiff(tree$tip.label, warm_tol$Species))

#check to make sure all species are included
warm_toltree$tip.label

#run the pgls function using maximum likelihood to compare warm_tol across groups
warm_tol.slim <- pgls(warm_tol~southern_limit, comparative.data(warm_toltree,warm_tol,"Species"),lambda = "ML")
summary(warm_tol.slim)
coef(warm_tol.slim)
plot(warm_tol$warm_tol~warm_tol$southern_limit, ylab = "Warming Tolerance", xlab = "Southern Limit")
abline(a=coef(warm_tol.slim)[1],b=coef(warm_tol.slim)[2])


warm_tol.habi <- pgls(warm_tol~habitat_final, comparative.data(warm_toltree,warm_tol,"Species"),lambda = "ML")
summary(warm_tol.habi)
phylANOVA(y=warm_tol$warm_tol,x=warm_tol$habitat_final, tree = warm_toltree,nsim=1000)

warm_tol_strat<-warm_tol[warm_tol$stratum!="",]
keep <- warm_tol_strat %>%
  group_by(stratum) %>%
  summarise(no_rows = length(stratum)) %>%
  filter(no_rows>1) %>%
  pull(stratum)
warm_tol_strat<-warm_tol_strat[warm_tol_strat$stratum %in% keep,]
warm_tol_strat$stratum<-factor(warm_tol_strat$stratum,levels = unique(warm_tol_strat$stratum)[order(unique(warm_tol_strat$stratum))])
warm_tol_strat<-warm_tol_strat[,c("Species","warm_tol","stratum")]
warm_tol.strat <- pgls(warm_tol~stratum, comparative.data(warm_toltree,warm_tol_strat,"Species"),lambda = "ML")
summary(warm_tol.strat)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, warm_tol_strat$Species))
phylANOVA(y=warm_tol_strat$warm_tol,x=warm_tol_strat$strat, tree = temp.tree,nsim=1000)


warm_tol_diet<-warm_tol[warm_tol$diet!="",]
keep <- warm_tol_diet %>%
  group_by(diet) %>%
  summarise(no_rows = length(diet)) %>%
  filter(no_rows>1) %>%
  pull(diet)
warm_tol_diet<-warm_tol_diet[warm_tol_diet$diet %in% keep,]
warm_tol_diet$diet<-factor(warm_tol_diet$diet,levels = unique(warm_tol_diet$diet)[order(unique(warm_tol_diet$diet))])
warm_tol_diet<-warm_tol_diet[,c("Species","warm_tol","diet")]
warm_tol.diet <- pgls(warm_tol~diet, comparative.data(warm_toltree,warm_tol_diet,"Species"),lambda = "ML")
summary(warm_tol.diet)
temp.tree<-drop.tip(tree, setdiff(tree$tip.label, warm_tol_diet$Species))
phylANOVA(y=warm_tol_diet$warm_tol,x=warm_tol_diet$diet, tree = temp.tree,nsim=1000)



##### global model ######
#run the pgls function using maximum likelihood to compare bmr across groups
##### rmr - upper thermal limit ######

#create new data frame tnz which includes species for which tnz was measured
tnz<-birds[complete.cases(birds$tnz),]

#pare down from bigger tree to include only species for which tnz was measured
tnztree<-drop.tip(tree, setdiff(tree$tip.label, tnz$Species))

#check to make sure all species are included
tnztree$tip.label
tnz.glob<-tnz[,c("Species","tnz","stratum","habitat_final","southern_limit","diet")]
tnz.glob<-tnz.glob[tnz.glob$diet!="piscivore",]
tnz.glob<-tnz.glob[tnz.glob$habitat_final!="edge"&tnz.glob$habitat_final!="open",]
tnz.glob$diet<-factor(tnz.glob$diet,levels = unique(tnz.glob$diet)[order(unique(tnz.glob$diet))])
tnz.glob$habitat_final<-factor(tnz.glob$habitat_final,levels = unique(tnz.glob$habitat_final)[order(unique(tnz.glob$habitat_final))])
tnz.glob$stratum<-factor(tnz.glob$stratum,levels = unique(tnz.glob$stratum)[order(unique(tnz.glob$stratum))])
tnz.glob_pgls <- pgls(tnz~southern_limit+habitat_final+stratum+diet, comparative.data(tnztree,tnz.glob,"Species"),lambda = "ML")
summary(tnz.glob_pgls)

