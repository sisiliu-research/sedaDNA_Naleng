# by Sisi Liu: Sisi Liu: sisi.liu@awi.de | sisi.liu.research@gmail.com (permanent address)
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1

# rda and venn for terristrial vegetation

#== setting
dir.create("step03_assigned/03_rda/outTables")
dir.create("step03_assigned/03_rda/outFigures")
econame="terrestrial"
ngs="APMG-5-10-28-34-35"
ecogroup="vegetation"

#== load package
library(readr)
library(vegan)
library(ggplot2)

#== load data
# envi data
envi=read.csv(paste0("step03_assigned/03_rda/", econame,"_envi.csv"))
names(envi)
envi=envi[c("temperature", "glaciers_area_approx", "Permafrost_catchment_approx", "herbivory", "land_use")]
yak.sqrt2=sqrt(sqrt(envi["herbivory"]))
human.sqrt2=sqrt(sqrt(envi["land_use"]))
envi=cbind(envi[1:3], yak.sqrt2, human.sqrt2)
# plants
spe=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-TerrestrialPlants-N100_3000-p1f5.csv"), row.names = 1)
sel.spe=sqrt(sqrt(spe))

#== Check the Axis lengths of DCA1, use CCA if Axis lengths of DCA1 >4
decorana(sel.spe)
#== check single explained% for all predictors
nams <- names(envi)
# Single variance
meta.single.var=NULL
for (i in 1:length(nams)) {
  form.i <- formula(paste("sel.spe ~", paste(nams[i])))
  print(rda(form.i, data=envi))
  # explained
  rda.explained=rda(form.i, data=envi)$CCA$eig/rda(form.i, data=envi)$tot.chi*100
  #
  print(anova(rda(form.i, data=envi)))
  #
  single.var=anova(rda(form.i, data=envi))
  call.model=rda(form.i, data=envi)$call
  single.var$call.model=as.character(call.model)[2]
  single.var$rda.explained=rda.explained
  meta.single.var=rbind(meta.single.var, single.var)
}
write.csv(meta.single.var, paste0("step03_assigned/03_rda/outTables/anova-RDA-", ngs, "-", econame, "-", ecogroup, "-Single_variance-beforeVIF.csv"))

#== VIF check inner correlation
envi.cca=cca(sel.spe~., envi)
vif_df=vif.cca(envi.cca)
write.csv(vif_df, paste0("step03_assigned/03_rda/outTables/vif-", ngs, "-", econame, "-", ecogroup, "-r1.csv"))
# remove temperature
envi.cca=cca(sel.spe~., envi[, -1])
vif_df=vif.cca(envi.cca)
write.csv(vif_df, paste0("step03_assigned/03_rda/outTables/vif-", ngs, "-", econame, "-", ecogroup, "-r2.csv"))
#== rda using envi factors with vif < 3
envi.sig=envi[, c("glaciers_area_approx", "Permafrost_catchment_approx", "herbivory", "land_use")]
nams <- list(c("glaciers_area_approx", "Permafrost_catchment_approx"), "herbivory", "land_use")

# Single variance
meta.single.var=NULL
for (i in 1:length(nams)) {
  print(nams[[i]])
  #
  if (i == 1) {
    form.i <- formula(paste("sel.spe ~", paste(nams[[i]][1], "+", nams[[i]][2])))
  } else {
    form.i <- formula(paste("sel.spe ~", paste(nams[[i]])))
  }
  #
  print(rda(form.i, data=envi.sig))
  # explained
  rda.explained=rda(form.i, data=envi.sig)$CCA$eig/rda(form.i, data=envi.sig)$tot.chi*100
  #
  print(anova(rda(form.i, data=envi.sig)))
  #
  single.var=anova(rda(form.i, data=envi.sig))
  call.model=rda(form.i, data=envi.sig)$call
  single.var$call.model=as.character(call.model)[2]
  single.var$rda.explained=rda.explained
  meta.single.var=rbind(meta.single.var, single.var)
}
write.csv(meta.single.var, paste0("step03_assigned/03_rda/outTables/anova-RDA-", ngs, "-", econame, "-", ecogroup, "-afterVIF-Single_variance.csv"))

# total
mod1=rda(sel.spe ~., data=envi.sig)
total=round(rda(sel.spe ~., data=envi.sig)$CCA$eig/rda(sel.spe ~., data=envi.sig)$tot.chi*100, digits = 1)
total.anova=anova(rda(sel.spe ~., data=envi.sig))
total.anova$totalExp=total[1]
write.csv(total.anova, paste0("step03_assigned/03_rda/outTables/anova-RDA-", ngs, "-", econame, "-", ecogroup, "-afterVIF-total_variance.csv"))

#== Variation partitioning | Venn plot
# p < 0.5 (as single factor)
X1=envi.sig[1]
X2=envi.sig[2]
X3=envi.sig[3]
X4=envi.sig[4]
# Variation partitioning 
spe.part <- varpart(sel.spe, X1 + X2, X3, X4)
# Partition the variation in community composition
indfract=spe.part$part$indfract
write.csv(indfract, paste0("step03_assigned/03_rda/outTables/anova-RDA-", ngs, "-", econame, "-", ecogroup, "-afterVIF-indfract.csv"))
#
pdf(file = paste0("step03_assigned/03_rda/outFigures/Venn-", ngs, "-", econame, "-", ecogroup, "-single-P0.5.pdf"), width = 15, height = 10)
plot(spe.part, digits = 2, bg = c("red", "blue", "orange"), 
     Xnames = c("Cryosphere\n(Glacier and permafrost extent)", "Herbivory", "Land use"), 
     id.size = 0.8, cex = 0.8)
dev.off() 
#
# Test of fraction: Conditionial / Unique variance
meta.cond.var=NULL
for (i in 1:length(nams)) {
  print(i)
  if (i == 1) {
    form.i <- formula(paste("sel.spe ~", paste(nams[[i]][1], "+", nams[[i]][2]),
                            "+ Condition(", paste(nams[-i], collapse = " + "), ")"))
    
  } else {
    form.i <- formula(paste("sel.spe ~", paste(nams[[i]]),
                            "+ Condition(", paste(unlist(nams[-i]), collapse = " + "), ")"))
  }
  
  print(anova(rda(form.i, data=envi.sig)))
  cond.var=anova(rda(form.i, data=envi.sig))
  call.model=rda(form.i, data=envi.sig)$call
  cond.var$call.model=as.character(call.model)[2]
  meta.cond.var=rbind(meta.cond.var, cond.var)
}
write.csv(meta.cond.var, paste0("step03_assigned/03_rda/outTables/anova-RDA-", ngs, "-", econame, "-", ecogroup, "-afterVIF-Conditionial_variance.csv"))
#== end ==

