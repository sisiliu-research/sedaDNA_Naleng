# by Sisi Liu: sisi.liu@awi.de; sisi.liu.research@gmail.com
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1

#==== README ====
# This script for rda and variation partitioning for terrestrial mammalian

# Please read: Legendre and Legendre, 2012, Numerical ecology | 11.1 Redundancy analysis (RDA) to deeply understand:
#     why species data (response variables) Y and exploratory data (envi data) X should be centred data?
#     why X should be standardized if it includes multiple unit dimensions?

# further read if you are interested in "What is the difference between variation partitioning and variance partitioning (Lepš & Šmilauer 2003, 2014, Multivariate Analysis of Ecological Data using CANOCO)?

#==== Load package ====
library(readr)
library(vegan)
library(ggplot2)

#==== Working path and folder ====
dir.create("11_ngsLCA_L30_post/04_rda/01Tables")
dir.create("11_ngsLCA_L30_post/04_rda/01Figures")

#==== Working NGS ====
ngs="APMG-5-10-28-34-35"
orgs="MAM"
neta=2450 # determined by network analysis
eco="terrestrial"
trans="rlog"

#==== Load envi data, Source: sheet = ENVI in Dataset 1====
envi=read.csv(paste0("11_ngsLCA_L30_post/04_rda/ENVI.csv"), row.names = 1)[c("Temperature", "Glaciers_area", "Permafrost_catchment", "PC1", "Land_use")]

# Standardize quantitative environmental data (non-centerd data with different unit dimensions should be standardized)
envi$Temperature <- decostand(envi$Temperature, method = "standardize")
envi$Glaciers_area <- decostand(envi$Glaciers_area, method = "standardize")
envi$Permafrost_catchment <- decostand(envi$Permafrost_catchment, method = "standardize")
envi$PC1 <- decostand(envi$PC1, method = "standardize")
envi$Land_use <- decostand(envi$Land_use, method = "standardize")

#==== Load rlog-count NGS data ====
if (is.null(neta)) {
  t.spe=read.csv(paste0("11_ngsLCA_L30_post/02_nomolize/00Tables/01_", orgs, "_", trans, ".csv"), row.names = 1)
} else {
  orgs=paste0(orgs, "-", neta)
  t.spe=read.csv(paste0("11_ngsLCA_L30_post/02_nomolize/00Tables/06_", orgs, "_", trans, ".csv"), row.names = 1)
}
#==== CCA or RDA ==== 
#Axis lengths of DCA1, use CCA if Axis lengths of DCA1 > 4
# decorana does work for negative values (e.g., results of rlog and vst, positive values are approximately equal to select common taxa)
if(!(trans %in% c("rlog", "rclr", "vst"))) {
  decorana(t.spe)
} else {
  pov_cols <- t.spe %>% 
    mutate_all(~ ifelse(. < 0, 0, .)) %>% 
    select_if(~ sum(.) > 0)
  decorana(pov_cols)
}
#==== Find sig. envi ==== 
#check explained % for each predictor
nams <- names(envi)
# Single variance
meta.single.var=NULL
for (i in 1:length(nams)) {
  form.i <- formula(paste("t.spe ~", paste(nams[i])))
  print(rda(form.i, data=envi))
  modi=rda(form.i, data=envi)
  # anova table and adjust R^2 %
  an.modi=data.frame(anova(modi), adj.r.squared = round(RsquareAdj(modi)$adj.r.squared*100, 3), call=as.character(modi$call)[2])
  meta.single.var=rbind(meta.single.var, an.modi)
}
write.csv(meta.single.var, paste0("11_ngsLCA_L30_post/04_rda/00Tables/", orgs, "_", trans, "_sRDA_", ngs, "-", eco, ".csv"))

#==== Do VIF ====
envi.rda=rda(t.spe~., envi)
vif_df=round(vif.cca(envi.rda))
write.csv(vif_df, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_", ngs, "-", eco,  "-r1.csv"))

# remove temperature
envi.rda=rda(t.spe~., envi[, -3])
vif_df=round(vif.cca(envi.rda))
write.csv(vif_df, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_", ngs, "-", eco,  "-r2.csv"))

# remove PC1
envi.rda=rda(t.spe~., envi[, c(1,2,5)])
vif_df=round(vif.cca(envi.rda))
write.csv(vif_df, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_", ngs, "-", eco,  "-r3.csv"))

#==== Joint envi with rda using envi with vif <= 3 ====
envi.sig=envi[, c("Glaciers_area", "Permafrost_catchment", "Land_use")]
nams <- list(c("Glaciers_area", "Permafrost_catchment"), "Land_use")
# Single variance
meta.single.var=NULL
for (i in 1:length(nams)) {
  print(nams[[i]])
  #
  if (i == 1) {
    form.i <- formula(paste("t.spe ~", paste(nams[[i]][1], "+", nams[[i]][2])))
  } else {
    form.i <- formula(paste("t.spe ~", paste(nams[[i]])))
  }
  #
  print(rda(form.i, data=envi.sig))
  modi=rda(form.i, data=envi)
  # anova table and adjust R^2 %
  an.modi=data.frame(anova(modi), adj.r.squared = round(RsquareAdj(modi)$adj.r.squared*100, 1), call=as.character(modi$call)[2])
  meta.single.var=rbind(meta.single.var, an.modi)
}
write.csv(meta.single.var, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_pjRDA_", ngs, "-", eco, ".csv"))

#==== Total explained % ====
mod1=rda(t.spe ~., data=envi.sig)
#test for adding all environmental variables
an.mod=anova(mod1)
an.adjR2=data.frame(adj.r.squared=round(RsquareAdj(mod1)$adj.r.squared*100), call=as.character(mod1$call)[2])
#tests for individual terms
an.axis=anova(mod1, by='axis')
an.margin=anova(mod1, by='margin')
# out
write.csv(an.mod, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_jRDA_", ngs, "-", eco, ".csv"))
write.csv(an.adjR2, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_jRDA_adjR2_", ngs, "-", eco, ".csv"))
write.csv(an.axis, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_jRDA_axis_", ngs, "-", eco, ".csv"))
write.csv(an.margin, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_vif_jRDA_margin_", ngs, "-", eco, ".csv"))

#==== Variation partitioning | Venn plot ====
# p < 0.5 (as single factor)
X1=envi.sig[1]
X2=envi.sig[2]
X3=envi.sig[3]
# Variation partitioning 
spe.part <- varpart(t.spe, X1 + X2, X3)
# Partition the variation in community composition
indfract=spe.part$part$indfract
write.csv(indfract, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_indfract_", ngs, "-", eco, ".csv"))
#
pdf(file = paste0("11_ngsLCA_L30_post/04_rda/01Figures/", orgs, "_", trans, "_venn_", ngs, "-", eco,  ".pdf"), width = 15, height = 10)
plot(spe.part, digits = 2, bg = c("red", "blue", "orange"), 
     Xnames = c("Cryosphere\n(Glacier and permafrost extent)", "Land use"), 
     id.size = 0.8, cex = 0.8)
dev.off() 

# Test of fraction: Conditionial / Unique variance
meta.cond.var=NULL
for (i in 1:length(nams)) {
  print(i)
  if (i == 1) {
    form.i <- formula(paste("t.spe ~", paste(nams[[i]][1], "+", nams[[i]][2]),
                            "+ Condition(", paste(nams[-i], collapse = " + "), ")"))
    
  } else {
    form.i <- formula(paste("t.spe ~", paste(nams[[i]]),
                            "+ Condition(", paste(unlist(nams[-i]), collapse = " + "), ")"))
  }
  
  print(anova(rda(form.i, data=envi.sig)))
  cond.var=anova(rda(form.i, data=envi.sig))
  call.model=rda(form.i, data=envi.sig)$call
  cond.var$call.model=as.character(call.model)[2]
  meta.cond.var=rbind(meta.cond.var, cond.var)
}
write.csv(meta.cond.var, paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_venn_", ngs, "-", eco, ".csv"))
# save image to submit
save.image(file = paste0("11_ngsLCA_L30_post/04_rda/01Tables/", orgs, "_", trans, "_VarPartitioning_", ngs, "_", eco, ".RData"))
#==== END ====










