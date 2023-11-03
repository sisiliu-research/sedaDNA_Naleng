# by Sisi Liu: sisi.liu@awi.de; sisi.liu.research@gmail.com
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 14.1

#==== Load package ====
library(readr)
library(dplyr)
library(DESeq2)
library(DEP)
library(ggplot2)

#==== Working path and folders ====
setwd("~/04_NLshotgun/bowtie2")
# generate or reproduce?
reproduce="FALSE"
if (reproduce) {
  load("11_ngsLCA_L30_post/02_nomolize/nomolizeEnvi.Rdata")
} else {
  dir.create("11_ngsLCA_L30_post/02_nomolize")
  dir.create("11_ngsLCA_L30_post/02_nomolize/01Tables")
  dir.create("11_ngsLCA_L30_post/02_nomolize/01Figures")
}

#==== Which NGS data ====
# metagenomics (shotgun sequencing)
ngs="APMG-5-10-28-34-35"
orgs="PLANT"

#==== Load NGS data, Source: Dataset 1 ====
# taxa, family, samples
# raw count
if (orgs == "PLANT") {
  df2=read.csv(paste0("11_ngsLCA_L30_post/00_data/Tables/RAW_COUNT_", orgs, ".csv"), check.names = F, row.names = 2)
  df2=df2[-c(1,2)]
  neta=2450
} else if (orgs == "MAM") {
  df2=read.csv(paste0("11_ngsLCA_L30_post/00_data/Tables/RAW_COUNT_", orgs, ".csv"), check.names = F, row.names = 2)
  df2=df2[-c(1,2)]
  neta=2450
} else {
  df2=read.csv(paste0("11_ngsLCA_L30_post/00_data/Tables/RAW_COUNT_", orgs, ".csv"), check.names = F, row.names = 2)
  df2=df2[-c(1,2)]
  neta=50
  
}
# species x samples
spe=df2

if(TRUE) {
  #==== Profile transformation ====
  # the name varies in fields, e.g., profile transformation in Ecology, TSS in Biology
  spe.prop=as.data.frame(prop.table(as.matrix(spe), margin = 2)*100)
  colSums(spe.prop)
  #write.csv(spe.prop, paste0("11_ngsLCA_L30_post/02_nomolize/01Tables/01_", orgs,"-prop.csv"))
  
  #==== Mean-variance relationship ====
  countData <- as.matrix(spe)
  condition <- factor(seq_len(ncol(spe)))
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
  #
  trans="raw"
  meanSdPlot(dds)
  #ggsave(file = paste0("11_ngsLCA_L30_post/02_nomolize/01Figures/", orgs, "_", trans, "_meanvar_DEP.pdf"), width = 4, height = 3)
  #
  trans="vst"
  vsd=varianceStabilizingTransformation(dds, blind = TRUE, fitType = "parametric")
  meanSdPlot(vsd)
  #ggsave(file = paste0("11_ngsLCA_L30_post/02_nomolize/01Figures/", orgs, "_", trans, "_meanvar_DEP.pdf"), width = 4, height = 3)
  #
  trans="rlog"
  rld=rlog(dds)
  meanSdPlot(rld)
  #ggsave(file = paste0("11_ngsLCA_L30_post/02_nomolize/01Figures/", orgs, "_", trans, "_meanvar_DEP.pdf"), width = 4, height = 3)
  
  #==== save ====
  #write.csv(t(assay(vsd)), paste0("11_ngsLCA_L30_post/02_nomolize/01Tables/01_", orgs, "_vst.csv"))
  #write.csv(t(assay(rld)), paste0("11_ngsLCA_L30_post/02_nomolize/01Tables/01_", orgs, "_rlog.csv"))
  
  #==== PCA ====
  pca=rda(t(assay(rld)))
  plot(pca, display="site")
  pca.scores.site=data.frame(scores(pca, choices=c(1,2))[["sites"]])
  # Get the percentage of variance explained by each PC
  pca_var_explained <- pca$CA$eig / sum(pca$CA$eig) * 100
  #write.csv(pca.scores.site, paste0("11_ngsLCA_L30_post/02_nomolize/01Tables/02_", orgs, "_rlog_PC12.csv"))
  #write.csv(pca_var_explained, paste0("11_ngsLCA_L30_post/02_nomolize/01Tables/03_", orgs, "_rlog_PCs_explained.csv"))
  
  #=== Subset based on common taxa ====
  tmp=as.data.frame(t(spe))
  tax=colnames(tmp[, colSums(tmp) >= neta])
  tmp=t(spe.prop)
  tmp=tmp[, colnames(tmp) %in% tax]
  #write.csv(tmp, paste0("11_ngsLCA_L30_post/02_nomolize/01Tables/04_", orgs,"-", neta, "_prop.csv"))
  
  tmp=t(assay(rld))
  tmp=tmp[, colnames(tmp) %in% tax]
  #write.csv(tmp, paste0("11_ngsLCA_L30_post/02_nomolize/01Tables/05_", orgs,"-", neta, "_rlog.csv"))
}
#==== END ====
