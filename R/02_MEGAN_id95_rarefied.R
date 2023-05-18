# by Sisi Liu: sisi.liu@awi.de | sisi.liu.research@gmail.com (permanent address)
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1

#== rarefy 

library(readr)
library(dplyr)
setwd("~/megan095")

#== create folders
dir.create("step03_assigned/02_rarefy/")
dir.create("step03_assigned/02_rarefy/outTables")
dir.create("step03_assigned/02_rarefy/outFigures")

#== setting
ngs="APMG-5-10-28-34-35"
ecogroup="TerrestrialPlants" # same procedures for other ecogroup

#== load
indf=read.csv(paste0("step03_assigned/outTables/", ngs,"-idc95-cellRanks-c2-clean-freq2-family-Viridiplantae-FaGeSpOthInTP-land.csv"), row.names = 1)
#indf=read.csv(paste0("step03_assigned/outTables/", ngs, "-idc95-cellRanks-c2-clean-freq2-family-MammaliaExHominidae-FaGeSpOthInTP-land.csv"), row.names = 1)
#indf=read.csv(paste0("step03_assigned/outTables/", ngs, "-idc95-cellRanks-c2-clean-freq2-family-aquaticGroups.csv"), row.names = 1)

#== view
length(unique(indf[indf$rank == "Species", "taxon"])) #96
length(unique(indf[indf$rank == "Genus", "taxon"])) #118
length(unique(indf[indf$rank == "Family", "taxon"])) #68

#== select cols and aggregate
subdf=indf[c("age", "taxid", "taxonReads", "taxon", "genus", "family", "rank")]
# sum up to genus level
ge=subdf[subdf$rank %in% c("Genus", "Species", "subgenus"), ]
agg_ge=aggregate(.~age+genus, ge[c("age", "genus", "taxonReads")], FUN=sum)
# attach family name
falist=unique(ge[c("genus", "family")])
agg_ge=left_join(agg_ge, falist, by = "genus")
# "subtribe", "subfamily", "tribe" to Family
fa=subdf[subdf$rank %in% c("Family", "subtribe", "subfamily", "tribe"), ]
agg_fa=aggregate(.~age+family, fa[c("age", "family", "taxonReads")], FUN=sum)
agg_fa$taxon=agg_fa$family
agg_fa=agg_fa[c("age", "taxon", "taxonReads", "family" )]
# rename
names(agg_ge)=c("age", "taxa", "reads", "family")
names(agg_fa)=c("age", "taxa", "reads", "family")
# merge
fgdf=rbind(agg_ge, agg_fa)
fgdf=fgdf[order(-fgdf$age, fgdf$taxa), ]
rownames(fgdf)=seq(1, dim(fgdf)[1], 1)
# long to wide
df=reshape(fgdf[, c("taxa", "family", "age", "reads")], idvar = c("taxa", "family"), timevar = "age", direction = "wide")
df[is.na(df)]=0
# change colname
names(df)=sub(".*reads.", "", colnames(df)) 

#==  determinate base count for resampling
library(vegan)
# age
age=unique(subdf$age)
# check sample size for rarefying community
rownames(df)=df[, 1]
df=df[-c(1,2)]
df=as.data.frame(t(df))
rownames(df)=age
# rarefy curve based on minimal count
S <- specnumber(df)
(raremax <- min(rowSums(df)))
Srare <- rarefy(df, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
pdf(file = paste0("step03_assigned/02_rarefy/outFigures/", ngs, "-", ecogroup, "_rarecurve_min.pdf"))
rarecurve(df, step = 20, sample = raremax, col = "blue", cex = 0.8, label = T)
title(paste0("Rare curve based on min rowSum: ", raremax, " for ", ecogroup))
dev.off()

# sequence
basecount=seq(1500, 3500, 100)
#basecount=seq(100, 250, 50) # mammalian
#basecount=seq(500, 2500, 100) # aquatic 
for (i in basecount) {
  # check if higher base count have higher species coverage
  pdf(file = paste0("step03_assigned/02_rarefy/outFigures/",  ngs, "-", ecogroup, "_rarecurve_", i, ".pdf"))
  rarecurve(df, step = 40, sample = i, col = "blue", cex = 0.8, label = T) # step = 20 (default) for mammalian and aquatic due to fewer total reads
  title(paste0("Rare curve based on median rowSum: ", i, " for ", ecogroup))
  dev.off()
}
#== select a base count
basecount=3000
#basecount=150 # mammalian
#basecount=2000 # aquatic

#== rarefy
# long to wide
df=reshape(fgdf[, c("taxa", "family", "age", "reads")], idvar = c("taxa", "family"), timevar = "age", direction = "wide")
df[is.na(df)]=0
specseq=df
names(specseq)[1]="Species"
names(specseq)[2]="Family"
# make Sample//species names unique
specseq$SampleUnique=make.unique(specseq$Species)
# add familyname to speciesname for later sorting to families
specseq$SampleUnique=paste0(specseq$Family,"_",specseq$SampleUnique)
str(specseq)
# get the sample
samplespresent=names(select(specseq, contains("reads.")))
# get the minimal total count per sample and total number of species per sample
metastat=NULL
for(cli in samplespresent) {
  # find minimum sample number == rarefaction num
  colnamesrepeats=names(specseq)[grep(cli,names(specseq))]
  subdf=specseq[,c("SampleUnique",colnamesrepeats)]
  # exclude zero species
  if(dim(subdf)[2]>2)
  {
    subdf=subdf[apply(subdf[,2:dim(subdf)[2]],1,sum)>0,]
    rarefytocounts=min(apply(subdf[,2:dim(subdf)[2]], 2, sum))
  } else
  {
    subdf=subdf[subdf[,2]>0,]
    rarefytocounts=sum(subdf[,2])
  }
  metastat=rbind(metastat, data.frame(Samples=cli, Rarefaction=rarefytocounts, TotalSpeciesNum=dim(subdf)[1]))
}
metastat
nsampleff=basecount
nrepeats=100 # the randomly sub-sampling will run 100 times
# generate taxa pool per sample
genrare=list()
allrepeatcolnames=samplespresent
for(yrcoli in allrepeatcolnames){
  print(yrcoli)
  allspec=specseq$SampleUnique[specseq[,yrcoli]>0]
  allspec_counts=specseq[,yrcoli][specseq[,yrcoli]>0]
  if(length(allspec)>0)
  {
    sampleeffort=list()
    for(nsampleeffi in nsampleff)
    {
      repeatsample=list()
      for(repi in 1:nrepeats)
      {
        repeatsample[[repi]]=sample(allspec,nsampleeffi,replace=TRUE, prob=allspec_counts/sum(allspec_counts))
      }
      sampleeffort[[which(nsampleff==nsampleeffi)]]=repeatsample
    }
    genrare[[which(allrepeatcolnames==yrcoli)]]=list(allspec,sampleeffort)
  } else
  {
    print(paste0("ERROR, no species to sample! ",yrcoli))
  }
}
names(genrare)=allrepeatcolnames
length(genrare)

# Rarefaction on Family level: this output is the absolute species richness on family level
famorignames=unique(specseq$Family)
famorignames[is.na(famorignames)]="NA"
familylevels=famorignames
totspec=NULL
totfam=NULL
for(li in 1:length(genrare)){
  print(li)
  for(li2 in 1:length(genrare[[li]][[2]]))
  {
    
    spectot=NULL
    spectot4fam=NULL
    for(repi in 1:nrepeats)
    {
      pei=unique(genrare[[li]][[2]][[li2]][[repi]])
      spectot=c(spectot,length(pei))
      spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split="_",pei),function(x)return(x[1]))),levels=familylevels)))
    }
    totspec=rbind(totspec, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),Nspecies=spectot))
    totfam=rbind(totfam, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
  }
}
str(totspec)
str(totfam)
# aggregate mean and confidence interval 95% for families
famdf=NULL
for(sampleidi in unique(totfam$T)){
  sampleidi
  totspecsub=totspec[totspec$T==sampleidi,"Nspecies"]
  famdf=rbind(famdf,data.frame(Sample=sampleidi, Family="ALL", NSpecMean=mean(totspecsub), NSpecCI95=sd(totspecsub)*1.96)) 
  # "ALL" is the rarefied total species richness/sample: for correlation coefficient calculation and plotting in Fig. 2
  for(familyi in names(totfam)[3:dim(totfam)[2]])
  {
    totfamsub=totfam[totfam$T==sampleidi,familyi]	
    famdf=rbind(famdf,data.frame(Sample=sampleidi, Family=familyi, NSpecMean=mean(totfamsub), NSpecCI95=sd(totfamsub)*1.96))
  }
}
str(famdf)

# Rarefaction on species level: to get the rarefied count per sequence per sample
specorignames=unlist(lapply(strsplit(split="_", unique(specseq$SampleUnique)),function(x)return(x[2])))
specorignames[is.na(specorignames)]="NA"
specieslevels=specorignames
taxa100=NULL
for(li in 1:length(genrare)){
  print(li)
  for(li2 in 1:length(genrare[[li]][[2]]))
  {
    spectot4fam=NULL
    for(repi in 1:nrepeats)
    {
      pei=genrare[[li]][[2]][[li2]][[repi]]
      spectot4fam=rbind(spectot4fam,table(factor(unlist(lapply(strsplit(split="_",pei),function(x)return(x[2]))),levels=specieslevels)))
    }
    taxa100=rbind(taxa100, data.frame(T=names(genrare)[li],SampleEff=length(genrare[[li]][[2]][[li2]][[repi]]),spectot4fam))
  }
}
str(taxa100)
rowSums(taxa100[,-(1:2)]) # Make sure the total counts in each sample equals to nsampleff

# aggregate  mean and confidence interval 95%
specdfmean=NULL
for(sampleidi in unique(taxa100$T))
{
  print(sampleidi)
  totspecsub=totspec[totspec$T==sampleidi,"Nspecies"]
  specdfmean=rbind(specdfmean,data.frame(Sample=sampleidi, Family="ALL", NSpecMean=mean(totspecsub), NSpecCI95=sd(totspecsub)*1.96))
  
  for(familyi in names(taxa100)[3:dim(taxa100)[2]])
  {
    totfamsub=taxa100[taxa100$T==sampleidi,familyi]	
    specdfmean=rbind(specdfmean,data.frame(Sample=sampleidi, Family=familyi, NSpecMean=mean(totfamsub), NSpecCI95=sd(totfamsub)*1.96))
  }
}
str(specdfmean)
#
raredf=specdfmean[specdfmean$Family !="ALL", ]
names(raredf)[2]=c("Species")
raredf=left_join(raredf, specseq[c(1,2)], by = "Species")
raredf=raredf[c("Species", "Family", "Sample", "NSpecMean", "NSpecCI95")]
names(raredf)=c("taxa", "family", "sample", "NSpecMean", "NSpecCI95")
write.csv(raredf, paste0("step03_assigned/02_rarefy/outTables/", ngs,"-hops-TP-Rarefy_TaxonLevel_N", nrepeats, "_", nsampleff, "_", ecogroup, ".csv"))
# long to wide
raredf_wide=reshape(raredf[1:4], idvar = c("taxa", "family"), timevar = "sample", direction = "wide")  
# change colname
names(raredf_wide)=sub(".*reads.", "", colnames(raredf_wide)) 
# %
raredf_wide_prob=cbind(raredf_wide[c(1,2)], as.data.frame(prop.table(as.matrix(raredf_wide[-c(1,2)]), margin = 2)*100))
colSums(raredf_wide_prob[-c(1,2)])
# save
write.csv(raredf_wide, paste0("step03_assigned/02_rarefy/outTables/", ngs,"-hops-TP-Rarefy_TaxonLevel_N", nrepeats, "_", nsampleff, "_", ecogroup, "_wide.csv"))
write.csv(raredf_wide_prob, paste0("step03_assigned/02_rarefy/outTables/", ngs,"-hops-TP-Rarefy_TaxonLevel_N", nrepeats, "_", nsampleff, "_", ecogroup, "_wide_prop.csv"))

#== select common taxa for profile ploting, PCA, RDA, and Network analysis
# max >= 1% and times >= 5
opt="p1f5"
df=raredf_wide_prob[-c(2)]
rownames(df)=df[, 1]
df[, 1]=NULL
df=as.data.frame(t(df))
#
max.abb <- apply(df, 2, max)
n.occ <- colSums(df > 0)
spp.want <- which(max.abb >= 1 & n.occ >= 5)
df.want=df[, spp.want] # 69
write.csv(df.want, paste0("step03_assigned/02_rarefy/outTables/", ngs, "-", ecogroup, "-N", nrepeats, "_", nsampleff, "-", opt, ".csv"))
#== rarefaction end == 

#== PCA analysis only for terrestrial plants (Fig. 2 J)
pca=rda(sqrt(sqrt(df.want)))
plot(pca)
pca.scores=scores(pca, scaling = 3)
pca.sites=pca.scores$site
write.csv(df.want, paste0("step03_assigned/02_rarefy/outTables/", ngs, "-", ecogroup, "-N", nrepeats, "_", nsampleff, "-", opt, "-PCA.csv"))
#== PCA end ==

