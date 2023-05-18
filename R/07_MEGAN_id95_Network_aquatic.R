# by Sisi Liu: Sisi Liu: sisi.liu@awi.de | sisi.liu.research@gmail.com (permanent address)
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1

# model-based network 

#== setting
dir.create("step03_assigned/04_network/")
dir.create("step03_assigned/04_network/outTables_aquatic")
dir.create("step03_assigned/04_network/outFigures_aquatic")
ngs="APMG-5-10-28-34-35"
econame="aquatic"

#== packages
library(ecoCopula)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(readr)
library(reshape2)
library(igraph)
library(ggraph)

#== partial graph
graph_from_partial<-function(partial){
  diag(partial)<-0
  vertex <- data.frame(name = colnames(partial))
  Edges=igraph::as_data_frame(igraph::graph_from_adjacency_matrix(partial,mode="undirected",weighted = TRUE))
  if(ncol(Edges)==3){
    colnames(Edges)[3]<- "partcor"
  }
  
  igraph::graph_from_data_frame(Edges, directed=FALSE, vertices=vertex)
}
# ecoCopula works by first estimating a model, and then using the fitted model to estimate direct associations.
# For binary and count data we use the manyglm function from the mvabund package to estimate the model.

#== load data
# envi
envi=read.csv(paste0("step03_assigned/03_rda/", econame, "_envi.csv"))
names(envi)
envi_df=envi[c("age", "glacier_mass", "kramer2010")]
human.sqrt2=sqrt(sqrt(envi_df["kramer2010"]))
envi_df=cbind(envi_df[1:2], human.sqrt2)
rownames(envi_df)=envi_df$age

# rarefied data
sp_01=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-hops-TP-Rarefy_TaxonLevel_N100_2000_AquaticGroups_wide.csv"), row.names = 2, check.names = F)
sp_01=as.data.frame(t(sp_01[-c(1,2)]))
sp_02=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-hops-TP-Rarefy_TaxonLevel_N100_2000_AquaticGroups_wide_prop.csv"), row.names = 2, check.names = F)
sp_02=as.data.frame(t(sp_02[-c(1,2)]))

# select common taxa
df=sp_02
max.abb <- apply(df, 2, max)
n.occ <- colSums(df > 0)
spp.want <- which(max.abb >= 1 & n.occ >= 5)
sp_02.want=sp_02[, spp.want]
sp=sp_01[, colnames(sp_01) %in% colnames(sp_02.want)]

# age
rownames(sp)=envi_df$age

#== output setting
occurrence="rp1_f5_assigned_rare2000"
mi="age_mass"

#== mod1 with rarefied counts
df_Y=round(sp, digits = 0)
df_m=envi_df
sp_mvabund=mvabund(df_Y)

options(scipen=999)
pdf(file = paste0("step03_assigned/04_network/outFigures_aquatic/mod1-meanvar_", ngs, "-", econame, "_", occurrence, ".pdf"), width = 5, height = 5)
meanvar.plot(sp_mvabund, table = T, las = 1, cex.axis = .5)
dev.off()
#You can clearly see that the species with high means (on the x axis) also have high variances (y axis).
# Thus, select family="negative.binomial"
mod1 <- manyglm(sp_mvabund ~ envi_df$age + envi_df$glacier_mass, family="negative.binomial")
# plots
pdf(file = paste0("step03_assigned/04_network/outFigures_aquatic/", "mod1-residual_", ngs, "-", econame, "_", occurrence, ".pdf"), width = 6, height = 5)
plot(mod1)
dev.off()
# anova
mod1.anova=anova(mod1, p.uni="adjusted")
mod1.table=mod1.anova[["table"]]
mod1.uni.p=mod1.anova[["uni.p"]]
mod1.uni.test=mod1.anova[["uni.test"]]
# anova output
write.csv(mod1.table, paste0("step03_assigned/04_network/outTables_aquatic/mod1-table_", ngs, "-", econame, "_", mi, "-occ-", occurrence, ".csv"))
write.csv(mod1.uni.p, paste0("step03_assigned/04_network/outTables_aquatic/mod1-uni-p_", ngs, "-", econame, "_", mi, "-occ-", occurrence, ".csv"))
write.csv(mod1.uni.test, paste0("step03_assigned/04_network/outTables_aquatic/mod1-uni-test_", ngs, "-", econame, "_", mi, "-occ-", occurrence, ".csv"))

#== igraph with rarefied count
output_cgr=cgr(mod1, method = "AIC", seed = 3, 
               #lambda=0, # estimate for time
               n.samp = 1500, n.lambda = 100)
#== out igraph setting
total_count="samp1500"

#== save
inf=output_cgr$best_graph$graph # optimal final network
raw_cov=output_cgr$raw$cov
inf_cov=raw_cov*inf
# partial cor (for graph generation)
part_cor_full=output_cgr$best_graph$part
part_cor_positive=part_cor_full
part_cor_positive[part_cor_positive<0]=0
diag(part_cor_positive)=-1
write.csv(inf, paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph.csv"))
write.csv(raw_cov, paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-raw-cov.csv"))
write.csv(inf_cov, paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-cov.csv"))
write.csv(part_cor_full, paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-full-partial-cov.csv"))
write.csv(part_cor_positive, paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-positive-partial-cov.csv"))
# graph
g=graph_from_adjacency_matrix(as.matrix(part_cor_positive), mode = "undirected", weighted = TRUE, diag = FALSE)
isolated.nodes = which(degree(g)==0)
g2=delete.vertices(g, isolated.nodes)
plot(g)
plot(g2)
#weights
weights=get.data.frame(g2)
weights$model=mi
write.csv(weights, paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-weights-model-", mi, "-occ-", occurrence, "-sumC-", total_count, "-positive.csv"))

#== find optimal 3 modules with max. repetitions (-> optimal 3 modules should be robust)
res=seq(0, 1, 0.01)
for (resi in res) {
  print(resi)
  all_mem_cluster=NULL
  for (i in 1:999) {
    #resi=0.15
    mem=cluster_louvain(g2, resolution = resi, weights = abs(E(g2)$weight))
    num_module=max(mem$membership)
    if (num_module == 3) {
      mem_cluster=data.frame(mem$membership, mem$names)
      mem_cluster=mem_cluster[order(mem_cluster$mem.membership), ]
      mem_cluster$repetitions=i
      all_mem_cluster=rbind(all_mem_cluster, mem_cluster)
    } else {
      next
    }
  }
  if (!is.null(all_mem_cluster)) {
    agg=aggregate(.~mem.names, all_mem_cluster[c("mem.membership", "mem.names")], FUN = sum)
    num_rep=dim(all_mem_cluster)[1]/vcount(g2)
    agg$Nmean=round(agg$mem.membership/num_rep, digits = 0)
    agg=agg[order(agg$Nmean, decreasing = T), ]
    write.csv(agg, paste0("step03_assigned/04_network/outTables_aquatic/", ngs,"-mod1-member_cluster3_", resi, "-", mi, "-occ-", occurrence, "-", total_count, ".csv"))
  } else {
    print(paste0("no 3 module under resolution:", resi))
  }
}
#== index calculation
# load data: partial cor.
part_cor_full=read.csv(paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-full-partial-cov.csv"),
                       row.names = 1)
part_cor_positive=read.csv(paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-positive-partial-cov.csv"),
                           row.names = 1)
# graph with partial full links
g_full=graph_from_adjacency_matrix(as.matrix(part_cor_full), mode = "undirected", weighted = TRUE, diag = FALSE)
isolated.nodes = which(degree(g_full)==0)
g_full=delete.vertices(g_full, isolated.nodes)
g_full_df=as.data.frame(get.edgelist(g_full))
names(g_full_df)=c("from", "to")

# graph with partial positive links
g=graph_from_adjacency_matrix(as.matrix(part_cor_positive), mode = "undirected", weighted = TRUE, diag = FALSE)
isolated.nodes = which(degree(g)==0)
g2=delete.vertices(g, isolated.nodes)

# load module
resi=0.62
cluster=read.csv(paste0("step03_assigned/04_network/outTables_aquatic/", ngs,"-mod1-member_cluster3_", resi, "-", mi, "-occ-", occurrence, "-", total_count, ".csv"), row.names = 1)
cluster=cluster[c(1, 3)]
names(cluster)=c("nodes", "memberships")
cluster$memberships=ifelse(cluster$memberships <= 1, 1, 
                           ifelse(cluster$memberships >1 & cluster$memberships <= 2, 2, 3))
metacluster=NULL
taxacluster=NULL
for (i in 1:max(cluster$memberships)) {
  cluster_sub=cluster[cluster$memberships == i, ]
  cluster_sub_nodes=cluster_sub$nodes
  #
  # full links
  full_cluster_from=g_full_df[grepl(paste0(cluster_sub_nodes, collapse = "|"), g_full_df$from), ]
  full_cluster_to=g_full_df[grepl(paste0(cluster_sub_nodes, collapse = "|"), g_full_df$to), ]
  full_cluster_df=rbind(full_cluster_from, full_cluster_to)
  
  # positive links within module
  part_cor_positive_c=part_cor_positive[colnames(part_cor_positive)%in% cluster_sub$nodes, ]
  part_cor_positive_c=part_cor_positive_c[, colnames(part_cor_positive_c) %in% cluster_sub$nodes]
  #
  # generate graph
  g_cluster=graph_from_adjacency_matrix(as.matrix(part_cor_positive_c), mode = "undirected", weighted = TRUE, diag = F)
  isolated.nodes = which(degree(g_cluster)==0)
  g2_cluster=delete.vertices(g_cluster, isolated.nodes)
  
  # index edges nodes
  edges_cluster=gsize(g2_cluster)
  nodes_cluster=gorder(g2_cluster)
  deg_cluster=degree(g2_cluster)
  
  # degree centralization
  center_deg_cluster=centr_degree(g2_cluster, normalized = TRUE)
  center_deg_res_cluster=as.data.frame(center_deg_cluster$res)
  center_deg_centralization_cluster=center_deg_cluster$centralization
  center_deg_max_cluster=center_deg_cluster$theoretical_max
  
  #Hubs and authorities (should be same if A == AT)
  hs_cluster <- as.data.frame(hub_score(g2_cluster, weights = E(g2_cluster)$weigth)$vector)
  names(hs_cluster)="hs_cluster"
  as_cluster <- as.data.frame(authority_score(g2_cluster, weights = E(g2_cluster)$weigth)$vector)
  names(as_cluster)="as_cluster"
  # out
  metacluster=rbind(metacluster, data.frame(cluster = i, edges=edges_cluster, nodes=nodes_cluster, 
                                            deg_centralization_cluster=center_deg_centralization_cluster, 
                                            positive_links_percent = edges_cluster/nrow(full_cluster_df), 
                                            full_cluster_links = nrow(full_cluster_df)))
  taxacluster=rbind(taxacluster, cbind(cluster_sub$memberships, hs_cluster, as_cluster, stren_cluster))
}
#== save
taxacluster=taxacluster[order(-taxacluster$`cluster_sub$memberships`, -taxacluster$hs_cluster), ]
write.csv(metacluster, paste0("step03_assigned/04_network/outTables_aquatic/", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-metacluster3-res", resi, "-deg.csv"))
write.csv(taxacluster, paste0("step03_assigned/04_network/outTables_aquatic/", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-taxacluster3-res", resi, "-deg.csv"))
# keystone species
keystone_sp=taxacluster[taxacluster$hs_cluster >= 0.8, ]
write.csv(keystone_sp, paste0("step03_assigned/04_network/outTables_aquatic/", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-keystone3-res", resi, ".csv"))

#== cluster and igraph
# for each eco group
raw_igraph=read.csv(paste0("step03_assigned/outTables/", ngs, "-hops-fiter-TP-FW-taxa-to-ecogroup.csv"), row.names = 1)
raw_igraph=raw_igraph[!(raw_igraph$ecogroup %in% c("Terrestrial_Plants", "Terrestrial_Mammals")), ]
raw_igraph=raw_igraph[c("taxa", "ecogroup")]
names(raw_igraph)=c("nodes", "group")
#
module=cluster
igraph=raw_igraph
module_class=left_join(module, igraph, by = "nodes")
#
module_group=NULL
for (i in unique(module_class$memberships)) {
  for (gi in unique(module_class$group)) {
    nodesN=dim(module_class[module_class$memberships == i & module_class$group == gi, ])[1]
    module_group=rbind(module_group, data.frame(memberships=i, group=gi, nodesN=nodesN))
  }
}
#
weights=read.csv(paste0("step03_assigned/04_network/outTables_aquatic/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-FW-weights-model-", mi, "-occ-", occurrence, "-sumC-", total_count, "-positive.csv"),
                 row.names = 1)
names(weights)=c("from_nodes", "to_nodes", "weight", "model")
names(module)=c("from_nodes", "from_community")
weights_module=left_join(weights, module, by = "from_nodes")
names(module)=c("to_nodes", "to_community")
weights_module=left_join(weights_module, module, by = "to_nodes")
#
weights_module$modules=ifelse(weights_module$from_community == weights_module$to_community, "within", "between")
write.csv(weights_module, paste0("step03_assigned/04_network/outTables_aquatic/", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi,"-occ-", occurrence, "-sumC-", total_count, "-weights_module3-res", resi,".csv"))
#
weights_module$from_to_module=paste0(weights_module$from_community, "-", weights_module$to_community)
#
within=weights_module[weights_module$modules == "within", ]
write.csv(within, paste0("step03_assigned/04_network/outTables_aquatic/", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi,"-occ-", occurrence, "-sumC-", total_count, "-within3-res", resi,".csv"))
#
within_df=as.data.frame(table(within$from_to_module))
#
igraph=raw_igraph
names(igraph)=paste0("from_", colnames(igraph))
modules_class=left_join(within, igraph, by = "from_nodes")
#
igraph=raw_igraph
names(igraph)=paste0("to_", colnames(igraph))
modules_class=left_join(modules_class, igraph, by = "to_nodes")
modules_class$from_to_group=paste0(modules_class$from_group, "-", modules_class$to_group)
#
edges_percent_df=NULL
for (i in unique(module_class$memberships)) {
  for (gi in unique(module_class$group)) {
    subdf=modules_class[modules_class$from_to_module == paste0(i, "-", i) & modules_class$from_to_group == paste0(gi, "-", gi), ]
    edges_percent=dim(subdf)[1]/metacluster[metacluster$cluster == i, "full_cluster_links"]
    edges_percent_df=rbind(edges_percent_df, data.frame(memberships=i, group=paste0(gi, "-", gi), edgesN = dim(subdf)[1], edges_percent=edges_percent, full_cluster_links = metacluster[metacluster$cluster == i, "full_cluster_links"]))
  }
}
finaltable=cbind(module_group, edges_percent_df)
write.csv(finaltable, paste0("step03_assigned/04_network/outTables_aquatic/", ngs, "-", econame, "Network-hops-fiter-TP-FW-model", mi,"-occ-", occurrence, "-sumC-", total_count, "-edges_percent3-res", resi,".csv"))
#== network analysis end ==


#== plotting modules in Fig. 3G
# ecogroup
ecogroup=read.csv(paste0("step03_assigned/outTables/", ngs, "-hops-fiter-TP-FW-taxa-to-ecogroup.csv"), row.names = 1)
# data
df1=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-AquaticGroups-N100_2000-p1f5.csv"), row.names = 1)
library(circlize)
library(igraph)
# attach group
df1=as.data.frame(t(df1))
df1$taxa=rownames(df1)
df1=left_join(df1, ecogroup, by = "taxa")
df1$sum=rowSums(df1[1:40])
#
# pre14 ka and post-14 ka, only for select colour
df1$pre14=rowSums(df1[, 1:10])/10
df1$post14=rowSums(df1[, 11:40])/30
df1$pre14_post14=df1$pre14-df1$post14
# check
pre_tx=df1[df1$pre14_post14>0, "taxa"]
pre_tx=pre_tx[-5]
post_tx=df1[df1$pre14_post14<0 & df1$ecogroup == "Aquatic_Bacteria", "taxa"]
post_tx
#
within=read.csv(paste0("step03_assigned/04_network/outTables_aquatic/", ngs, "-aquaticNetwork-hops-fiter-TP-FW-modelage_mass-occ-rp1_f5_assigned_rare2000-sumC-samp1500-within3-res0.62.csv"), row.names = 1)
#== moduels 1, 2, and 3
mi="1-1"
chord=within[within$from_to_module==mi, c("from_nodes", "to_nodes", "weight")]
# plot
if (TRUE) {
  # nodes in netwrok
  nodes=data.frame(taxa=df1$taxa)
  nodes=left_join(nodes, ecogroup, "taxa")
  nodes=nodes[grepl("Aquatic", nodes$ecogroup), ]
  nodes$lifeform=ifelse(nodes$taxa %in% pre_tx, "01_Cyanobacteria",
                        ifelse(nodes$taxa %in% post_tx, "03_Cyanobacteria", nodes$lifeform))
  
  nodes=nodes[order(nodes$ecogroup, nodes$lifeform, nodes$taxa), ]
  nodes$ID=as.character(seq(1:nrow(nodes)))
  nodes=nodes%>%select("ID", everything())
  nodes$ID=paste0(nodes$ID, ".", nodes$taxa)
  #
  igraph=nodes
  names(igraph)[2]=names(chord)[1]
  chord=left_join(chord, igraph, by = "from_nodes")
  names(igraph)[2]=names(chord)[2]
  chord=left_join(chord, igraph, by = "to_nodes")
  #
  names(igraph)[2]= "Species"
  names(igraph)[3]= "group"
  # get nodes
  within_cluster_nodes=unique(c(chord$from_nodes, chord$to_nodes))
  chord_igraph=igraph[igraph$Species %in% within_cluster_nodes, c("ID", "Species", "group")]
  unique(chord_igraph$group)
  #
  sisi_color=c(rgb(41, 64, 146, maxColorValue = 255), # glacier bacteria
               # glacier algae,
               rgb(55, 125, 115, maxColorValue = 255),
               # plant
               rgb(78, 125, 37, maxColorValue = 255),
               # fish
               rgb(223, 132, 50, maxColorValue = 255),
               # mammals
               rgb(235, 82, 45, maxColorValue = 255),
               # freshwater bacteria,
               rgb(91, 45, 95, maxColorValue = 255),
               # others
               rgb(40, 40, 40, maxColorValue = 255))
  # add color
  if (mi == "1-1") {
    chord_igraph$color=ifelse(chord_igraph$group == "02_Cyanobacteria", sisi_color[6], 
                              ifelse(chord_igraph$group == "01_Algae", sisi_color[2], 
                                     ifelse(chord_igraph$group == "05_Fishes", sisi_color[4], 
                                            ifelse(chord_igraph$group == "03_Plants", sisi_color[3],
                                                   ifelse(chord_igraph$group == "01_Cyanobacteria", sisi_color[1],
                                                          ifelse(chord_igraph$group == "03_Cyanobacteria", sisi_color[1], sisi_color[5]))))))
  } else {
    chord_igraph$color=ifelse(chord_igraph$group == "02_Cyanobacteria", sisi_color[6], 
                              ifelse(chord_igraph$group == "01_Algae", sisi_color[2], 
                                     ifelse(chord_igraph$group == "05_Fishes", sisi_color[4], 
                                            ifelse(chord_igraph$group == "03_Plants", sisi_color[3],
                                                   ifelse(chord_igraph$group == "01_Cyanobacteria", sisi_color[6],
                                                          ifelse(chord_igraph$group == "03_Cyanobacteria", sisi_color[6], sisi_color[5]))))))
  }
  # color id
  glacier_color=chord_igraph$color
  names(glacier_color)=chord_igraph$ID
  text.order=chord_igraph$ID
  chord_plot=chord[c("ID.x", "ID.y", "weight")]
  #par(cex = 0.5, mar = c(0, 0, 0, 0))
  circos.clear()
  circos.par(start.degree = 90)
  if (mi != "3-3") {
    chordDiagram(chord_plot, directional=F, transparency = 0.5, 
                 grid.col=glacier_color, 
                 #big.gap = 140, # for post-3.6
                 order = text.order)
  } else {
    chordDiagram(chord_plot, directional=F, transparency = 0.5, 
                 grid.col=glacier_color, 
                 big.gap = 140, # for post-3.6
                 order = text.order)
  }
}
#== plotting end ==














