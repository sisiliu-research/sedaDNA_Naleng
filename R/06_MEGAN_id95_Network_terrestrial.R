# by Sisi Liu: Sisi Liu: sisi.liu@awi.de | sisi.liu.research@gmail.com (permanent address)
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1

# model-based network

#== setting
dir.create("step03_assigned/04_network/")
dir.create("step03_assigned/04_network/outTables_land")
dir.create("step03_assigned/04_network/outFigures_land")
ngs="APMG-5-10-28-34-35"
econame="Terrestrial"

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
envi=read.csv("step03_assigned/03_rda/terrestrial_envi.csv")
names(envi)
envi_df=envi[c("age", "Permafrost_catchment_approx", "glaciers_area_approx")]
rownames(envi_df)=envi_df$age
# plant
sp_01=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-hops-TP-Rarefy_TaxonLevel_N100_3000_TerrestrialPlants_wide.csv"), row.names = 2, check.names = F)
sp_01=as.data.frame(t(sp_01[-c(1,2)]))
sp_02=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-hops-TP-Rarefy_TaxonLevel_N100_3000_TerrestrialPlants_wide_prop.csv"), row.names = 2, check.names = F)
sp_02=as.data.frame(t(sp_02[-c(1,2)]))
# mammals
sp_03=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-hops-TP-Rarefy_TaxonLevel_N100_150_TerrestrialMamm_wide.csv"), row.names = 2, check.names = F)
sp_03=as.data.frame(t(sp_03[-c(1,2)]))
sp_04=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-hops-TP-Rarefy_TaxonLevel_N100_150_TerrestrialMamm_wide_prop.csv"), row.names = 2, check.names = F)
sp_04=as.data.frame(t(sp_04[-c(1,2)]))

#== select common plant taxa
df=sp_02
max.abb <- apply(df, 2, max)
n.occ <- colSums(df > 0)
spp.want <- which(max.abb >= 1 & n.occ >= 5)
sp_02.want=sp_02[, spp.want]
sp_01.want=sp_01[, colnames(sp_01) %in% colnames(sp_02.want)]
# select common mammalian taxa
df=sp_04
max.abb <- apply(df, 2, max)
n.occ <- colSums(df > 0)
spp.want <- which(max.abb >= 1 & n.occ >= 5)
sp_04.want=sp_04[, spp.want]
sp_03.want=sp_03[, colnames(sp_03) %in% colnames(sp_04.want)]
sp=cbind(sp_01.want, sp_03.want)
# uniform age
rownames(sp)=envi_df$age

#== setting of option
occurrence="rp1_f5_assigned_p3000_m150"
mi="age_cry"

#== mod1 with rarefied counts
df_Y=round(sp, digits = 0)
df_m=envi_df
sp_mvabund=mvabund(df_Y)

options(scipen=999)
pdf(file = paste0("step03_assigned/04_network/outFigures_land/mod1-meanvar_", ngs, "-", econame, "_", occurrence, ".pdf"), width = 5, height = 5)
meanvar.plot(sp_mvabund, table = T, las = 1, cex.axis = .5)
dev.off()
#You can clearly see that the species with high means (on the x axis) also have high variances (y axis).
# Thus, select family="negative.binomial"
mod1 <- manyglm(sp_mvabund ~ envi_df$age + envi_df$Permafrost_catchment_approx + envi_df$glaciers_area_approx,
                family="negative.binomial")
# plots
pdf(file = paste0("step03_assigned/04_network/outFigures_land/", "mod1-residual_", ngs, "-", econame, "_", occurrence, ".pdf"), width = 6, height = 5)
plot(mod1)
dev.off()

# anova
mod1.anova=anova(mod1, p.uni="adjusted")
mod1.table=mod1.anova[["table"]]
mod1.uni.p=mod1.anova[["uni.p"]]
mod1.uni.test=mod1.anova[["uni.test"]]
# anova output
write.csv(mod1.table, paste0("step03_assigned/04_network/outTables_land/mod1-table_", ngs, "-", econame, "_", mi, "-occ-", occurrence,".csv"))
write.csv(mod1.uni.p, paste0("step03_assigned/04_network/outTables_land/mod1-uni-p_", ngs, "-", econame, "_", mi, "-occ-", occurrence, ".csv"))
write.csv(mod1.uni.test, paste0("step03_assigned/04_network/outTables_land/mod1-uni-test_", ngs, "-", econame, "_", mi, "-occ-", occurrence, ".csv"))

#== igraph with rarefied count
# graph
output_cgr=cgr(mod1, method = "AIC", seed = 3, 
               #lambda=0, # estimate for time
               n.samp = 1500, n.lambda = 100)
#== out setting of igraph
total_count="samp1500"
#== save
inf=output_cgr$best_graph$graph # optimal final network
raw_cov=output_cgr$raw$cov
inf_cov=raw_cov*inf
part_cor_full=output_cgr$best_graph$part
part_cor_positive=part_cor_full
part_cor_positive[part_cor_positive<0]=0
diag(part_cor_positive)=-1
# save
write.csv(inf, paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph.csv"))
write.csv(raw_cov, paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-raw-cov.csv"))
write.csv(inf_cov, paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-cov.csv"))
write.csv(part_cor_full, paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-full-partial-cov.csv"))
write.csv(part_cor_positive, paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-positive-partial-cov.csv"))

#== graph with partial positive links
g=graph_from_adjacency_matrix(as.matrix(part_cor_positive), mode = "undirected", weighted = TRUE, diag = FALSE)
isolated.nodes = which(degree(g)==0)
g2=delete.vertices(g, isolated.nodes)
plot(g)
plot(g2)
#== weights
weights=get.data.frame(g2)
weights$model=mi
write.csv(weights, paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-weights-model-", mi, "-occ-", occurrence, "-sumC-", total_count, "-positive.csv"))

#== find optimal 2 modules with max. repetitions and repeats of res (-> optimal 2 modules should be robust)
res=seq(0, 1, 0.01)
for (resi in res) {
  print(resi)
  all_mem_cluster=NULL
  for (i in 1:999) {
    #resi=0.2
    mem=cluster_louvain(g2, resolution = resi, weights = abs(E(g2)$weight))
    num_module=max(mem$membership)
    if (num_module == 2) {
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
    num_rep=dim(all_mem_cluster)[1]/68
    agg$Nmean=round(agg$mem.membership/num_rep, digits = 0)
    agg=agg[order(agg$Nmean, decreasing = T), ]
    write.csv(agg, paste0("step03_assigned/04_network/outTables_land/", ngs, "-mod1-member_cluster_", resi, "-", mi, "-occ-", occurrence, "-", total_count, ".csv"))
  } else {
    print(paste0("no 2 module under resolution:", resi))
  }
}
#== index calculation
# load data: partial cor.
part_cor_full=read.csv(paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-full-partial-cov.csv"), row.names = 1)
part_cor_positive=read.csv(paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-cgr-best-graph-positive-partial-cov.csv"), row.names = 1)

# graph with partial full links
g_full=graph_from_adjacency_matrix(as.matrix(part_cor_full), mode = "undirected", weighted = TRUE, diag = FALSE)
isolated.nodes = which(degree(g_full)==0)
g_full=delete.vertices(g_full, isolated.nodes)
g_full_df=as.data.frame(get.edgelist(g_full))
names(g_full_df)=c("from", "to")

#= graph with partial positive links
g=graph_from_adjacency_matrix(as.matrix(part_cor_positive), mode = "undirected", weighted = TRUE, diag = FALSE)
isolated.nodes = which(degree(g)==0)
g2=delete.vertices(g, isolated.nodes)

# load module:
resi=0.38
cluster=read.csv(paste0("step03_assigned/04_network/outTables_land/", ngs, "-mod1-member_cluster_", resi, "-", mi, "-occ-", occurrence, "-", total_count, ".csv"), row.names = 1)
cluster=cluster[c(1, 3)]
names(cluster)=c("nodes", "memberships")
metacluster=NULL
taxacluster=NULL
for (i in 1:max(cluster$memberships)) {
  cluster_sub=cluster[cluster$memberships == i, ]
  cluster_sub_nodes=cluster_sub$nodes
  # full links
  full_cluster_from=g_full_df[grepl(paste0(cluster_sub_nodes, collapse = "|"), g_full_df$from), ]
  full_cluster_to=g_full_df[grepl(paste0(cluster_sub_nodes, collapse = "|"), g_full_df$to), ]
  full_cluster_df=rbind(full_cluster_from, full_cluster_to)
  
  # positive links within module
  part_cor_positive_c=part_cor_positive[colnames(part_cor_positive)%in% cluster_sub$nodes, ]
  part_cor_positive_c=part_cor_positive_c[, colnames(part_cor_positive_c) %in% cluster_sub$nodes]
  #
  # index edges nodes
  g_cluster=graph_from_adjacency_matrix(as.matrix(part_cor_positive_c), mode = "undirected", weighted = TRUE, diag = FALSE)
  isolated.nodes = which(degree(g_cluster)==0)
  g2_cluster=delete.vertices(g_cluster, isolated.nodes)
  # plotting
  pdf(file = paste0("step03_assigned/04_network/outFigures_land/", ngs, "-mod1-member_cluster_", i, "_res",resi, "-", mi, "-occ-", occurrence, "-", total_count, ".pdf"), width = 10, height = 8)
  plot(g2_cluster)
  dev.off()
  
  # graph characters
  edges_cluster=gsize(g2_cluster)
  nodes_cluster=gorder(g2_cluster)
  deg_cluster=degree(g2_cluster)
  
  # degree,number of ties
  center_deg_cluster=centr_degree(g2_cluster, normalized = TRUE)
  center_deg_res_cluster=as.data.frame(center_deg_cluster$res)
  center_deg_centralization_cluster=center_deg_cluster$centralization
  center_deg_max_cluster=center_deg_cluster$theoretical_max
  
  # Hubs and authorities (should be same if A == AT)
  hs_cluster <- as.data.frame(hub_score(g2_cluster, weights = E(g2_cluster)$weigth)$vector)
  names(hs_cluster)="hs_cluster"
  as_cluster <- as.data.frame(authority_score(g2_cluster, weights = E(g2_cluster)$weigth)$vector)
  names(as_cluster)="as_cluster"
  
  #
  metacluster=rbind(metacluster, data.frame(cluster = i, edges=edges_cluster, nodes=nodes_cluster, 
                                            deg_centralization_cluster=center_deg_centralization_cluster, 
                                            positive_links_percent = edges_cluster/dim(full_cluster_df)[1], 
                                            full_cluster_links = dim(full_cluster_df)[1]))
  taxacluster=rbind(taxacluster, cbind(cluster_sub$memberships, hs_cluster, as_cluster, stren_cluster))
}
#== save
taxacluster=taxacluster[order(-taxacluster$`cluster_sub$memberships`, -taxacluster$hs_cluster), ]
write.csv(metacluster, paste0("step03_assigned/04_network/outTables_land/", ngs, "-", econame, "Network-hops-filter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-metacluster-res", resi, "-deg.csv"))
write.csv(taxacluster, paste0("step03_assigned/04_network/outTables_land/", ngs, "-", econame, "Network-hops-filter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-taxacluster-res", resi, "-deg.csv"))
# keystone species
keystone_sp=taxacluster[taxacluster$hs_cluster >= 0.8, ]
write.csv(keystone_sp, paste0("step03_assigned/04_network/outTables_land/", ngs, "-", econame, "Network-hops-filter-TP-model", mi, "-occ-", occurrence, "-sumC-", total_count,"-keystone-res", resi, ".csv"))

#== for each ecogroup
raw_igraph=read.csv(paste0("step03_assigned/outTables/", ngs, "-hops-filter-TP-FW-taxa-to-ecogroup.csv"), row.names = 1)
raw_igraph=raw_igraph[raw_igraph$ecogroup %in% c("Terrestrial_Plants", "Terrestrial_Mammals"), ]
raw_igraph=raw_igraph[c("taxa", "ecogroup")]
names(raw_igraph)=c("nodes", "group")

#== cluster and igraph
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
weights=read.csv(paste0("step03_assigned/04_network/outTables_land/mod1-", ngs, "-", econame, "Network-hops-fiter-TP-weights-model-", mi, "-occ-", occurrence, "-sumC-", total_count, "-positive.csv"), row.names = 1)
names(weights)=c("from_nodes", "to_nodes", "weight", "model")
names(module)=c("from_nodes", "from_community")
weights_module=left_join(weights, module, by = "from_nodes")
names(module)=c("to_nodes", "to_community")
weights_module=left_join(weights_module, module, by = "to_nodes")
#
weights_module$modules=ifelse(weights_module$from_community == weights_module$to_community, "within", "between")
write.csv(weights_module, paste0("step03_assigned/04_network/outTables_land/", ngs, "-", econame, "Network-hops-fiter-TP-model", mi,"-occ-", occurrence, "-sumC-", total_count, "-weights_module-res", resi,".csv"))
#
weights_module$from_to_module=paste0(weights_module$from_community, "-", weights_module$to_community)
#
within=weights_module[weights_module$modules == "within", ]
write.csv(within, paste0("step03_assigned/04_network/outTables_land/", ngs, "-", econame, "Network-hops-fiter-TP-model", mi,"-occ-", occurrence, "-sumC-", total_count, "-within-res", resi,".csv"))

#== index for module
within_df=as.data.frame(table(within$from_to_module))
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
write.csv(finaltable, paste0("step03_assigned/04_network/outTables_land/", ngs, "-", econame, "Network-hops-fiter-TP-model", mi,"-occ-", occurrence, "-sumC-", total_count, "-edges_percent-res", resi,".csv"))
#== network analysis end ==


#== plotting modules in Fig. 2L
library(circlize)
library(igraph)

# ecogroup
ecogroup=read.csv(paste0("step03_assigned/outTables/", ngs, "-hops-filter-TP-FW-taxa-to-ecogroup.csv"), row.names = 1)
# plants
df1=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-TerrestrialPlants-N100_3000-p1f5.csv"), row.names = 1)
df2=read.csv(paste0("step03_assigned/02_rarefy/outTables/", ngs, "-TerrestrialMamm-N100_150-p1f5.csv"), row.names = 1)

# nodes in netwrok
nodes=data.frame(taxa=c(names(df1), names(df2)))
nodes=left_join(nodes, ecogroup, "taxa")
nodes=nodes[grepl("Terrestrial", nodes$ecogroup), ]
nodes=nodes[order(nodes$ecogroup, nodes$lifeform, nodes$taxa), ]
nodes$ID=as.character(seq(1:nrow(nodes)))
nodes$ID=paste0(nodes$ID, ".", nodes$taxa)
#
sisi_color=c(rgb(242, 159, 64, maxColorValue = 255), # forbs
             # gr
             rgb(196, 204, 68, maxColorValue = 255),
             # wood
             rgb(157, 198, 175, maxColorValue = 255),
             # wide growth
             rgb(254, 253, 197, maxColorValue = 255),
             # mammals forb
             rgb(216, 68, 41, maxColorValue = 255),
             # mammals wood
             rgb(233, 187, 197, maxColorValue = 255),
             # wide
             rgb(242, 220, 178, maxColorValue = 255))

#== load data: Line 281 
within=read.csv(past0("step03_assigned/04_network/outTables_land/", ngs, "-TerrestrialNetwork-hops-fiter-TP-modelage_cry-occ-rp1_f5_assigned_p3000_m150-sumC-samp1500-within-res0.38.csv"), row.names = 1)

#== show module 1 and 2
chord=within[within$from_to_module=="1-1", c("from_nodes", "to_nodes", "weight")]
#chord=within[within$from_to_module=="2-2", c("from_nodes", "to_nodes", "weight")]
igraph=nodes
igraph=igraph%>%select("ID", everything())
names(igraph)[2]=names(chord)[1]
chord=left_join(chord, igraph, by = "from_nodes")
names(igraph)[2]=names(chord)[2]
chord=left_join(chord, igraph, by = "to_nodes")
#
names(igraph)[2]= "Species"
names(igraph)[3]="group"
# get nodes
within_cluster_nodes=unique(c(chord$from_nodes, chord$to_nodes))
chord_igraph=igraph[igraph$Species %in% within_cluster_nodes, c("ID", "Species", "group")]
unique(chord_igraph$group)
chord_igraph$color=ifelse(chord_igraph$group == "06_Forbs", sisi_color[1], 
                          ifelse(chord_igraph$group == "07_Graminoid", sisi_color[2], 
                                 ifelse(chord_igraph$group == "09_Wide growth form", sisi_color[4], 
                                        ifelse(chord_igraph$group == "08_Woody", sisi_color[3],
                                               ifelse(chord_igraph$group == "10_Steppe and meadow", sisi_color[5], 
                                                      ifelse(chord_igraph$group == "11_Woodland", sisi_color[6], sisi_color[7]))))))
glacier_color=chord_igraph$color
names(glacier_color)=chord_igraph$ID
text.order=chord_igraph$ID
chord_plot=chord[c("ID.x", "ID.y", "weight")]
#par(cex = 0.5, mar = c(0, 0, 0, 0))
circos.clear()
circos.par(start.degree = 90)
chordDiagram(chord_plot, directional=F, transparency = 0.5, 
             grid.col=glacier_color,
             order = text.order)
#== plotting end ==










