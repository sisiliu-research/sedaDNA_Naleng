# by Sisi Liu: sisi.liu@awi.de; sisi.liu.research@gmail.com
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS 14.1

#==== Refs and Tutorial ====
# Popovic, G. C., Hui, F. K. C., and Warton, D. I. (2018). A general algorithm for covariance modeling of discrete data. Journal of Multivariate Analysis 165, 86–100. doi: 10.1016/j.jmva.2017.12.002.
# Popovic, G. C., Warton, D. I., Thomson, F. J., Hui, F. K. C., and Moles, A. T. (2019). Untangling direct species associations from indirect mediator species effects with graphical models. Methods in Ecology and Evolution 10, 1571–1583. doi: 10.1111/2041-210X.13247.
# https://cran.r-project.org/web/packages/ecoCopula/vignettes/the_basics.html

#==== load packages ====
library(readr)
library(dplyr)
library(reshape)
library(DESeq2)
library(mvabund)
library(ggplot2)
library(ecoCopula)
library(igraph)

#==== output ====
dir.create("11_ngsLCA_L30_post/03_offset/")
dir.create("11_ngsLCA_L30_post/03_offset/00Tables")
dir.create("11_ngsLCA_L30_post/03_offset/00Figures")

#==== load ngs, Source: Dataset 1 ====
indf0=read.csv(paste0("11_ngsLCA_L30_post/00_data/Tables/RAW_COUNT_PLANT.csv"), check.names = F, row.names = 2)
indf1=read.csv(paste0("11_ngsLCA_L30_post/00_data/Tables/RAW_COUNT_MAM.csv"), check.names = F, row.names = 2)
indf=rbind(indf0, indf1)

# convert to samples x taxa
rdf=as.data.frame(t(indf[-c(1, 2)]))

#==== load envi data, Source: sheet = ENVI in Dataset 1  ====
eco="terrestrial"
envi=read.csv(paste0("11_ngsLCA_L30_post/04_rda/ENVI_", eco, ".csv"), row.names = 1)[c("Glaciers_area", "Permafrost_catchment", "Land_use")]
rownames(envi)=rownames(rdf)

# find size factor
spe=as.data.frame(t(rdf))
countData <- as.matrix(spe)
condition <- factor(seq_len(ncol(spe)))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), design = ~ condition)
dds <- estimateSizeFactors(dds, type = c("ratio"))
esf=sizeFactors(dds)
envi$sizeFactors=esf

#==== initial H0 vs. H1 models ====
sp_mvabund=mvabund(rdf)
mod0 = stackedsdm(sp_mvabund, ~ 1 + offset(log(sizeFactors)), data=envi, family = "negative.binomial")
mod1 = stackedsdm(sp_mvabund, ~ Glaciers_area + Permafrost_catchment + Land_use + offset(log(sizeFactors)), data=envi, family = "negative.binomial")

# latent factor analysis
ord_mod0=cord(mod0, seed=123)
ord_mod1=cord(mod1, seed=123)

# biplot 
pdf(paste0("11_ngsLCA_L30_post/03_offset/00Figures/cord_bioplot.pdf"), width = 22, height = 16)  # Open PDF file
par(mfrow=c(1,2))
# samples
ice<-ifelse(envi$age > 14, "blue", ifelse(envi$age < 14 & envi$age > 3.6, "red", "skyblue")) 
# H0
plot(ord_mod0, site.col = ice, biplot=T, cex.lab = .5)
title(main = "Null hypothesis (H0)\n1 + offset(log(sizeFactors))", cex.main = 0.5)
# H1
plot(ord_mod1, site.col = ice, biplot=T, cex.lab = .5)
title(main = "Hypothesis (H1)\nGlaciers extent + Permafrost area + Land use + offset(log(sizeFactors))", cex.main = 0.5)
dev.off()

#==== determine abundance limit (until generating the marginal network) and res (identical for H0 and H1) ====
abun_effect=data.frame(taxa=colnames(rdf), total_abun=colSums(rdf))
abun_seq=seq(0, 5000, 50)
# computing inefficiently for j < 750
for (j in abun_seq) {
  print(j)
  abun_tax=abun_effect[abun_effect$total_abun > j, ]
  # subset
  sub_df=rdf[, colnames(rdf) %in% abun_tax$taxa]
  sub_sp_mvabund=mvabund(sub_df)
  # mods
  mod2=stackedsdm(sub_sp_mvabund, ~ 1 + offset(log(sizeFactors)), data=envi, family = "negative.binomial")
  mod3=stackedsdm(sub_sp_mvabund, ~ Glaciers_area + Permafrost_catchment + Land_use + offset(log(sizeFactors)), data=envi, family = "negative.binomial")
  
  #== can fit a completely dense marginal network?
  graph_mod2 = cgr(mod2, 
                   #lambda=0, # use lambda=0 for estimate time, normally < 5 min
                   method = "BIC", seed = 123, n.samp = 500, n.lambda = 100)
  
  graph_mod3 = cgr(mod3, graph_mod4$all_graphs$lambda.opt, 
                   method = "BIC", seed = 123, n.samp = 500)
  
  #== can find 2 modules?
  # H0:
  output_cgr=graph_mod2
  occ=paste0("graph_mod2_a", j)
  #
  inf=output_cgr$best_graph$graph # optimal final network
  raw_cov=output_cgr$raw$cov
  inf_cov=raw_cov*inf
  best_cov=output_cgr$best_graph$cov
  part_cor_full=output_cgr$best_graph$part
  part_cor_positive=part_cor_full
  part_cor_positive[part_cor_positive<0]=0
  diag(part_cor_positive)=-1
  #
  partial=part_cor_full
  diag(partial)<-0
  vertex <- data.frame(name = colnames(partial))
  Edges=igraph::as_data_frame(igraph::graph_from_adjacency_matrix(partial,mode="undirected",weighted = TRUE))
  if(ncol(Edges)==3){
    colnames(Edges)[3]<- "partcor"
  }
  write.csv(Edges, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-Edges.csv"))
  # positive links
  Edges_pov=Edges[Edges$partcor > 0, ]
  best_pov_g=igraph::graph_from_data_frame(Edges_pov, directed=FALSE, vertices=vertex)
  plot(best_pov_g)
  
  isolated = which(degree(best_pov_g)==0)
  g2 = delete.vertices(best_pov_g, isolated)
  plot(g2)
  
  # keep resi with 2 modules
  res=seq(0, 1, 0.01)
  for (resi in res) {
    print(resi)
    all_mem_cluster=NULL
    all_modularity=NULL
    for (m in 1:999) {
      #resi=0.2
      set.seed(123)
      mem=cluster_louvain(g2, resolution = resi)
      
      num_module=max(mem$membership)
      if (num_module == 2) {
        
        mem_cluster=data.frame(mem$membership, mem$names)
        mem_cluster=mem_cluster[order(mem_cluster$mem.membership), ]
        mem_cluster$repetitions=m
        all_mem_cluster=rbind(all_mem_cluster, mem_cluster)
        
      } else {
        next
      }
    }
    if (!is.null(all_mem_cluster)) {
      agg=aggregate(.~mem.names, all_mem_cluster[c("mem.membership", "mem.names")], FUN = sum)
      num_rep=dim(all_mem_cluster)[1]/nrow(agg)
      agg$Nmean=round(agg$mem.membership/num_rep, digits = 0)
      agg=agg[order(agg$Nmean, decreasing = T), ]
      write.csv(agg, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-member_cluster-resi-", resi,".csv"))
      
    } else {
      print(paste0("no 2 module under resolution:", resi))
    }
  }
  
  # H1:
  output_cgr=graph_mod3
  occ=paste0("graph_mod3_a", j)
  
  # optimal final network
  inf=output_cgr$best_graph$graph 
  raw_cov=output_cgr$raw$cov
  inf_cov=raw_cov*inf
  best_cov=output_cgr$best_graph$cov
  part_cor_full=output_cgr$best_graph$part
  part_cor_positive=part_cor_full
  part_cor_positive[part_cor_positive<0]=0
  diag(part_cor_positive)=-1
  
  # save positive and negative edges
  partial=part_cor_full
  diag(partial)<-0
  vertex <- data.frame(name = colnames(partial))
  Edges=igraph::as_data_frame(igraph::graph_from_adjacency_matrix(partial,mode="undirected",weighted = TRUE))
  if(ncol(Edges)==3){
    colnames(Edges)[3]<- "partcor"
  }
  write.csv(Edges, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-Edges.csv"))
  
  # positive links
  Edges_pov=Edges[Edges$partcor > 0, ]
  best_pov_g=igraph::graph_from_data_frame(Edges_pov, directed=FALSE, vertices=vertex)
  isolated = which(degree(best_pov_g)==0)
  g2 = delete.vertices(best_pov_g, isolated)
  # keep resi with 2 modules
  res=seq(0, 1, 0.01)
  for (resi in res) {
    print(resi)
    all_mem_cluster=NULL
    all_modularity=NULL
    for (m in 1:999) {
      #resi=0.2
      set.seed(123)
      mem=cluster_louvain(g2, resolution = resi)
      num_module=max(mem$membership)
      if (num_module == 2) {
        mem_cluster=data.frame(mem$membership, mem$names)
        mem_cluster=mem_cluster[order(mem_cluster$mem.membership), ]
        mem_cluster$repetitions=m
        all_mem_cluster=rbind(all_mem_cluster, mem_cluster)
      } else {
        next
      }
    }
    if (!is.null(all_mem_cluster)) {
      agg=aggregate(.~mem.names, all_mem_cluster[c("mem.membership", "mem.names")], FUN = sum)
      num_rep=dim(all_mem_cluster)[1]/nrow(agg)
      agg$Nmean=round(agg$mem.membership/num_rep, digits = 0)
      agg=agg[order(agg$Nmean, decreasing = T), ]
      write.csv(agg, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-member_cluster-resi-", resi,".csv"))
    } else {
      print(paste0("no 2 module under resolution:", resi))
    }
  }
  
}
# output: a completely dense marginal network can be generated when j = 2450, identical resi of H0 and H1 is found also.

#==== Network with common taxa ====
i=2450
abun_tax=abun_effect[abun_effect$total_abun > i, ]
# subset
sub_df=rdf[, colnames(rdf) %in% abun_tax$taxa]
sub_sp_mvabund=mvabund(sub_df)
# mods
mod4=stackedsdm(sub_sp_mvabund, ~ 1 + offset(log(sizeFactors)), data=envi, family = "negative.binomial")
mod5=stackedsdm(sub_sp_mvabund, ~ Glaciers_area + Permafrost_catchment + Land_use + offset(log(sizeFactors)), data=envi, family = "negative.binomial")

# H0: Dunn-Smyth Residuals
mod4_res=as.data.frame(t(as.data.frame(residuals(mod4))))
mod4_predict=as.data.frame(t(as.data.frame(predict(mod4))))
mod4_res$taxa=rownames(mod4_res)
mod4_predict$taxa=rownames(mod4_predict)
# long
mod4_res=melt(mod4_res, "taxa")
mod4_predict=melt(mod4_predict, "taxa")
mod4_fits=cbind(mod4_res, mod4_predict[, 3])
names(mod4_fits)=c("taxa", "variable", "value", "predict")
# plot
p1=ggplot(mod4_fits, aes(x=predict, y=value, group=taxa)) +
  scale_x_continuous("Linear predictor values") +
  scale_y_continuous("Dunn-Smyth Residuals") +
  geom_point(aes(color=taxa), show.legend = F) +
  theme_bw()+
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x.bottom = element_text(size = 8, color = "black", angle = 0))
# H1: Dunn-Smyth Residuals
mod5_res=as.data.frame(t(as.data.frame(residuals(mod5))))
mod5_predict=as.data.frame(t(as.data.frame(predict(mod5))))
mod5_res$taxa=rownames(mod5_res)
mod5_predict$taxa=rownames(mod5_predict)
# long
mod5_res=melt(mod5_res, "taxa")
mod5_predict=melt(mod5_predict, "taxa")
mod5_fits=cbind(mod5_res, mod5_predict[, 3])
names(mod5_fits)=c("taxa", "variable", "value", "predict")
# plot
p2=ggplot(mod5_fits, aes(x=predict, y=value, group=taxa)) +
  scale_x_continuous("Linear predictor values") +
  scale_y_continuous("Dunn-Smyth Residuals") +
  geom_point(aes(color=taxa), show.legend = F) +
  theme_bw()+
  theme(axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x.bottom = element_text(size = 8, color = "black", angle = 0))
# Histogram
hist(residuals(mod4))
hist(residuals(mod5))
# Q-Q
qqnorm(residuals(mod4))
qqnorm(residuals(mod5))

#==== Fitting Gaussian copula graphical lasso to co-occurrence data ====
graph_mod4 = cgr(mod4, 
                 #lambda=0, # use lambda=0 for estimate time, normally < 5 min
                 method = "BIC", seed = 123, n.samp = 500, n.lambda = 100)

graph_mod5 = cgr(mod5, graph_mod4$all_graphs$lambda.opt, 
                 method = "BIC", seed = 123, n.samp = 500)

# save the best graph
output_cgr=graph_mod5
occ=paste0("graph_mod5_a", i)
#
inf=output_cgr$best_graph$graph 
raw_cov=output_cgr$raw$cov
inf_cov=raw_cov*inf
best_cov=output_cgr$best_graph$cov
part_cor_full=output_cgr$best_graph$part
part_cor_positive=part_cor_full
part_cor_positive[part_cor_positive<0]=0
diag(part_cor_positive)=-1

# save positive and negative edges
partial=part_cor_full
diag(partial)<-0
vertex <- data.frame(name = colnames(partial))
Edges=igraph::as_data_frame(igraph::graph_from_adjacency_matrix(partial,mode="undirected",weighted = TRUE))
if(ncol(Edges)==3){
  colnames(Edges)[3]<- "partcor"
}
write.csv(Edges, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-Edges.csv"))

# positive links
Edges_pov=Edges[Edges$partcor > 0, ]
best_pov_g=igraph::graph_from_data_frame(Edges_pov, directed=FALSE, vertices=vertex)
isolated = which(degree(best_pov_g)==0)
g2 = delete.vertices(best_pov_g, isolated)
# 
res=0.69 # the lowest resi of H1 is identical to H0's.
for (resi in res) {
  print(resi)
  all_mem_cluster=NULL
  all_modularity=NULL
  for (m in 1:999) {
    #resi=0.2
    set.seed(123)
    mem=cluster_louvain(g2, resolution = resi)
    num_module=max(mem$membership)
    if (num_module == 2) {
      mem_cluster=data.frame(mem$membership, mem$names)
      mem_cluster=mem_cluster[order(mem_cluster$mem.membership), ]
      mem_cluster$repetitions=m
      all_mem_cluster=rbind(all_mem_cluster, mem_cluster)
    } else {
      next
    }
  }
  if (!is.null(all_mem_cluster)) {
    agg=aggregate(.~mem.names, all_mem_cluster[c("mem.membership", "mem.names")], FUN = sum)
    num_rep=dim(all_mem_cluster)[1]/nrow(agg)
    agg$Nmean=round(agg$mem.membership/num_rep, digits = 0)
    agg=agg[order(agg$Nmean, decreasing = T), ]
    write.csv(agg, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-member_cluster-resi-", resi,".csv"))
  } else {
    print(paste0("no 2 module under resolution:", resi))
  }
}

#==== degree and keystone ====
  cluster=read.csv(paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-member_cluster-resi-", resi, ".csv"), row.names = 1)
  cluster=cluster[c(1, 3)]
  names(cluster)=c("nodes", "memberships")
  metacluster=NULL
  taxacluster=NULL
  for (k in 1:max(cluster$memberships)) {
    print(k)
    cluster_sub=cluster[cluster$memberships == k, ]
    cluster_sub_nodes=cluster_sub$nodes
    
    # full edges
    full_cluster_from=Edges[Edges$from %in% cluster_sub_nodes, ]
    full_cluster_from_to=full_cluster_from[full_cluster_from$to %in% cluster_sub_nodes, ]
    
    # pov within
    Edges_pov_from=Edges_pov[Edges_pov$from %in% cluster_sub_nodes, ]
    Edges_pov_from_to=Edges_pov_from[Edges_pov_from$to %in% cluster_sub_nodes, ]
    
    # graph with pov edges
    g=graph_from_data_frame(Edges_pov_from_to, directed=FALSE, vertices=vertex)
    isolated.nodes = which(degree(g)==0)
    mg=delete.vertices(g, isolated.nodes)
    
    # number of edges
    edges_cluster=gsize(mg)
    # number of nodes (taxa in this study)
    nodes_cluster=gorder(mg)
    # node's edges
    deg_cluster=degree(mg)
    
    # Centralize a graph
    center_deg_cluster=centr_degree(mg, normalized = TRUE)
    center_deg_res_cluster=as.data.frame(center_deg_cluster$res)
    center_deg_centralization_cluster=center_deg_cluster$centralization
    
    # Hubs
    hs_cluster <- as.data.frame(hub_score(mg, weights = E(mg)$weigth)$vector)
    names(hs_cluster)="hs_cluster"
  
    # normalize positive links by number of nodes within module
    metacluster=rbind(metacluster, data.frame(cluster = k, edges=edges_cluster, nodes=nodes_cluster, 
                                              deg_centralization_cluster=center_deg_centralization_cluster, 
                                              positive_links_percent = edges_cluster/dim(full_cluster_from_to)[1], 
                                              full_cluster_links = dim(full_cluster_from_to)[1]))
    taxacluster=rbind(taxacluster, cbind(cluster_sub$memberships, hs_cluster))
  }
  taxacluster=taxacluster[order(-taxacluster$`cluster_sub$memberships`, -taxacluster$hs_cluster), ]
  write.csv(metacluster, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-resi-", resi, "-metacluster.csv"))
  write.csv(taxacluster, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-resi-", resi, "-taxacluster.csv"))
  
  # keystone species
  keystone_sp=taxacluster[taxacluster$hs_cluster >= 0.8, ]
  write.csv(keystone_sp, paste0("11_ngsLCA_L30_post/03_offset/00Tables/", occ, "-resi-", resi, "-keystone.csv"))
#== END ==  
  
  
  
  
  
  
  


