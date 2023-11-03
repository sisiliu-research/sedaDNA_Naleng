# by Sisi Liu: sisi.liu@awi.de | sisi.liu.research@gmail.com (permanent address)
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1

# Permafrost simulation: GLM: 22 - 0 ka
library(readr)
library(raster)
library(sp)
library(rgdal)
library(sf)
library(raster)
library(dismo)
library(rasterVis)
library(ggplot2)
#
setwd("/permafrost_glm")
dir.create("predict")
#== load temperature/source data: https://github.com/StefanKruse/R_PastElevationChange
paleo_temp=read.csv2("tempRCP45_only500yrsteps_linearinterpolated_22ka.csv")
#> head(paleo_temp)
#yr.
#1    0,0.225
#2   500,0.38
#3   1000,0.8
#4 1500,0.145
#5 2000,0.165
#6 2500,0.235

#=============== Build present-day GLM model =============
# load dem data
#source data are downloaded from: 
# 1. bio - https://www.worldclim.org/data/worldclim21.html; 
# 2. Permafrost - https://tc.copernicus.org/articles/11/2527/2017/; 
# 3. Hengduan dem - SRTM-Downloader plugin in QGIS software https://qgis.org/de/site/
bio1_sp_aea_30m=raster("~/AEA_stack/bio1_dem_aea_WGS84.tif")
permafrost_map_30m=raster("~/AEA_stack/permafrost_dem_aea_WGS84.tif")
hd_dem_crop_aea=raster("~/AEA_stack/hengduan_dem_aea_WGS84.tif")

#== prepare present-day temperature and permafrost dataframe
stack_modern=stack(list(bio1_sp_aea_30m, permafrost_map_30m))
stack_modern_df=as.data.frame(stack_modern)
# remove NaN
stack_modern_df_NaN=stack_modern_df[complete.cases(stack_modern_df), ]
# names
names(stack_modern_df_NaN)=c("maat", "Perma_Distr_map")

#== build a model, here an example with glm
m_glm=glm(Perma_Distr_map~., data = stack_modern_df_NaN)
#summary(m_glm)
#Call:
#  glm(formula = Perma_Distr_map ~ ., data = stack_modern_df_NaN)

#Deviance Residuals: 
#  Min        1Q    Median        3Q       Max  
#-0.42652  -0.15359  -0.07977   0.03921   0.97359  

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.697e-01  4.030e-05    4210   <2e-16 ***
#  maat        -4.699e-02  1.357e-05   -3463   <2e-16 ***
---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#(Dispersion parameter for gaussian family taken to be 0.06179602)
#Null deviance: 3951644  on 51952356  degrees of freedom
#Residual deviance: 3210449  on 51952355  degrees of freedom
#AIC: 2803395
#Number of Fisher Scoring iterations: 2
  
#============= present-day m_glm end: y=1.697e-01 - -4.699e-02x with p-value <2e-16 *** ===========

#============= predict past permafrost within Naleng's catchment based on m_glm =============

#== prepare dem of Naleng's catchment
# load Naleng Shapefile and reprojection
nalengcatch = readOGR("Naleng-DEM/Naleng-DEM.shp")
nalengcatch_aea=spTransform(nalengcatch, crs(hd_dem_crop_aea))
# load Naleng dem and cropped 
nalengdem_aea_30m=raster("nalengdem_aea_30m.tif")
nalengdem_aea_30m=raster::mask(raster::crop(nalengdem_aea_30m, raster::extent(nalengcatch_aea)), nalengcatch_aea)

#== bio1 of Naleng
naleng_base_t=raster::mask(raster::crop(bio1_sp_aea_30m, raster::extent(nalengcatch_aea)), nalengcatch_aea)
#plot(naleng_base_t)
#title("BIO1")

#== add paleo anomaly
naleng_temp_adj=list()
for (i in 1:dim(paleo_temp)[1]) {
  print(i)
  anomaly=paleo_temp[i, 2]
  buff_temp_adj=naleng_base_t+anomaly
  names(buff_temp_adj)=paste0("yearBP_", paleo_temp[i, 1])
  naleng_temp_adj[[i]]=buff_temp_adj
}
stackr = stack(naleng_temp_adj)
plot(stackr)
# palaeo temp and remove NA
stack_paleo=as.data.frame(stackr, xy=TRUE)
stack_paleo_df=stack_paleo[complete.cases(stack_paleo), ]

#== predict
meta_paleo_glm_mean=stack_paleo_df[c("x", "y")]
for (i in 3:dim(stack_paleo_df)[2]) {
  ka_layer=colnames(stack_paleo_df[i])
  print(ka_layer)
  #
  sub_ka=stack_paleo_df[c(ka_layer)]
  names(sub_ka)=c("maat")
  paleo_glm=predict.glm(m_glm, sub_ka, type="response", se.fit = T)
  #
  df_paleo_glm=data.frame(mean=paleo_glm$fit, CI95=paleo_glm$se.fit, ka_layer=ka_layer)
  write.csv(df_paleo_glm, paste0("predict/", ka_layer, ".csv"))
  #
  fit_mean=as.data.frame(df_paleo_glm[, 1])
  names(fit_mean)=paste0("glm_permafrost.", ka_layer)
  meta_paleo_glm_mean=cbind(meta_paleo_glm_mean, fit_mean)
}
#== take the coordinates and name the predicted layers
meta_paleo_glm_mean=stack_paleo_df[c("x", "y")] # x and y are coordinates of grids
for (i in 1:dim(paleo_temp)[1]) {
  print(i)
  ka_layer=paste0("yearBP_", paleo_temp[i, 1])
  df_paleo_glm=read.csv(paste0("predict/", ka_layer, ".csv"), row.names = 1)
  fit_mean=as.data.frame(df_paleo_glm[, 1])
  names(fit_mean)=paste0("glm_permafrost.", ka_layer)
  meta_paleo_glm_mean=cbind(meta_paleo_glm_mean, fit_mean)
}
paleo_perma=meta_paleo_glm_mean  
# convert dataframe to stack
coordinates(paleo_perma) = ~ x + y
gridded(paleo_perma) = TRUE
stack_paleo_perma=stack(paleo_perma)
crs(stack_paleo_perma)="+proj=aea +lat_0=30 +lon_0=95 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
#plot(stack_paleo_perma, maxnl=19)
#writeRaster(stack_paleo_perma, filename = "stack_paleo_perma_30m_glm_22-0ka.tif", overwrite = TRUE)

#============ predict palaeo permafrost within catchment end ==============

#======= find the threshold of permafrost presence by comparing preidcted permafrost at 0 ka and presenty-day

#== load Paleo_permafrost_glm
stack_paleo_perma_30m_glm=stack("stack_paleo_perma_30m_glm_22-0ka.tif")

#== name layers 
years=seq(0, 22000, 500)
for (i in 1:45) {
  print(names(stack_paleo_perma_30m_glm)[i])
  names(stack_paleo_perma_30m_glm)[i] = paste0("year_", years[i])
  print(names(stack_paleo_perma_30m_glm)[i])
}

#== reverse layers 22-0 ka
stack_perma_rev=stack_paleo_perma_30m_glm[[45]]
for (i in 44:1) {
  print(i)
  stack_perma_rev=stack(stack_perma_rev,stack_paleo_perma_30m_glm[[i]])
}
names(stack_perma_rev)
#plot(stack_perma_rev)

#== prepare 0 ka and present-day
# 0 ka
ka0=stack_perma_rev[[45]]
# present-day
nalengdem_permafrost_map_30m=raster::mask(raster::crop(permafrost_map_30m, raster::extent(nalengcatch_aea)), nalengcatch_aea)
plot(nalengdem_permafrost_map_30m)
title("Permafrost map within Naleng catchment, v.2017", cex.main = 1)
modern=nalengdem_permafrost_map_30m
# stack
ka0_modern=stack(list(modern, ka0))
ka0_modern_df=as.data.frame(ka0_modern, xy=TRUE)
ka0_modern_df=ka0_modern_df[complete.cases(ka0_modern_df), ]

#== y~x
vslm=lm(year_0~permafrost_dem_aea_WGS84, ka0_modern_df)
anova(vslm)
# y = 0.1579 + 0.0646x
# plot
ggplot(data = ka0_modern_df, aes(x = permafrost_dem_aea_WGS84, y = year_0)) +
  scale_y_continuous(breaks = seq(0.05, 0.25, 0.0005)) +
  geom_point(size=2) +
  geom_quantile(quantiles = 0.5) #Median
# taking x = 0.9694 

#== check predicted permafrost of ka 0
vlimit=0.22052
Fplus_map=ka0
values(Fplus_map)[values(Fplus_map) < vlimit] = NA
#
p_Fplus_map=gplot(Fplus_map) + 
  scale_x_continuous(name = "Lon") +
  scale_y_continuous(name = "Lat") +
  geom_tile(aes(fill = value)) +
  facet_wrap(~ variable) +
  geom_polygon(data = nalengcatch_aea, aes(long, lat), fill = NA, color = "black", size = 0.25) +
  scale_fill_gradientn(colours = rev(terrain.colors(225))) +
  coord_equal() +
  theme(axis.title = element_text(size = 5),
        axis.text = element_text(size = 4),
        title = element_text(size = 5),
        legend.text = element_text(size = 5),
        strip.text = element_text(size = 4.5))
p_Fplus_map

#== cut-off all paleo layers
vlimit=vlimit
stack_perma_rev_adj=stack_perma_rev
values(stack_perma_rev_adj)[values(stack_perma_rev_adj) < vlimit] = 0
values(stack_perma_rev_adj)[values(stack_perma_rev_adj) >= vlimit] = 1
plot(stack_perma_rev_adj)
writeRaster(stack_perma_rev_adj, filename = paste0("stack_perma_rev_adj_v", vlimit, ".tif"), overwrite = TRUE)
#======== end ===========
          






