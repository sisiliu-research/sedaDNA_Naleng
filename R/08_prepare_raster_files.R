# by Sisi Liu: Sisi Liu: sisi.liu@awi.de | sisi.liu.research@gmail.com (permanent address)
# R version 4.2.2 (2022-10-31)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.1

library(raster)
library(sp)
library(rgdal)
library(sf)
library(raster)
library(dismo)

#== make extent/lake Naleng crs
lat = 31.11333300
lon = 99.75516700
byextent=extent(c(lon-0.5, lon+0.5, lat-0.5, lat+0.5))

#== load data
bio1=raster("raster/wc2.1_30s_bio_1-WGS84-res-bio.tif")
merged_dem=raster("raster/merged_dem.tif")
permafrost_map=raster("raster/Perma_Distr_map.tif")

#== check
crs(bio1)
crs(merged_dem)
crs(permafrost_map)

#== crop BIO1
t_bio1 = crop(bio1, byextent)
t_bio1_sp = t_bio1/10 
t_bio1_sp #crs: +proj=longlat +datum=WGS84 +no_defs

#== crop SRTM 30 m
merged_dem_crop=crop(merged_dem, byextent) #crs: +proj=longlat +datum=WGS84 +no_defs

#== re-projection to crs = 102025|"+proj=aea +lat_0=30 +lon_0=95 +lat_1=15 +lat_2=65 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
bio1_sp_aea=projectRaster(t_bio1_sp, crs = CRS(SRS_string = "ESRI:102025"))
merged_dem_aea=projectRaster(merged_dem_crop, crs = CRS(SRS_string = "ESRI:102025"))

#== resample to 30 m resolution
bio1_sp_aea_30m=resample(bio1_sp_aea, merged_dem_aea, method='bilinear')
permafrost_map_30m=resample(permafrost_map, merged_dem_aea, method='bilinear')

#== save/the three raster files are provided under raster/
writeRaster(bio1_sp_aea_30m, filename = "bio1_dem_aea_WGS84.tif", overwrite = TRUE)
writeRaster(merged_dem_aea, filename = "hengduan_dem_aea_WGS84.tif", overwrite = TRUE)
writeRaster(permafrost_map_30m, filename = "permafrost_dem_aea_WGS84.tif", overwrite = TRUE)
#======== end =======








