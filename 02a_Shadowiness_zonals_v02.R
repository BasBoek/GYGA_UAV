# Bastiaen Boekelo, March 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Create zonal statistics before visualization to inspect how shadowiness and shininess influence the VI response.

rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

library(rgdal)
library(raster)
library(MASS)    # The following objects are masked from 'package:raster': area, select
source("functions/zonal_stats_field.R")
source("functions/global_variables.R")


for(COMBI in 1:length(combis)){
  
  PR_FI     <- combis[COMBI]
  
  inputs    <- list.files(paste0(wd, "2_Intermediate"), pattern=PR_FI)
  
  if(length(inputs) > 0){
    ras_names   <- grep('_Mask.tif$', inputs, value=TRUE)
    points_name <- grep('.shp$', inputs, value=TRUE)
    
    ortho_name  <- grep('Ortho', ras_names, value=TRUE)
    NDVI_name   <- grep('NDVI', ras_names, value=TRUE)
    NDRE_name   <- grep('NDRE', ras_names, value=TRUE)
    
    if(length(ortho_name) == 1 & length(NDVI_name) == 1 & length(NDRE_name) == 1 & length(points_name) == 1 ){
      print("preprocess rasters")
      
      # Input rasters
      ortho   <- stack(paste0(wd, "2_Intermediate/", ortho_name))
      NDRE    <- raster(paste0(wd, "2_Intermediate/", NDRE_name))
      NDVI    <- raster(paste0(wd, "2_Intermediate/", NDVI_name))

      # Preprocess mask raster
      
      red_norm  <- ( ortho[[1]] - cellStats(ortho[[1]], "mean", na.rm=T) ) / cellStats(ortho[[1]], "sd", na.rm=T)
      NIR_norm  <- ( ortho[[5]] - cellStats(ortho[[5]], "mean", na.rm=T) ) / cellStats(ortho[[5]], "sd", na.rm=T)
      refl_avg  <- ( red_norm + NIR_norm ) / 2
      avg_pos   <- refl_avg - cellStats(refl_avg, "min", na.rm=T)

      # # PREVENTING SCRIPT FROM CRASHING:
      # NDRE@file@blockrows <- NDVI@file@blockrows
      # NDRE@file@blockcols <- NDVI@file@blockcols
      # 
      # avg_pos@file@blockrows <- NDVI@file@blockrows
      # avg_pos@file@blockcols <- NDVI@file@blockcols
      # avg_pos@file@driver    <- NDVI@file@driver
      
      # Input palm points
      points      <- readOGR(paste0(wd, "2_Intermediate"), substr(points_name, 1, nchar(points_name)-4) )
      points      <- spTransform(points, crs(ortho) )
      
      ################
      ### ANALYSIS ###
      ################
      
      print("Zonal statistics")
      
      # Zonal statistics
      for(BUF in 1:length(buffers)){
        
        buf_dist <- buffers[BUF]
        
        print(paste0(PR_FI, ": NDRE, ", buf_dist, "m"))
        zonals_NDRE   <- zonal_stats_field(PR_FI, NDRE, avg_pos, points, buffer=buf_dist, brks)
        write.csv(zonals_NDRE, paste0(wd, "2_Intermediate/02_ZonalStats_", PR_FI, "_NDRE_", buf_dist, "m.csv"), row.names=F)
        
        print(paste0(PR_FI, ": NDVI, ", buf_dist, "m"))
        zonals_NDVI   <- zonal_stats_field(PR_FI, NDVI, avg_pos, points, buffer=buf_dist, brks)
        write.csv(zonals_NDVI, paste0(wd, "2_Intermediate/02_ZonalStats_", PR_FI, "_NDVI_", buf_dist, "m.csv"), row.names=F)
        
      }
      
    } else {
      
      print(paste0(PR_FI, ": not complete?"))
      
    } 
    
  } else {
    print(paste0(PR_FI, ": no input at all"))
  }
}




