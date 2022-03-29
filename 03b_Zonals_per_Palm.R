# Bastiaen Boekelo, March 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Create zonal statistics per tree and write .csv's to 2_Intermediate folder

##########################################
#### CHECK WD's & ADAPT IF NECESSARY #####
##########################################


rm(list=ls()) 

# Set Script wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

library(rgdal)
library(raster)
library(MASS)    # The following objects are masked from 'package:raster': area, select
source("functions/global_variables.R")
source("functions/03_zonal_stats_tree.R")

buffers   <- c(0.5, 1, 1.5, 2, 2.5)

all_palms <- readOGR(paste0(wd, "2_Intermediate"), "03_palm_centerpoints")
combis    <- unique(all_palms@data$PR_FI)
combis    <- combis[grepl("CK", combis)] # Remove CK for mismatch


for(COMBI in 1:length(combis)){
  
  PR_FI     <- combis[COMBI]
  
  inputs    <- list.files(paste0(wd, "2_Intermediate"), pattern=PR_FI)
  
  if(length(inputs) > 0){
    ras_names   <- grep('_Mask.tif$', inputs, value=TRUE)
    ortho_name  <- grep('Ortho', ras_names, value=TRUE)
    NDVI_name   <- grep('NDVI', ras_names, value=TRUE)
    NDRE_name   <- grep('NDRE', ras_names, value=TRUE)
    
    if(length(ortho_name) == 1 & length(NDVI_name) == 1 & length(NDRE_name) == 1){
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
      
      # Input palm points
      palms_F      <- all_palms[all_palms@data$PR_FI == PR_FI,]
      palms_F      <- spTransform(palms_F, crs(ortho) )
      
      ################
      ### ANALYSIS ###
      ################
      
      print("Zonal statistics")
      
      # Zonal statistics
      for(PALM_ID in 1:nrow(palms_F)){
        
        palm <- palms_F[PALM_ID,]
        palm_name <- palms_F@data$TRT_ID[PALM_ID]
        
        for(BUF in 1:length(buffers)){
          
          for(BRK in 1:(length(brks_all)-1)){
            i <- 0
            
            mask_ras <- round( avg_pos  / cellStats(avg_pos,  "max", na.rm=T) * 10000)
            
            brk1 <- brks_all[BRK]
            brk2 <- brks_all[BRK+1]
            
            brks      <- c(brk1, brk2)
            brk_vals  <- as.vector(quantile(mask_ras, brks))
            
            mask_ras[mask_ras <- mask_ras < brk_vals[1] | mask_ras >= brk_vals[2] ] <- NA
            mask_ras[mask_ras > 0] <- 1
            # mask_ras@file@datanotation    <- NDVI@file@datanotation
            
            NDRE_mask   <- mask(NDRE, mask_ras)
            NDVI_mask   <- mask(NDVI, mask_ras)
            
            buf_dist <- buffers[BUF]
            
            print(paste0(PR_FI, " ", as.character(palm@data$TRT_ID), " ", buf_dist, "m_brk_", brk2))
            zonals_NDRE_1palm_1buf_1brk   <- zonal_stats_tree(PR_FI, NDRE_mask, palm, buffer=buf_dist, brks, palm_name)
            zonals_NDVI_1palm_1buf_1brk   <- zonal_stats_tree(PR_FI, NDVI_mask, palm, buffer=buf_dist, brks, palm_name)
            zonals_1palm_1buf_1brk        <- rbind(zonals_NDRE_1palm_1buf_1brk, zonals_NDVI_1palm_1buf_1brk)
            
            if(BRK == 1){
              zonals_1palm_1buf <- zonals_1palm_1buf_1brk
            } else {
              zonals_1palm_1buf <- rbind(zonals_1palm_1buf, zonals_1palm_1buf_1brk)
            }
          }
          if(BUF == 1){
            zonals_all <- zonals_1palm_1buf
          } else {
            zonals_all <- rbind(zonals_all, zonals_1palm_1buf)
          }
        }
        write.csv(zonals_all, paste0(wd, "2_Intermediate/03_ZonalStatsPerPalm_", PR_FI, "__", palm_name, ".csv"), row.names=F)
      }
    } else {
      
      print(paste0(PR_FI, ": not complete?"))
      
    } 
    
  } else {
    print(paste0(PR_FI, ": no input at all"))
  }
}






