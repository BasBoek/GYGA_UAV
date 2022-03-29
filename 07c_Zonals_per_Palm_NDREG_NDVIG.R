# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Create zonal statistics per tree

rm(list=ls())    # Clean script, clean head.

# Set Script wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

library(rgdal)
library(raster)
library(MASS)    # The following objects are masked from 'package:raster': area, select
source("functions/global_variables.R")
source("functions/07_zonal_stats_tree.R")

# First number is break for shadows based on mean of RGB
# Second number is break for shines based on NIR
# DL_brks <- list( c(0.00, 2.00) * 32768,
#                  c(0.00, 0.80) * 32768,
#                  c(0.00, 0.70) * 32768,
#                  c(0.01, 2.00) * 32768,
#                  c(0.01, 0.80) * 32768,
#                  c(0.01, 0.70) * 32768,
#                  c(0.02, 2.00) * 32768,
#                  c(0.02, 0.80) * 32768,
#                  c(0.02, 0.70) * 32768)

DL_brks <- list( c(0.00, 2.00) * 32768,  # Nothing removed
                 c(0.00, 0.70) * 32768,  # Shines removed
                 c(0.02, 2.00) * 32768,  # Shadows removed
                 c(0.02, 0.70) * 32768)  # Shines and shadows removed


buffers   <- c(0.5, 1, 1.5, 2.5, 3.5)

VIs       <- c("NDREG", "NDVIG", "NDREGC", "NDVIGC")
nr_VIs    <- length(VIs)

all_palms <- readOGR(paste0(wd, "2_Intermediate"), "03_palm_centerpoints")
nr_rows   <- length(all_palms) * length(DL_brks) * length(buffers) * nr_VIs # nr_palms * nr_breaks * nr_buffers * nr_VIs
combis    <- unique(all_palms@data$PR_FI)



combis    <- combis[grepl("JB_F5", combis)]



# COMBI   <- 1
# PALM_ID <- 1
# BUF     <- 1
# PALM_ID <- 1
# BRK_SET <- 1

for(COMBI in 1:length(combis)){
  
  PR_FI     <- combis[COMBI]
  
  inputs    <- list.files(paste0(wd, "2_Intermediate"), pattern=PR_FI)
  
  i <- 0
  df <- data.frame(PR_FI = character(), 
                   palm_id = character(),
                   VI = character(),
                   Buffer_m = character(), 
                   Break_cat = numeric(), 
                   pix_nr = numeric(), 
                   pix_mean = numeric(), 
                   pix_sd = numeric())
  df[1:nr_rows,] <- NA
  
  
  if(length(inputs) > 0){
    ras_names     <- grep('_Mask.tif$', inputs, value=TRUE)
    ortho_name    <- grep('Ortho', ras_names, value=TRUE)
    NDREG_name    <- grep('NDREG_', ras_names, value=TRUE)
    NDVIG_name    <- grep('NDVIG_', ras_names, value=TRUE)
    NDREGC_name   <- grep('NDREGC', ras_names, value=TRUE)
    NDVIGC_name   <- grep('NDVIGC', ras_names, value=TRUE)
    
    # if(length(ortho_name) == 1 & length(NDVI_name) == 1 & length(NDRE_name) == 1){
    if(length(ortho_name) == 1 & length(NDREG_name) == 1 & length(NDREGC_name) == 1 & length(NDVIG_name) == 1 & length(NDVIGC_name) == 1){
      print("preprocess rasters")
      
      # Input rasters
      ortho    <- stack(paste0(wd, "2_Intermediate/", ortho_name))
      NDREG    <- raster(paste0(wd, "2_Intermediate/", NDREG_name))
      NDVIG    <- raster(paste0(wd, "2_Intermediate/", NDVIG_name))
      NDREGC   <- raster(paste0(wd, "2_Intermediate/", NDREGC_name))
      NDVIGC   <- raster(paste0(wd, "2_Intermediate/", NDVIGC_name))
      
      # Create raster masks based on NDVI > 7250 and 
      
      RGB_avg <- (ortho[[1]] + ortho[[2]] + ortho[[3]]) / 3
      
      max_RGB <- cellStats(RGB_avg,    "max", na.rm=T)
      min_NIR <- cellStats(ortho[[5]], "min", na.rm=T)
      
      # Input palm points
      palms_F      <- all_palms[all_palms@data$PR_FI == PR_FI,]
      palms_F      <- spTransform(palms_F, crs(ortho) )
      
      ################
      ### ANALYSIS ###
      ################
      
      print("Zonal statistics")
      
      # Zonal statistics
      for(BUF in 1:length(buffers)){
        buf_dist <- buffers[BUF]
        
        
        for(PALM_ID in 1:nrow(palms_F)){
          palm      <- palms_F[PALM_ID,]
          palm_name <- palms_F@data$TRT_ID[PALM_ID]
          
          for(BRK_SET in 1:length(DL_brks)){
            
            brk_D <- DL_brks[[BRK_SET]][1]
            brk_L <- DL_brks[[BRK_SET]][2]
            
            if(brk_D < max_RGB  &  brk_L > min_NIR){
              
              mask1       <- RGB_avg
              mask1[mask1 < brk_D] <- NA
              
              mask2       <- ortho[[5]] # NIR band
              mask2[mask2 > brk_L] <- NA
              
              mask_ras <- mask1 + mask2
              mask_ras[mask_ras > 0] <- 1
              
              print(paste0(PR_FI, " ", as.character(palm@data$TRT_ID), " ", buf_dist, "m_brk_", BRK_SET))
              
              NDREG_mask                          <- mask(NDREG, mask_ras)
              zonals_NDREG_1palm_1buf_1brk        <- zonal_stats_tree(PR_FI, NDREG_mask, palm, buf_dist, BRK_SET, palm_name)
              i <- i + 1
              df[i,] <- zonals_NDREG_1palm_1buf_1brk
              
              NDREGC_mask                          <- mask(NDREGC, mask_ras)
              zonals_NDREGC_1palm_1buf_1brk        <- zonal_stats_tree(PR_FI, NDREGC_mask, palm, buf_dist, BRK_SET, palm_name)
              i <- i + 1
              df[i,] <- zonals_NDREGC_1palm_1buf_1brk
              
              NDVIG_mask                          <- mask(NDVIG, mask_ras)
              zonals_NDVIG_1palm_1buf_1brk        <- zonal_stats_tree(PR_FI, NDVIG_mask, palm, buf_dist, BRK_SET, palm_name)
              i <- i + 1
              df[i,] <- zonals_NDVIG_1palm_1buf_1brk
              
              NDVIGC_mask                          <- mask(NDVIGC, mask_ras)
              zonals_NDVIGC_1palm_1buf_1brk        <- zonal_stats_tree(PR_FI, NDVIGC_mask, palm, buf_dist, BRK_SET, palm_name)
              i <- i + 1
              df[i,] <- zonals_NDVIGC_1palm_1buf_1brk
              
            }
          }
        }
      }
    } else {
      
      print(paste0(PR_FI, ": not complete?"))
      
    } 
    
  } else {
    print(paste0(PR_FI, ": no input at all"))
  }
  df <- df[rowSums(is.na(df)) != ncol(df), ]
  write.csv(df, paste0(wd, "2_Intermediate/07_ZonalStatsPerPalm_", PR_FI, "_", paste0(VIs, collapse="_"), ".csv"), row.names=F)
  
}







