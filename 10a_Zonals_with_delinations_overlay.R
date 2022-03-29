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
VIs       <- c("NDVI", "NDRE")

all_palms <- readOGR(paste0(wd, "2_Intermediate"), "03_palm_centerpoints")


nr_rows   <- length(all_palms) * length(DL_brks) * length(buffers) * length(VIs) # nr_palms * nr_breaks * nr_buffers * nr_VIs

PALM_ID   <- 1
PR_FI    <- "CK_F4"
VI_index <- 1

create_zonals <- function(PR_FI){
  
  inputs <- list.files(paste0(wd, "2_Intermediate/01a2_Crowns"), pattern=PR_FI)
  
  for(VI_index in 1:length(VIs)){
    str_VI <- VIs[VI_index]
    
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
      ras_names   <- grep('_Crown.tif$', inputs, value=TRUE)
      ortho_name  <- grep('_Ortho_', ras_names, value=TRUE)
      VI_name     <- grep(paste0("_", str_VI, "_"), ras_names, value=TRUE)
      
      if(length(ortho_name) == 1 & length(VI_name) == 1){
        print("preprocess rasters")
        
        # Input rasters
        ortho <- stack(paste0(wd, "2_Intermediate/01a2_Crowns/", ortho_name))
        VI    <- raster(paste0(wd, "2_Intermediate/01a2_Crowns/", VI_name))

        # Create raster masks based on NDVI > 7250 and 
        
        RGB_avg <- (ortho[[1]] + ortho[[2]] + ortho[[3]]) / 3
        
        max_RGB <- cellStats(RGB_avg,    "max", na.rm=T)
        min_NIR <- cellStats(ortho[[5]], "min", na.rm=T)
        
        # Input palm points
        palms_F      <- all_palms[all_palms@data$PR_FI == PR_FI,]
        palms_F      <- spTransform(palms_F, crs(ortho) )
        plot(palms_F)
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
                
                print(paste0(PR_FI, " ", str_VI, " ", as.character(palm@data$TRT_ID), " ", buf_dist, "m_brk_", BRK_SET))
                
                VI_mask                      <- mask(VI, mask_ras)
                palm_zonals  <- zonal_stats_tree(PR_FI, VI_mask, palm, buf_dist, BRK_SET, palm_name)
                i <- i + 1
                df[i,] <- palm_zonals
                
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
    
    write.csv(df, paste0(wd, "2_Intermediate/10_ZonalStatsPerCrown_", PR_FI, "_", str_VI, ".csv"), row.names=F)
    
  }
}

#create_zonals("CK_F4")
create_zonals("RI_F4")






