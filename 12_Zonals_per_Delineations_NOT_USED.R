
rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd  <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

library(rgdal)
library(raster)
library(MASS)    # The following objects are masked from 'package:raster': area, select
source("functions/global_variables.R")
source("functions/14_zonal_stats_tree.R")


DL_brks <- list( c(0.00, 2.00) * 32768,  # Nothing removed
                 c(0.00, 0.70) * 32768,  # Shines  removed
                 c(0.02, 2.00) * 32768,  # Shadows removed
                 c(0.02, 0.70) * 32768)  # Shines and shadows removed


PR_FI     <- "JB_F7"

tif       <- stack(paste0("../Data/1_Input/UAV/JB/2_Orthophoto_JB/5cm_Micasense_JB_F7.tif"))

all_palms <- readOGR("../Data/1_Input/UAV/JB/6_Crown_delineations_JB", layer ="Palm_Crown_Projection_Micasense_UNL_JB_F7")
all_palms <- spTransform(all_palms, crs(tif))

BLU   <- tif[[1]]
GREEN <- tif[[2]]
RED   <- tif[[3]]
RE    <- tif[[4]]
NIR   <- tif[[5]]
NDVI      <- round((NIR - RED) / (NIR + RED) * 1000)
NDRE      <- round((RE  - RED) / (RE  + RED) * 1000)
WDRVI     <- round((0.1 * NIR - RED) / (0.1 * NIR + RED) * 1000)
MSAVI     <- round(( (NIR - RED) * 1.5 ) / ( NIR + RED + 0.5 ) * 1000) # using L = 0.5 as compromise (https://wiki.landscapetoolbox.org/doku.php/remote_sensing_methods:modified_soil-adjusted_vegetation_index)
NRED      <- round(RED / (GREEN + NIR + RED) * 1000)
NGRN      <- round(GREEN / (GREEN + NIR + RED) * 1000)
NNIR      <- round(NIR / (GREEN + NIR + RED) * 1000)

# Green Weight
GRN_WGHT  <- GREEN
GRN_WGHT[GRN_WGHT > quantile(GRN_WGHT, 0.95)] <- 1
GRN_WGHT  <- sqrt(GRN_WGHT)
GRN_WGHT[GRN_WGHT > 250] <- NA
GRN_WGHT  <- GRN_WGHT / cellStats(GRN_WGHT,  "sum")

NDVI_palm <- mask(NDVI, all_palms[1,])

# (1) Derive mean values 
NDVI_avg   <- cellStats(NDVI,    "mean")
NDRE_avg   <- cellStats(NDRE,    "mean")
WDRVI_avg  <- cellStats(WDRVI,   "mean")
MSAVI_avg  <- cellStats(MSAVI,   "mean")
NRED_avg   <- cellStats(NRED,    "mean")
NGRN_avg   <- cellStats(NGRN,    "mean")
NNIR_avg   <- cellStats(NNIR,    "mean")

NDVI_GW    <- cellStats(NDVI  * GRN_WGHT,  "sum")
NDRE_GW    <- cellStats(NDRE  * GRN_WGHT,  "sum")
WDRVI_GW   <- cellStats(WDRVI * GRN_WGHT,  "sum")
MSAVI_GW   <- cellStats(MSAVI * GRN_WGHT,  "sum")
NRED_GW    <- cellStats(NRED  * GRN_WGHT,  "sum")
NGRN_GW    <- cellStats(NGRN  * GRN_WGHT,  "sum")
NNIR_GW    <- cellStats(NNIR  * GRN_WGHT,  "sum")

nr_rows   <- length(all_palms) * length(DL_brks)

VIs       <- c("NDVI", "NDRE", "NDVIG", "NDREG", "NDVIGC", "NDREGC", "MSAVI", "WDRVI")

COMBI     <- 1
PALM_ID   <- 1

# Create function looping over input data,  

for(VI_index in 1:length(VIs)){
  str_VI <- VIs[VI_index]
  
  
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
    ras_names   <- grep('_Mask.tif$', inputs, value=TRUE)
    ortho_name  <- grep('_Ortho_', ras_names, value=TRUE)
    VI_name     <- grep(paste0("_", str_VI, "_"), ras_names, value=TRUE)
    
    if(length(ortho_name) == 1 & length(VI_name) == 1){
      print("preprocess rasters")
      
      # Input rasters
      ortho <- stack(paste0(wd, "2_Intermediate/", ortho_name))
      VI    <- raster(paste0(wd, "2_Intermediate/", VI_name))
      
      # Create raster masks based on NDVI > 7250 and 
      
      RGB_avg <- (ortho[[1]] + ortho[[2]] + ortho[[3]]) / 3
      GREEN   <- ortho[[2]]
      
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
  
  write.csv(df, paste0(wd, "2_Intermediate/10_ZonalStatsPerPalm_", PR_FI, "_", str_VI, ".csv"), row.names=F)
  
  
  
}

