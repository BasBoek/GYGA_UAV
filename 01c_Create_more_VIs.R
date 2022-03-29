# Bastiaen Boekelo, April 2021
# Nebraska Project
# Calculate different VIs

library(raster)
library(rgdal)

rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

PR_FI_files <- list.files(paste0(wd, "2_Intermediate/01b_singlebands"), pattern = ".tif$")
fields      <- unique(substr(list.files(paste0(wd, "2_Intermediate/01b_singlebands"), pattern = ".tif$"), 10, 14))


PR_FI <- "CK_F7"
for(PR_FI in fields){
  
  PR_FI_sel <- PR_FI_files[grepl(PR_FI, PR_FI_files)]
  
  RED    <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[3]))
  NIR    <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[5]))
  RE     <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[4]))
  GREEN  <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[2]))
  

  print(PR_FI)
  # BLUE   <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[1]))
  
  # WDRVI  <- round((0.1 * NIR - RED) / (0.1 * NIR + RED) * 1000)
  # MSAVI  <- round(( (NIR - RED) * 1.5 ) / ( NIR + RED + 1.5 ) * 1000) # using L = 0.5 as compromise (https://wiki.landscapetoolbox.org/doku.php/remote_sensing_methods:modified_soil-adjusted_vegetation_index)
  # NDVIG  <- round((NIR - RED) / (NIR + RED) * log(1 + GREEN / 32768 * 100) * 1000 - 1000) # With sqrt better differentiation, but with log less dependent on total brightness of image
  # NDREG  <- round((RE  - RED) / (RE  + RED) * log(1 + GREEN / 32768 * 100) * 1000 - 1000) # With sqrt better differentiation, but with log less dependent on total brightness of image
  # NDVIG  <- round((NIR - RED) / (NIR + RED) * log10(10 + GREEN / 32768 * 100) * 1000 - 1000) # With sqrt better differentiation, but with log less dependent on total brightness of image
  # NDREG  <- round((RE  - RED) / (RE  + RED) * log10(10 + GREEN / 32768 * 100) * 1000 - 1000) + 400 # With sqrt better differentiation, but with log less dependent on total brightness of image
  
  
  # Not including MSAVI2 since it is very similar to MSAVI and seems to differentiate less between palm trees and soil than MSAVI 
  #MSAVI2 <- (  2 * NIR + 1 - sqrt( (2 * NIR + 1)^2 - 8*(NIR - RED))  ) / 2
  
  # NDVIGC <- NDVIG
  # NDVIGC[NDVIGC < 0] <- NA
  # NDREGC <- NDREG
  # NDREGC[NDREGC < 0] <- NA
  # 

  # writeRaster(WDRVI, paste0(wd, "2_Intermediate/01_WDRVI_", PR_FI, "_Mask.tif"), format= "GTiff", datatype= 'INT2S', overwrite=T)
  # writeRaster(MSAVI, paste0(wd, "2_Intermediate/01_MSAVI_", PR_FI, "_Mask.tif"), format= "GTiff", datatype= 'INT2S', overwrite=T)
  # writeRaster(NDVIG, paste0(wd, "2_Intermediate/01_NDVIG_", PR_FI, "_Mask.tif"), format= "GTiff", datatype= 'INT2S', overwrite=T)
  # writeRaster(NDREG, paste0(wd, "2_Intermediate/01_NDREG_", PR_FI, "_Mask.tif"), format= "GTiff", datatype= 'INT2S', overwrite=T)
  # writeRaster(NDVIGC, paste0(wd, "2_Intermediate/01_NDVIGC_", PR_FI, "_Mask.tif"), format= "GTiff", datatype= 'INT2S', overwrite=T)
  # writeRaster(NDREGC, paste0(wd, "2_Intermediate/01_NDREGC_", PR_FI, "_Mask.tif"), format= "GTiff", datatype= 'INT2S', overwrite=T)
  
}

# Calculate the overall weighed (by sd) averages 

greens <- c()
nirs   <- c()
reds   <- c()
greens_sd <- c()
nirs_sd   <- c()
reds_sd   <- c()

fields_temp <- fields[fields != "CK_F3" & fields != "CK_F5"] # Remove these fields for this analysis because weirdly high reflectances

for(i in 1:length(fields_temp)){
  
  PR_FI <- fields_temp[i]
  print(PR_FI)
  PR_FI_sel <- PR_FI_files[grepl(PR_FI, PR_FI_files)]
  
  RED    <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[3]))
  NIR    <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[5]))
  RE     <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[4]))
  GREEN  <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[2]))
  
  greens[i] <- round(cellStats(GREEN,   "mean", na.rm=T) / 32768 * 100,2)
  nirs[i]    <- round(cellStats(NIR,    "mean", na.rm=T) / 32768 * 100,2)
  reds[i]    <- round(cellStats(RED,    "mean", na.rm=T) / 32768 * 100,2)
  
  greens_sd[i]    <- round(cellStats(GREEN,    "sd", na.rm=T) / 32768 * 100,2)
  nirs_sd[i]    <- round(cellStats(NIR,    "sd", na.rm=T) / 32768 * 100,2)
  reds_sd[i]    <- round(cellStats(RED,    "sd", na.rm=T) / 32768 * 100,2)
  
}
weighted.mean(greens, 1/greens_sd)
weighted.mean(nirs, 1/nirs_sd)
weighted.mean(reds, 1/reds_sd)

mean(greens_sd)

# So if for every sd deviating from mean the weigh should be * 0.33
# 

# Green weighing proposition: (sd_green = 2, so if alpha=1, then with every 2% reflection deviating from 4.7, 
# The weigh will be (1 + 2 =) 3 times less
# round(1 / (1 + (green_pix - 4.7)^2 ) * 100)


all_palms <- readOGR(paste0(wd, "2_Intermediate"), "03_palm_centerpoints")
PR_FI <- "SS_F7"
palm <- all_palms[all_palms@data$PR_FI == PR_FI &all_palms@data$TRT_ID == "BMP_12",]

#pdf("test_4.7.pdf")
for(i in 1:length(fields)){
  
  # PR_FI <- fields[i]
  print(PR_FI)
  
  PR_FI_sel <- PR_FI_files[grepl(PR_FI, PR_FI_files)]

  GREEN   <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[2])) / 32768 * 100
  NIR     <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[5])) / 32768 * 100
  RED     <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[3])) / 32768 * 100
  RE      <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[4])) / 32768 * 100
    
  palm    <- spTransform(palm, crs(NIR) )
  palmbuf <- buffer(palm, 4)
  
  WDRVI        <- round((0.1 * NIR - RED) / (0.1 * NIR + RED) * 1000)
  WDRVI_mask   <- crop(WDRVI, palmbuf)
  NDVI         <- (NIR - RED) / (NIR + RED)
  NDVI_mask    <- crop(NDVI,  palmbuf)
  NDRE         <- (RE - RED) / (RE + RED)
  NDRE_mask    <- crop(NDRE,  palmbuf)
  GREEN_mask   <- crop(GREEN, palmbuf)
  #GREEN_weight <- (cellStats(GREEN_mask, "sd", na.rm=T)*1) / ( (cellStats(GREEN_mask, "sd", na.rm=T)*1) + (GREEN_mask - (cellStats(GREEN_mask, "mean", na.rm=T) - 0.5 * cellStats(GREEN_mask, "sd", na.rm=T)))^2 ) 
  #GREEN_weight <- sqrt(1 / ( 1 + (GREEN_mask - (cellStats(GREEN_mask, "mean", na.rm=T) + 0 * cellStats(GREEN_mask, "sd", na.rm=T)))^2 ))
  GREEN_weight <- (1 / ( 1 + (GREEN_mask - (cellStats(GREEN_mask, "mean", na.rm=T) ))^2 ))^0.75
  GREEN_weight <- GREEN_weight / cellStats(GREEN_weight, "sum", na.rm=T)
  plot(GREEN_weight)
  
  plot(WDRVI_mask)
  WDRVI_weighed <- WDRVI_mask * GREEN_weight
  cellStats(WDRVI_weighed, "sum", na.rm=T) 
  cellStats(WDRVI_mask, "mean", na.rm=T) 
  
  plot(NDVI_mask)
  NDVI_weighed <- NDVI_mask * GREEN_weight
  cellStats(NDVI_weighed, "sum", na.rm=T) 
  cellStats(NDVI_mask, "mean", na.rm=T) 
  
  plot(NDRE_mask)
  NDRE_weighed <- NDRE_mask * GREEN_weight
  cellStats(NDRE_weighed, "sum", na.rm=T) 
  cellStats(NDRE_mask, "mean", na.rm=T) 

  
  ds
}
# Conclusions:
# Green weighing does not matter for NDVI and NDRE, but it certainly does for WDRVI
# Probably more proper selection of relevant pixels. Parameters (0.75 and the 1's might be tweaked, but let's first see how this goes)





