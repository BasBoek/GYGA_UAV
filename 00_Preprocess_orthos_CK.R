# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Write single band rasters to intermediate folder


##########################################
#### CHECK WD's & ADAPT IF NECESSARY #####
##########################################

      
# Set Script and Data wd
rm(list=ls())  
   
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

library(raster)

ras_names <- list.files(paste0(wd, "1_Input/UAV/CK/2b_Reflectance_CK"), pattern=".tif$", full.names = F)

for(i in 1:length(ras_names)){
  ras <- stack(paste0(wd, "1_Input/UAV/CK/2b_Reflectance_CK/", ras_names[i]))
  ras <- round(ras * 32768) # Comparable with original format, see Micasense site
  
  writeRaster(ras, paste0(wd, "1_Input/UAV/CK/2_Orthophotos_CK/", ras_names[1]), format="GTiff", overwrite=TRUE, datatype='INT4S')

}

