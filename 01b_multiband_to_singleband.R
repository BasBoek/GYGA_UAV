# Bastiaen Boekelo, April 2021
# Nebraska Project
# Transform the calculated masked orthomosaics to singleband images

library(raster)

rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

ras_files  <- list.files(paste0(wd, "2_Intermediate"), pattern = ".tif$")
ras_files  <- ras_files[grepl("01_Ortho_", ras_files)]

for(FILENR in 1:length(ras_files)){
  ras_name <- ras_files[FILENR]
  ras <- stack(paste0(wd, "2_Intermediate/", ras_name))
  for(i in 1:5){
    writeRaster(ras[[i]], paste0('../Data/2_Intermediate/01b_singlebands/', substr(ras_name, 1, nchar(ras_name)-4), "_b", i, ".tif"), "GTiff", overwrite=T)
  }
}




