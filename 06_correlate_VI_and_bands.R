# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Inspect how VI's correspond with orthophotos

rm(list=ls())

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

library(raster)

# Load data
file_VI    <- paste0(wd, "2_Intermediate/01_NDVI_RI_F4_Mask.tif")
file_ortho <- paste0(wd, "2_Intermediate/01_Ortho_RI_F4_Mask.tif")

ras_NDVI  <- raster(file_VI)
ras_ortho <- stack(file_ortho)
ras       <- ras_ortho

ras_ortho[[1]] <- ( ras[[1]] - cellStats(ras[[1]], "mean", na.rm=T) ) / cellStats(ras[[1]], "sd", na.rm=T)
ras_ortho[[2]] <- ( ras[[2]] - cellStats(ras[[2]], "mean", na.rm=T) ) / cellStats(ras[[2]], "sd", na.rm=T)
ras_ortho[[3]] <- ( ras[[3]] - cellStats(ras[[3]], "mean", na.rm=T) ) / cellStats(ras[[3]], "sd", na.rm=T)
ras_ortho[[4]] <- ( ras[[4]] - cellStats(ras[[4]], "mean", na.rm=T) ) / cellStats(ras[[4]], "sd", na.rm=T)
ras_ortho[[5]] <- ( ras[[5]] - cellStats(ras[[5]], "mean", na.rm=T) ) / cellStats(ras[[5]], "sd", na.rm=T)

plot(ras_NDVI, ras_ortho[[1]], maxpixels=200000)#, ylim = c(0, 20000))
plot(ras_NDVI, ras_ortho[[2]], maxpixels=200000)#, ylim = c(0, 20000))
plot(ras_NDVI, ras_ortho[[3]], maxpixels=200000)#, ylim = c(0, 20000))
plot(ras_NDVI, ras_ortho[[4]], maxpixels=200000)#, ylim = c(0, 20000))
plot(ras_NDVI, ras_ortho[[5]], maxpixels=200000)#, ylim = c(0, 30000))



