# Bastiaen Boekelo, April 2021
# Nebraska Project
# Inspect some statistics 

##########################################
#### CHECK WD's & ADAPT IF NECESSARY #####
##########################################


library(raster)
library(rasterpdf)

dev.off(which = grDevices::dev.cur())
rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

fields      <- unique(substr(list.files(paste0(wd, "2_Intermediate/01b_singlebands"), pattern = ".tif$"), 10, 14))


pdf(paste0(wd, "3_Output/01d_Stats_VIs_v2.pdf"))
plot.new()

PR_FI_files <- list.files(paste0(wd, "2_Intermediate"), pattern = ".tif$")
for(PR_FI in fields){

  print(PR_FI)
  
  sel   <- PR_FI_files[grepl( PR_FI, PR_FI_files)]
  MSAVI <- raster(paste0(wd, "2_Intermediate/", sel[1]))
  NDRE  <- raster(paste0(wd, "2_Intermediate/", sel[2]))
  NDREG <- raster(paste0(wd, "2_Intermediate/", sel[3]))
  NDVI  <- raster(paste0(wd, "2_Intermediate/", sel[4]))
  NDVIG <- raster(paste0(wd, "2_Intermediate/", sel[5]))
  WDRVI <- raster(paste0(wd, "2_Intermediate/", sel[7]))
  
  plot(MSAVI)
  title(sel[1])
  plot(NDRE)
  title(sel[2])
  plot(NDREG)
  title(sel[3])
  plot(NDVI)
  title(sel[4])
  plot(NDVIG)
  title(sel[5])
  plot(WDRVI)
  title(sel[7])
  
}

dev.off()


# Amelioration: Make Green first relative to RGB before applying a weight make relative again



# plot(NDREG, NDRE, title(sel[1]))
# plot(NDVIG, NDVI)










