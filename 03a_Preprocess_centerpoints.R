# Bastiaen Boekelo, March 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Read palm points and write in consistent way to Data/02_Intermediate folder, for proper naming convention


##########################################
#### CHECK WD's & ADAPT IF NECESSARY #####
##########################################


rm(list=ls())  # Clean script <- clean mind

library(rgdal)
source("functions/global_variables.R")
source("functions/add_zero_to_single_characters_in_vector.R")

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd   <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

# List shapefiles to use
shps <- list.files(paste0(wd, "2_Intermediate/"), pattern="shp$")
shps <- shps[grepl("01_Points_", shps)] # Make sure only right input points
shps <- shps[!grepl("CK_F4", shps)] # Remove Field 4 (not present in ground measurements and no TREAT_NO combination) 

# Preprocess per file
for(SHP_NR in 1:length(shps)){
  
  PR_FI_name <- substr(shps[SHP_NR], 11, 15)
  print(PR_FI_name)
  
  shp  <- readOGR(paste0(wd, "2_Intermediate"), substr(shps[SHP_NR], 1 , nchar(shps[SHP_NR]) - 4))
  
  # Create consistent columns, try to catch all exceptions
  if("NO_TS"      %in% colnames(shp@data)){      shp@data <- add_zero(shp@data, "NO_TS") }       
  if("N0_TS"      %in% colnames(shp@data)){      shp@data <- add_zero(shp@data, "N0_TS")  }     
  if("SPL_POINTS" %in% colnames(shp@data)){      shp@data <- add_zero(shp@data, "SPL_POINTS") } 
  if("TRIAL"      %in% colnames(shp@data)){     shp@data[, "TREATMENT"] <- shp@data[, "TRIAL"] }
  
  
  shp@data[,"TRT_ID"] <- paste(sep="_", shp@data[,"TREATMENT"], shp@data[,"NR_ID"])
  shp@data[,"PR_FI"] <- PR_FI_name
  shp@data <- shp@data[,c("TRT_ID", "PR_FI")]
  
  if(SHP_NR == 1){
    shp_all <- shp
  } else {
    shp_all <- rbind(shp_all, shp)
  }
}

# Write all in one shapefile
writeOGR(shp_all, paste0(wd, "2_Intermediate"), "03_palm_centerpoints", driver="ESRI Shapefile", overwrite_layer=TRUE)





