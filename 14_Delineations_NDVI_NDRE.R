# Bastiaen Boekelo, June 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Extract info NDVI and NDRE

library(rgdal)
library(raster)

rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd  <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

CK_F4_shp <- readOGR("../Data/1_Input/UAV/CK/6_Crown_Delineations_CK", "Palm_Crown_Micasense_CK_F4")
JB_F7_shp <- readOGR("../Data/1_Input/UAV/JB/6_Crown_Delineations_JB", "Palm_Crown_Projection_Micasense_UNL_JB_F7")
RI_F4_shp <- readOGR("../Data/1_Input/UAV/RI/6_Crown_Delineations_RI", "Palm_Crown_Projection_Micasense_UNL_RI_F4")
SS_F5_shp <- readOGR("../Data/1_Input/UAV/SS/6_Crown_Delineations_SS", "Palm_Crown_Projection_Micasense_UNL_SS_F5")

CK_F4_trt <- readOGR("../Data/1_Input/UAV/CK/5_Treatment_Boundary_CK", "Palm_Crown_Micasense_CK_F4")
JB_F7_trt <- readOGR("../Data/1_Input/UAV/JB/5_Treatment_Boundary_JB", "Treatment_Boundary_JB_F7")
RI_F4_trt <- readOGR("../Data/1_Input/UAV/RI/5_Treatment_Boundary_RI", "Palm_Crown_Projection_Micasense_UNL_RI_F4")
SS_F5_trt <- readOGR("../Data/1_Input/UAV/SS/5_Treatment_Boundary_SS", "Palm_Crown_Projection_Micasense_UNL_SS_F5")


ras   <- raster("../Data/2_Intermediate/13_Treatment_NNIR_JB_F7.tif")
plot(ras)

test  <- mask(ras, JB_F7_shp)
plot(test)

write.csv(CK_F4@data, "../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_CK_F4.csv", row.names=F)
write.csv(JB_F7@data, "../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_JB_F7.csv", row.names=F)
write.csv(RI_F4@data, "../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_RI_F4.csv", row.names=F)
write.csv(SS_F5@data, "../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_SS_F5.csv", row.names=F)


