# Bastiaen Boekelo, May 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Inspect correlations between created spatial variables and ground measurement

library(rgdal)
library(raster)

rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd  <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

provinces <- c("CK", "JB", "RI", "SS")
fields    <- paste0("F", 1:8)

combis    <- expand.grid(provinces, fields)
combis    <- paste0(combis$Var1, "_", combis$Var2)
combis    <- sort(combis)
combis    <- combis[combis != "CK_F7"]
combis    <- combis[combis != "CK_F8"]

# Existing combis
temp   <- list.files(paste0(wd, "2_Intermediate/"), pattern = "13_Treatment_weight")
temp   <- substr(temp, 21,25)
combis <- setdiff(combis, temp)


FI_tifs   <- list.files(paste0(wd, "2_Intermediate/01a3_FieldMask/"), pattern = ".tif$")
FI_shps   <- list.files(paste0(wd, "2_Intermediate/01a3_FieldMask/"), pattern = ".shp$")

df <- data.frame(PR_FI    = character(), # Field name
                 treatment= character(), # REF  / BMP
                 metric   = character(), # mean / sd
                 blue     = numeric(),
                 green    = numeric(),
                 red      = numeric(),
                 RE       = numeric(),
                 NIR      = numeric(),
                 NDVI     = numeric(),
                 NDRE     = numeric(),
                 WDRVI    = numeric(),
                 MSAVI    = numeric(),
                 NRED     = numeric(),
                 NGRN     = numeric(),
                 NNIR     = numeric(),
                 GWNDVI   = numeric(),
                 GWNDRE   = numeric(),
                 GWWDRVI  = numeric(),
                 GWMSAVI  = numeric(),
                 GWNRED   = numeric(),
                 GWNGRN   = numeric(),
                 GWNNIR   = numeric())





treats <- c("BMP", "REF")
i <- 1
for(i in 1:length(combis)){
  
  PR_FI     <- combis[i]
  tif_name  <- FI_tifs[grepl(PR_FI, FI_tifs)]
  shp_name  <- FI_shps[grepl(PR_FI, FI_shps)]
  
  if(length(tif_name) > 0 & length(shp_name) > 0){
    
    print(PR_FI)

    tif       <- stack(paste0(wd, "2_Intermediate/01a3_FieldMask/" , tif_name))
    shp       <- readOGR(paste0(wd, "2_Intermediate/01a3_FieldMask"), substr(shp_name, 1, nchar(shp_name)-4))
    shp       <- spTransform(shp, crs(tif))
    
    BLUE      <- tif[[1]]
    GREEN     <- tif[[2]]
    RED       <- tif[[3]]
    RE        <- tif[[4]]
    NIR       <- tif[[5]]
    
    # VIs
    NDVI      <- round((NIR - RED) / (NIR + RED) * 1000)
    NDRE      <- round((RE  - RED) / (RE  + RED) * 1000)
    WDRVI     <- round((0.1 * NIR - RED) / (0.1 * NIR + RED) * 1000)
    MSAVI     <- round(( (NIR - RED) * 1.5 ) / ( NIR + RED + 0.5 ) * 1000) # using L = 0.5 as compromise (https://wiki.landscapetoolbox.org/doku.php/remote_sensing_methods:modified_soil-adjusted_vegetation_index)
    NRED      <- round(RED / (GREEN + NIR + RED) * 1000)
    NGRN      <- round(GREEN / (GREEN + NIR + RED) * 1000)
    NNIR      <- round(NIR / (GREEN + NIR + RED) * 1000)
    
    # Green Weight
    GRN_WGHT <- GREEN
    GRN_WGHT[GRN_WGHT > quantile(GRN_WGHT, 0.95)] <- 1
    GRN_WGHT <- sqrt(GRN_WGHT)
    
    # Write rasters
    writeRaster(NDVI,     paste0('../Data/2_intermediate/13_Treatment_NDVI_', PR_FI, ".tif") ,  "GTiff", overwrite=T)
    writeRaster(NDRE,     paste0('../Data/2_intermediate/13_Treatment_NDRE_', PR_FI, ".tif"),   "GTiff", overwrite=T)
    writeRaster(WDRVI,    paste0('../Data/2_intermediate/13_Treatment_WDRVI_', PR_FI, ".tif"),  "GTiff", overwrite=T)
    writeRaster(MSAVI,    paste0('../Data/2_intermediate/13_Treatment_MSAVI_', PR_FI, ".tif"),  "GTiff", overwrite=T)
    writeRaster(NRED,     paste0('../Data/2_intermediate/13_Treatment_NRED_', PR_FI, ".tif"),   "GTiff", overwrite=T)
    writeRaster(NGRN,     paste0('../Data/2_intermediate/13_Treatment_NGRN_', PR_FI, ".tif"),   "GTiff", overwrite=T)
    writeRaster(NNIR,     paste0('../Data/2_intermediate/13_Treatment_NNIR_', PR_FI, ".tif"),   "GTiff", overwrite=T)
    writeRaster(GRN_WGHT, paste0('../Data/2_intermediate/13_Treatment_weight_', PR_FI, ".tif"), "GTiff", overwrite=T)
    
    for(TRT in treats){
      
      print(TRT)
      
      # Clip  
      shp_trt       <- shp[shp$TREATMENT == TRT,]
      
      blue_avg  <- cellStats(BLUE,  "mean")
      green_avg <- cellStats(GREEN, "mean")
      red_avg   <- cellStats(RED,   "mean")
      RE_avg    <- cellStats(RE,    "mean")
      NIR_avg   <- cellStats(NIR,   "mean")
      
      blue_std  <- cellStats(BLUE,  "sd")
      green_std <- cellStats(GREEN, "sd")
      red_std   <- cellStats(RED,   "sd")
      RE_std    <- cellStats(RE,    "sd")
      NIR_std   <- cellStats(NIR,   "sd")
      
      
      # Calculate variables
      NDVI_mask     <- mask(NDVI, shp_trt)
      NDRE_mask     <- mask(NDRE, shp_trt)
      WDRVI_mask    <- mask(WDRVI, shp_trt)
      MSAVI_mask    <- mask(MSAVI, shp_trt)
      NRED_mask     <- mask(NRED, shp_trt)
      NGRN_mask     <- mask(NGRN, shp_trt)
      NNIR_mask     <- mask(NNIR, shp_trt)

      GRN_WGHT_mask <- mask(GRN_WGHT, shp_trt)
      GRN_WGHT_mask <- GRN_WGHT_mask / cellStats(GRN_WGHT_mask,  "sum")

      # (1) Derive mean values 
      NDVI_avg   <- cellStats(NDVI_mask,    "mean")
      NDRE_avg   <- cellStats(NDRE_mask,    "mean")
      WDRVI_avg  <- cellStats(WDRVI_mask,   "mean")
      MSAVI_avg  <- cellStats(MSAVI_mask,   "mean")
      NRED_avg   <- cellStats(NRED_mask,    "mean")
      NGRN_avg   <- cellStats(NGRN_mask,    "mean")
      NNIR_avg   <- cellStats(NNIR_mask,    "mean")
      
      NDVI_GW    <- cellStats(NDVI_mask  * GRN_WGHT_mask,  "sum")
      NDRE_GW    <- cellStats(NDRE_mask  * GRN_WGHT_mask,  "sum")
      WDRVI_GW   <- cellStats(WDRVI_mask * GRN_WGHT_mask,  "sum")
      MSAVI_GW   <- cellStats(MSAVI_mask * GRN_WGHT_mask,  "sum")
      NRED_GW    <- cellStats(NRED_mask  * GRN_WGHT_mask,  "sum")
      NGRN_GW    <- cellStats(NGRN_mask  * GRN_WGHT_mask,  "sum")
      NNIR_GW    <- cellStats(NNIR_mask  * GRN_WGHT_mask,  "sum")
      
      row1 <- c(PR_FI, "mean", TRT, blue_avg, green_avg, red_avg, RE_avg, NIR_avg,
                NDVI_avg,  NDRE_avg , WDRVI_avg , MSAVI_avg , NRED_avg , NGRN_avg , NNIR_avg , 
                NDVI_GW ,  NDRE_GW ,  WDRVI_GW ,  MSAVI_GW ,  NRED_GW ,  NGRN_GW ,  NNIR_GW )
      
      # (2) Derive sd values 
      NDVI_sd      <- cellStats(NDVI_mask,    "sd")
      NDRE_sd      <- cellStats(NDRE_mask,    "sd")
      WDRVI_sd     <- cellStats(WDRVI_mask,   "sd")
      MSAVI_sd     <- cellStats(MSAVI_mask,   "sd")
      NRED_sd      <- cellStats(NRED_mask,    "sd")
      NGRN_sd      <- cellStats(NGRN_mask,    "sd")
      NNIR_sd      <- cellStats(NNIR_mask,    "sd")
      
      NDVI_GW_sd   <- cellStats(NDVI_mask * GRN_WGHT_mask,  "sd")
      NDRE_GW_sd   <- cellStats(NDRE_mask * GRN_WGHT_mask,  "sd")
      WDRVI_GW_sd  <- cellStats(WDRVI_mask * GRN_WGHT_mask, "sd")
      MSAVI_GW_sd  <- cellStats(MSAVI_mask * GRN_WGHT_mask, "sd")
      NRED_GW_sd   <- cellStats(NRED_mask * GRN_WGHT_mask,  "sd")
      NGRN_GW_sd   <- cellStats(NGRN_mask * GRN_WGHT_mask,  "sd")
      NNIR_GW_sd   <- cellStats(NNIR_mask * GRN_WGHT_mask,  "sd")
      
      row2 <- c(PR_FI, "sd", TRT, blue_std, green_std, red_std, RE_std, NIR_std,
                NDVI_sd,      NDRE_sd ,     WDRVI_sd ,     MSAVI_sd ,     NRED_sd ,     NGRN_sd ,     NNIR_sd , 
                NDVI_GW_sd ,  NDRE_GW_sd ,  WDRVI_GW_sd ,  MSAVI_GW_sd ,  NRED_GW_sd ,  NGRN_GW_sd ,  NNIR_GW_sd )
      
      # Combine rows and write in table
      rows <- as.data.frame(t(cbind(row1, row2)))
      names(rows) <- names(df)
      write.csv(rows, paste0("../Data/0_Temp/13_Field_Treatment_stats_", PR_FI, "_", TRT, ".csv"), row.names=F)

      df   <- rbind(df, rows)
      
    }
  }
}
# write.csv(df, paste0("../Data/2_Intermediate/13_Field_Treatment_stats.csv"), row.names=F)
# 
# field_stats <- list.files(paste0(wd, "0_Temp"), pattern="13_Field_Treatment_stats_")
# df          <- do.call(rbind, lapply(paste0(wd, "0_Temp/", field_stats), read.csv))








