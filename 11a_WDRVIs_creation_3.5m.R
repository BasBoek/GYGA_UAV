# Bastiaen Boekelo, April 2021
# Nebraska Project
# Calculate different VIs

library(rgdal)
library(raster)

rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

PR_FI_files <- list.files(paste0(wd, "2_Intermediate/01b_singlebands"), pattern = ".tif$")
fields      <- unique(substr(list.files(paste0(wd, "2_Intermediate/01b_singlebands"), pattern = ".tif$"), 10, 14))

all_palms <- readOGR(paste0(wd, "2_Intermediate"), "03_palm_centerpoints")

filter <- matrix(c(1,2,2,2,1,
                   2,3,3,3,2,
                   2,3,3,3,2,
                   2,3,3,3,2,
                   1,2,2,2,1), nrow=5, byrow=TRUE) 
filter <- matrix(c(1,1,1,
                   1,0,1,
                   1,1,1), nrow=3, byrow=TRUE) 


df <- data.frame(PR_FI    = character(), 
                 palm_id  = character(),
                 WDRVI    = numeric(),
                 GWWDRVI  = numeric(),
                 NWWDRVI  = numeric(),
                 NIR_FI   = numeric(),
                 RED_FI   = numeric(),
                 GREEN_PA = numeric())
df[1:length(all_palms),] <- NA
i <- 0
FI_nr <- 1
PR_FI <- "CK_F1"
PALM_ID <- "REF_15"
for(FI_nr in 1:length(fields)){
  
  # Define field 
  PR_FI       <- fields[FI_nr]
  PR_FI_sel   <- PR_FI_files[grepl(PR_FI, PR_FI_files)]
  
  # Load rasters
  GREEN     <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[2])) / 32768 * 100
  NIR       <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[5])) / 32768 * 100
  RED       <- raster(paste0(wd, "2_Intermediate/01b_singlebands/", PR_FI_sel[3])) / 32768 * 100
  GREEN_AVG <- cellStats(GREEN,  "mean", na.rm=T) 
  NIR_AVG   <- cellStats(NIR,    "mean", na.rm=T) 
  RED_AVG   <- cellStats(RED,    "mean", na.rm=T) 
  
  # Load palms belonging to rasters
  palms_PR_FI <- all_palms[all_palms$PR_FI == PR_FI,]
  
  if(length(palms_PR_FI) != 0){
    
    palms_PR_FI <- spTransform(palms_PR_FI, crs(NIR) )
    
    # Calculate for every palm
    for(PALM_ID in palms_PR_FI$TRT_ID){
      i <- i + 1
      print(paste(PR_FI, PALM_ID, i))
      
      palm          <- palms_PR_FI[palms_PR_FI$TRT_ID == PALM_ID,]
      palmbuf       <- buffer(palm, 3.5)
      
     #WDRVI         <- round((0.1 * NIR - RED) / (0.1 * NIR + RED) * 1000)
      WDRVI         <- round((0.05 * NIR - RED) / (0.05 * NIR + RED) * 1000)
      WDRVI_mask    <- crop(WDRVI, palmbuf)
      
      # Green weight
      GREEN_mask    <- crop(GREEN, palmbuf)
      GRN_std       <- cellStats(GREEN_mask, "sd", na.rm=T)
      GRN_avg       <- cellStats(GREEN_mask, "mean", na.rm=T)
      #GREEN_weight  <- ( (GRN_std/1.5) / ( GRN_std/1.5 + (GREEN_mask - GRN_avg)^2 )^1)
      GREEN_weight  <- GRN_std /  ( GRN_std + (GREEN_mask - GRN_avg)^2 )
      #GREEN_weight  <- sqrt( (1) / ( 1 + (GREEN_mask - (GRN_avg - GRN_std/1)) )^2)
      #GREEN_weight  <- ( GRN_std/1.5 / ( (GREEN_mask - GRN_avg - GRN_std/2)^2 )^1)
      #GREEN_weight    <- focal(GREEN_weight, filter, fun="sum", na.rm=T)
      weight_GRN    <- GREEN_weight / cellStats(GREEN_weight, "sum", na.rm=T)
      #plot(GREEN_weight)

      # NIR weight
      NIR_mask      <- crop(NIR, palmbuf)
      NIR_std       <- cellStats(NIR_mask, "sd", na.rm=T)
      NIR_avg       <- cellStats(NIR_mask, "mean", na.rm=T)
      #NIR_weight    <- ( (NIR_std/1.5) / ( NIR_std/1.5 + (NIR_mask - NIR_std)^2 )^1)
      NIR_weight    <- ( (NIR_std) / ( NIR_std + (NIR_mask - (NIR_avg - NIR_std*0.75))^2 )^1)
      #NIR_weight    <- focal(NIR_weight, filter, fun="sum", na.rm=T)
      weight_NIR    <- NIR_weight / cellStats(NIR_weight, "sum", na.rm=T)
      #plot(weight_NIR)

      # Variables
      WDRVI_avg     <- cellStats(WDRVI_mask,    "mean", na.rm=T)

      WDRVI_Gw      <- WDRVI_mask * weight_GRN
      GWWDRVI       <- cellStats(WDRVI_Gw, "sum",  na.rm=T)

      WDRVI_Nw      <- WDRVI_mask * weight_NIR
      NWWDRVI       <- cellStats(WDRVI_Nw, "sum",  na.rm=T)



      df[i,] <- c(PR_FI, PALM_ID, WDRVI_avg, GWWDRVI, NWWDRVI, NIR_AVG, RED_AVG, GRN_avg)
    }
  }
}
write.csv(df, "../Data/3_Output/11_WDRVIs_a005__3.5m.csv", row.names = F)

