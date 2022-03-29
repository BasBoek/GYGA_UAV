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
                 NRED     = numeric(),
                 NGRN     = numeric(),
                 NNIR     = numeric())
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
  
  # Load palms belonging to rasters
  palms_PR_FI <- all_palms[all_palms$PR_FI == PR_FI,]
  
  if(length(palms_PR_FI) != 0){
    
    palms_PR_FI <- spTransform(palms_PR_FI, crs(NIR) )
    
    # Calculate for every palm
    for(PALM_ID in palms_PR_FI$TRT_ID){
      i <- i + 1
      print(paste(PR_FI, PALM_ID, i))
      
      palm          <- palms_PR_FI[palms_PR_FI$TRT_ID == PALM_ID,]
      palmbuf       <- buffer(palm, 4)
      
      NRED          <- round(RED / (GREEN + NIR + RED) * 1000)
      NRED_mask     <- crop(NRED, palmbuf)
      NRED          <- cellStats(NRED_mask, "mean",  na.rm=T) 
      
      NGRN          <- round(GREEN / (GREEN + NIR + RED) * 1000)
      NGRN_mask     <- crop(NGRN, palmbuf)
      NGRN          <- cellStats(NGRN_mask, "mean",  na.rm=T) 
      
      NNIR          <- round(NIR / (GREEN + NIR + RED) * 1000)
      NNIR_mask     <- crop(NNIR, palmbuf)
      NNIR          <- cellStats(NNIR_mask, "mean",  na.rm=T) 

      
      df[i,] <- c(PR_FI, PALM_ID, NRED, NGRN, NNIR)
    }
  }
}
write.csv(df, "../Data/3_Output/11_band_norms.csv", row.names = F)



