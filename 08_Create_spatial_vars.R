# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Create different spatial variables on field and tree level

rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

source("functions/aggregator.R")
library(maditr) # dcast
library(purrr) # map (during strsplit)

###################################
## 48 variables to create per VI ##
###################################

## 4 Different ways to post-process:
# Weighed
# Weighed     & Cut off first and last 10%
# Not weighed
# Not weighed & Cut off first and last 10%

## 12 Different areas for which we calculate zonal statistics:
# All circles:    Donuts with 0.5 m:    Donuts with 1.0 m:
# 0.5 m           1.0 m                 
# 1.0 m           1.5 m                 
# 1.5 m           2.5 m                 2.5 m
# 2.5 m           3.5 m                 3.5 m
# 3.5 m               


# Function that combines differently aggregated date (shell around the function 'aggregator')
combine_aggs <- function(df, VOI, GOI, functions){
  agg1       <- aggregator(df, VOI, GOI, functions[1])
  agg2       <- aggregator(df, VOI, GOI, functions[2])
  agg        <- cbind(agg1[,2:(ncol(agg1))], agg2[,VOI])
  names(agg) <- c(GOI, paste0(VOI, "__", functions[1]), paste0(VOI, "__", functions[2]))
  return(agg)
}

# Which donuts to create
donuts <- list(
  c(0.5, 1.0),
  c(1.0, 1.5),
  c(1.5, 2.5),
  c(2.5, 3.5)
)

# Load and preprocess palm data
palm_files        <- list.files(paste0(wd, "2_Intermediate"), pattern = "07_ZonalStatsPerPalm")
palms_all         <- do.call(rbind, lapply(paste0(wd, "2_Intermediate/", palm_files), read.csv))
palms_all         <- palms_all[rowSums(is.na(palms_all)) != ncol(palms_all), ]
palms_all$ID      <- paste0(palms_all$PR_FI, "_", substr(palms_all$palm_id ,1, 4), "P", substr(palms_all$palm_id,5,6))

#palms_all$ID      <- paste0(palms_all$PR_FI, "_P", palms_all$palm_id)
palms_all$shape   <- palms_all$Buffer_m
palms_all$PR_FI   <- NULL
palms_all$palm_id <- NULL

VIs      <- unique(palms_all$VI)

brk_cats <- unique(palms_all$Break_cat)
#buffers  <- unique(palms_all$Buffer_m)
IDs      <- unique(palms_all$ID)
nr_rows  <- length(brk_cats) * length(IDs)

CAT   <- 1
BUF   <- 1
ID_NR <- 1
D_NR  <- 1



############
### ALL DONUT CALCULATIONS
############


donut_vars <- function(VI){
  i <- 0
  for(CAT in 1:length(brk_cats)){ # For all break categories
    brk_cat <- brk_cats[CAT]
    print(paste(VI, brk_cat))
    for(ID_NR in 1:length(IDs)){ # For every palm tree
      ID <- IDs[ID_NR]
      
      # Select all buffers of a palm
      palm <- palms_all[palms_all$VI == VI & palms_all$Break_cat == brk_cat & palms_all$ID == ID, ]
      
      for(D_NR in 1:length(donuts)){
        i <- i + 1
        tiny  <- donuts[D_NR][[1]][1]
        big   <- donuts[D_NR][[1]][2]
        shape <- paste0(big, "-", tiny)
        
        sel_big        <- palm[palm$Buffer_m == big,]
        sel_tiny       <- palm[palm$Buffer_m == tiny,]
        
        if(nrow(sel_big) == 1 & nrow(sel_tiny) == 1){
          
          pix_nr_donut   <- sel_big$pix_nr - sel_tiny$pix_nr
          pix_mean_donut <- (( sel_big$pix_mean * sel_big$pix_nr ) - (  sel_tiny$pix_mean * sel_tiny$pix_nr)) /  pix_nr_donut
          
          if(i == 1){
            new             <- palm
            new[1:nr_rows,] <- NA
            new[i,] <- c(VI, NA, brk_cat, pix_nr_donut, pix_mean_donut, NA, ID, shape)
            
          } else {
            new[i,] <- c(VI, NA, brk_cat, pix_nr_donut, pix_mean_donut, NA, ID, shape)
          }
          
        } else {
          new[i,] <- c(VI, NA, brk_cat, NA, NA, NA, ID, shape) # No pix_mean & pix_nr when donut could not be made (by absence of at least one buffer)
        }
      }
    }
  }
  return(new)
}



############
### COMBINE ALL INFO
############

for(VI_NR in 1:length(VIs)){
  
  VI  <- VIs[VI_NR]
  
  new <- donut_vars(VI)
  
  if(VI_NR == 1){
    palms <- rbind(palms_all, new)
  } else {
    palms <- rbind(palms, new)
  }
  
}
palms$Buffer_m  <- NULL

palms$PR    <- unlist(map(strsplit(palms$ID, split = "_"), 1))
palms$FI    <- unlist(map(strsplit(palms$ID, split = "_"), 2))
palms$TRT   <- unlist(map(strsplit(palms$ID, split = "_"), 3))

palms$var   <- paste0(palms$VI, "_cat", palms$Break_cat)

shapes      <- unique(palms$shape)
for(SHAPE in shapes){
  palms_shp  <- palms[palms$shape == SHAPE,]
  df         <- dcast(palms_shp, ID + PR + FI + NR + TRT ~ var, value.var = "pix_mean")
  write.csv(df, paste0(wd, "3_Output/08_spatial_variables_shape_", SHAPE, ".csv"), row.names=F)
}














