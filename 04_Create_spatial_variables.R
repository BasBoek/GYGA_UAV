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
# 0.5 m           1.0 m                 1.5 m
# 1.0 m           1.5 m                 2.0 m
# 1.5 m           2.0 m                 2.5 m
# 2.0 m           2.5 m               
# 2.5 m               


# Function that combines differently aggregated date (shell around the function 'aggregator')
combine_aggs <- function(df, VOI, GOI, functions){
  agg1       <- aggregator(df, VOI, GOI, functions[1])
  agg2       <- aggregator(df, VOI, GOI, functions[2])
  agg        <- cbind(agg1[,2:(ncol(agg1))], agg2[,VOI])
  names(agg) <- c(GOI, paste0(VOI, "__", functions[1]), paste0(VOI, "__", functions[2]))
  return(agg)
}

VIs  <- c("NDVI", "NDRE")
cuts <- c(0.1, 0.3)


for(CUT in 1:length(cuts)){
  for(VI_NR in 1:length(VIs)){
    VI <- VIs[VI_NR]
    
    
    print(paste0("Process ", as.character(cuts[CUT] * 100), ":"))
    
    # Load and preprocess field data
    field_files  <- list.files(paste0(wd, "2_Intermediate"), pattern = "ZonalStats_")
    field_VIs    <- field_files[grepl(VI, field_files)]
    field_df     <- do.call(rbind, lapply(paste0(wd, "2_Intermediate/", field_VIs), read.csv))
    field_df$VI  <- VI # Can be removed after complete run of all scripts
    field_df_cut <- field_df[field_df$Break1 >= (0+cuts[CUT]) & field_df$Break2 <= (1-cuts[CUT]),]
    
    # Load and preprocess palm data
    palm_files        <- list.files(paste0(wd, "2_Intermediate"), pattern = "ZonalStatsPerPalm")
    palm_df           <- do.call(rbind, lapply(paste0(wd, "2_Intermediate/", palm_files), read.csv))
    palm_df           <- palm_df[palm_df$VI == VI,]
    palm_df$VI        <- VI # Can be removed after complete run of all scripts
    palm_df$ID        <- paste0(palm_df[,"Buffer_m"], "_", palm_df[,"Break1"])
    palm_df$barcenter <- palm_df$Break2 - (palm_df$Break2 - palm_df$Break1)/2
    palm_df$barwidth  <- palm_df$Break2 - palm_df$Break1
    
    
    ############
    ### 1) ALL NON-DONUT CALCULATIONS
    ############
    
    # Create weight for every buffer - pixel break combination
    # (inverted normalized coefficient of variation based on the mean mean and mean sd of pixels across the fields)
    FI_pixmean_cut            <- combine_aggs(field_df_cut, "pix_mean", c("Buffer_m", "Break1"), c("mean", "sd"))
    
    buffers <- unique(field_df$Buffer_m)
    
    print(paste0("(1) Calculate weighing based on field shadows and shines for ", VI))
    for(i in 1:length(buffers)){
      BUF <- unique(buffers)[i]
      field_temp     <- field_df[field_df$Buffer_m == BUF,]
      field_temp_cut <- field_df_cut[field_df_cut$Buffer_m == BUF,]
      
      FI_pixmean                <- combine_aggs(field_temp, "pix_mean", c("Buffer_m", "Break1"), c("mean", "sd"))
      FI_pixmean$ID             <- paste0(FI_pixmean[,"Buffer_m"], "_", FI_pixmean[,"Break1"])
      FI_pixmean$pix_cov        <- FI_pixmean$pix_mean__sd / FI_pixmean$pix_mean__mean
      FI_pixmean$weight         <- 1000 - ( FI_pixmean$pix_cov - min(FI_pixmean$pix_cov) ) / max(FI_pixmean$pix_cov - min(FI_pixmean$pix_cov) ) * 1000
      
      FI_pixmean_cut            <- combine_aggs(field_temp_cut, "pix_mean", c("Buffer_m", "Break1"), c("mean", "sd"))
      FI_pixmean_cut$ID         <- paste0(FI_pixmean_cut[,"Buffer_m"], "_", FI_pixmean_cut[,"Break1"])
      FI_pixmean_cut$pix_cov    <- FI_pixmean_cut$pix_mean__sd / FI_pixmean_cut$pix_mean__mean
      FI_pixmean_cut$weight     <- 1000 - ( FI_pixmean_cut$pix_cov - min(FI_pixmean_cut$pix_cov) ) /  max(FI_pixmean_cut$pix_cov - min(FI_pixmean_cut$pix_cov) ) * 1000
      
      if(i == 1){
        FI_pix <- FI_pixmean
        FI_cut <- FI_pixmean_cut
      } else {
        FI_pix <- rbind(FI_pix, FI_pixmean)
        FI_cut <- rbind(FI_cut, FI_pixmean_cut)
      }
      
    }
    # plot(FI_cut$Break1, FI_cut$weight, col=FI_cut$Buffer_m)
    
    agg_uncut <- FI_pix[,c("ID", "weight")]
    agg_cut   <- FI_cut[,c("ID", "weight")]
    
    # Add calculated weighs (based on field statistics) to the palms 
    palms_uncut <- merge(agg_uncut, palm_df, by="ID") # (in case that) palm buffers are not yet been done for 3m and 3.5m, data is omitted here
    palms_cut   <- merge(palm_df, agg_cut, by="ID", all.x=F)
    
    # Combine in 1 df
    palms_uncut$w_type <- "all"
    palms_cut$w_type   <- "sel"
    palms <- rbind(palms_uncut, palms_cut)
    
    
    print(paste("(2) Calculate weighed (based on step 1) and unweighed means for every palm CIRCLE for", VI))
    
    palms$mean_weighed   <- NA
    palms$mean_unweighed <- NA
    i <- 0
    for(TYPE in unique(palms$w_type)){
      for(PALM in unique(palms$palm_id)){
        for(PR_FI in unique(palms$PR_FI)){
          for(BUF in unique(palms$Buffer_m)){
            
            i <- i + 1
            sel <- palms[palms$PR_FI == PR_FI &
                           palms$palm_id == PALM &
                           palms$w_type == TYPE &
                           palms$Buffer_m == BUF,]
            
            sel$mean_weighed   <- (sel$pix_mean * sel$pix_nr * sel$weight / sum(sel$pix_nr[sel$pix_nr > 0]) ) / sum(sel$weight[sel$pix_nr > 0])
            sel$mean_unweighed <- (sel$pix_mean * sel$pix_nr / sum(sel$pix_nr[sel$pix_nr > 0]) )
            
            palms[palms$PR_FI == PR_FI &
                    palms$palm_id == PALM &
                    palms$w_type == TYPE &
                    palms$Buffer_m == BUF, c("mean_weighed", "mean_unweighed")] <- sel[,c("mean_weighed", "mean_unweighed")]
            
          }
        }
      }
    }
    
    # Calculate statistics with weighed and unweighed values 
    palms_w <- aggregator(palms, "mean_weighed",   c("w_type", "PR_FI", "palm_id", "VI", "Buffer_m"), "sum")
    palms_u <- aggregator(palms, "mean_unweighed", c("w_type", "PR_FI", "palm_id", "VI", "Buffer_m"), "sum")
    
    plot(palms_u$mean_unweighed[palms_u$w_type=="all"], palms_w$mean_weighed[palms_w$w_type=="all"], 
         col = rgb(red = 0, green = 0, blue = 0, alpha = 0.15))
    plot(palms_u$mean_unweighed[palms_u$w_type=="sel"], palms_w$mean_weighed[palms_w$w_type=="sel"], 
         col = rgb(red = 0, green = 0, blue = 0, alpha = 0.15))
    
    hist(palms_u$mean_unweighed[palms_u$w_type=="all"])
    hist(palms_u$mean_unweighed[palms_u$w_type=="sel"])
    hist(palms_w$mean_weighed[palms_w$w_type=="all"])
    hist(palms_w$mean_weighed[palms_w$w_type=="sel"])
    
    
    ############
    ### ALL DONUT CALCULATIONS
    ############
    
    donuts <- list(
      c(0.5, 1.0),
      c(0.5, 1.5),
      c(0.5, 2.0),
      c(0.5, 2.5),
      c(1.0, 1.5),
      c(1.0, 2.0),
      c(1.0, 2.5)
    )
    
    palm_buffers <- unique(palms$Buffer_m)
    
    print(paste("(3) Calculate weighed and unweighed means (based on step 1 and 2) for every palm DONUT for", VI))
    
    for(i in 1:length(donuts)){
      
      tiny <- donuts[i][[1]][1]
      big  <- donuts[i][[1]][2]
      
      new  <- palms
      new$pix_mean_donut       <- NA
      new$pix_nr_donut         <- NA
      new$weight_donut         <- NA
      new$mean_donut_weighed   <- NA
      new$mean_donut_unweighed <- NA
      
      for(TYPE in unique(new$w_type)){
        for(PALM in unique(new$palm_id)){
          for(PR_FI in unique(new$PR_FI)){
            
            sel <- new[new$PR_FI == PR_FI &
                         new$palm_id == PALM &
                         new$w_type == TYPE,]
            sel_big  <- sel[sel$Buffer_m == big,]
            sel_tiny <- sel[sel$Buffer_m == tiny,]
            
            pix_nr_donut   <- sel_big$pix_nr - sel_tiny$pix_nr
            pix_mean_donut <- (( sel_big$pix_mean * sel_big$pix_nr ) - (  sel_tiny$pix_mean * sel_tiny$pix_nr)) /  pix_nr_donut
            weight_donut   <- (( sel_big$weight * pi*big^2 ) + (  sel_tiny$weight * pi*tiny^2)) /  ( pi*big^2 + pi*tiny^2 )
            
            mean_donut_weighed   <- (pix_mean_donut * pix_nr_donut * weight_donut / sum(pix_nr_donut[pix_nr_donut > 0])) / sum(weight_donut[weight_donut > 0])
            mean_donut_unweighed <- pix_mean_donut * pix_nr_donut / sum(pix_nr_donut[pix_nr_donut > 0])
            
            pix_nr_donut         <- rep(pix_nr_donut, length(palm_buffers))         # later on we'll remove doubles
            pix_mean_donut       <- rep(pix_mean_donut, length(palm_buffers))       # later on we'll remove doubles
            weight_donut         <- rep(weight_donut, length(palm_buffers))         # later on we'll remove doubles
            mean_donut_weighed   <- rep(mean_donut_weighed, length(palm_buffers))   # later on we'll remove doubles
            mean_donut_unweighed <- rep(mean_donut_unweighed, length(palm_buffers)) # later on we'll remove doubles
            
            new[ new$PR_FI == PR_FI & new$palm_id == PALM & new$w_type == TYPE, "pix_mean_donut"]       <- pix_mean_donut
            new[ new$PR_FI == PR_FI & new$palm_id == PALM & new$w_type == TYPE, "pix_nr_donut"]         <- pix_nr_donut
            new[ new$PR_FI == PR_FI & new$palm_id == PALM & new$w_type == TYPE, "weight_donut"]         <- weight_donut
            new[ new$PR_FI == PR_FI & new$palm_id == PALM & new$w_type == TYPE, "mean_donut_weighed"]   <- mean_donut_weighed
            new[ new$PR_FI == PR_FI & new$palm_id == PALM & new$w_type == TYPE, "mean_donut_unweighed"] <- mean_donut_unweighed
            
          }
        }
      }
      new[,c("ID", "Buffer_m", "Break_val1", "Break_val2", "pix_nr", "pix_mean", "pix_sd", "weight", "mean_weighed", "mean_unweighed")] <- NULL
      new$donut <- paste0(big, "__", tiny)
      new       <- new[!duplicated(new), ] # Remove double rows..
      
      if(i == 1){
        palm_donuts <- new
      } else {
        palm_donuts <- rbind(palm_donuts, new)
      }
    }
    
    # Calculate statistics with weighed and unweighed values of the donuts
    donuts_w <- aggregator(palm_donuts, "mean_donut_weighed",   c("w_type", "PR_FI", "palm_id", "VI", "donut"), "sum")
    donuts_u <- aggregator(palm_donuts, "mean_donut_unweighed", c("w_type", "PR_FI", "palm_id", "VI", "donut"), "sum")
    
    plot(donuts_u$mean_donut_unweighed[donuts_u$w_type=="all"], donuts_w$mean_donut_weighed[donuts_w$w_type=="all"], 
         col = rgb(red = 0, green = 0, blue = 0, alpha = 0.15))
    
    plot(donuts_u$mean_donut_unweighed[donuts_u$w_type=="sel"], donuts_w$mean_donut_weighed[donuts_w$w_type=="sel"], 
         col = rgb(red = 0, green = 0, blue = 0, alpha = 0.15))
    
    
    
    ############
    ### COMBINE ALL INFO
    ############
    
    
    names(palms_u)[names(palms_u)   == 'Buffer_m'] <- 'shape'
    names(palms_w)[names(palms_w)   == 'Buffer_m'] <- 'shape'
    names(donuts_u)[names(donuts_u) == 'donut']    <- 'shape'
    names(donuts_w)[names(donuts_w) == 'donut']    <- 'shape'
    
    palms_u$weighed  <- "no"
    palms_w$weighed  <- "yes"
    donuts_u$weighed <- "no"
    donuts_w$weighed <- "yes"
    
    names(palms_u)[names(palms_u)   == 'mean_unweighed']       <- 'output'
    names(palms_w)[names(palms_w)   == 'mean_weighed']         <- 'output'
    names(donuts_u)[names(donuts_u) == 'mean_donut_unweighed'] <- 'output'
    names(donuts_w)[names(donuts_w) == 'mean_donut_weighed']   <- 'output'
    
    vars_VI <- rbind(palms_u, palms_w, donuts_u, donuts_w)
    
    if(VI_NR == 1){
      vars <- vars_VI
    } else {
      vars <- rbind(vars, vars_VI)
    }
    
    # Transform and write to file
    vars$PR  <- unlist(map(strsplit(vars$PR_FI, split = "_"), 1))
    vars$FI  <- unlist(map(strsplit(vars$PR_FI, split = "_"), 2))
    vars$TRT <- unlist(map(strsplit(vars$palm_id, split = "_"), 1))
    vars$NR  <- paste0("P", unlist(map(strsplit(vars$palm_id, split = "_"), 2)))
    vars$ID  <- paste(vars$PR, vars$FI, vars$TRT, vars$NR, sep="_")
    vars$var <- paste0(vars$VI, "___shape_", vars$shape, "___weight_", vars$weighed)
    df       <- dcast(vars, ID + PR + FI + TRT + NR ~ var, value.var = "output", mean)
    
    write.csv(df, paste0(wd, "3_Output/04_spatial_variables_cut", as.character(cuts[CUT]*100), ".csv"), row.names=F)
    
  }
  
}



# # Cleaning?
# rm(new, agg_uncut, agg_cut, palms_cut, FI_pixmean_cut, FI_pixmean, field_temp, FI_cut, FI_pix,
#    palms_uncut, palm_donuts, field_temp_cut, sel, sel_big, sel_tiny, 
#    vars_VI)










