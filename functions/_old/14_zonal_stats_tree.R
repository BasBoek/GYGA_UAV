library(rgeos) 

# VI_mask <- NDREG_mask
# palms_F <- palm
# buffer  <- buf_dist 

zonal_stats_tree <- function(VI_mask, VI_name, palms_F, BRK_SET, palm_name){
  
  RESULT_NR <- 0
  
  # Calculate per palm tree
  for(PALM_ID in 1:length(palms_F)){
    
    palm_crown     <- palms_F[palms_F@data$TRT_ID == palm_name,]

    # Buffer polygons
    img_crop    <- mask(VI_mask, palm_crown)
    
    #GREEN_crop  <- mask(GREEN_mask, palm_buf)

    # Only proceed when values are present
    checkmax    <- mean(maxValue(VI_mask)) # if 0 --> do not make zonal stats
    if(checkmax > 0){

      pixels    <- as.matrix(img_crop)
      npixels   <- length(pixels) - sum(is.na(pixels)) 
      b_mean    <- mean(pixels, na.rm=T)
      b_sd      <- sd(pixels, na.rm=T)
      result    <- as.data.frame(t(as.data.frame(c(
                                                   palm_name,
                                                   VI_name,
                                                   BRK_SET,
                                                   npixels, 
                                                   b_mean, 
                                                   b_sd))))
      
    } else {
      result <- as.data.frame(t(c(palm_name, VI_name, BRK_SET, NA, NA, NA)))
    }
  } 
  names(result) <- c(
                     "palm_id",
                     "VI",
                     "Break_cat", 
                     "pix_nr", 
                     "pix_mean", 
                     "pix_sd")
  rownames(result) <- NULL
  numcols <- c(4:ncol(result))
  result[numcols] <- sapply(result[numcols],as.numeric)
  
  
  return(result)
}  
