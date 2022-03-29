library(rgeos) 


# brks_all   <- c(0.50, 0.65, 0.80)

# Create percentile breaks
brks_all   <- c(0.00, 0.01, 0.02, 
                0.05, 0.10, 0.15,
                0.20, 0.35, 0.50,
                0.65, 0.80, 0.85,
                0.90, 0.95, 0.98,
                0.99, 1.00)

zonal_stats_tree <- function(PR_FI, VI_mask, palms_F, buffer, brks, palm_name){
  
  RESULT_NR <- 0
  
  for(BRK in 1:(length(brks_all)-1)){
    
    #Calculate per palm tree
    for(PALM_ID in 1:length(palms_F)){
      
      palm_pt   <- palms_F[palms_F@data$TRT_ID == palm_name,]
      
      # Buffer polygons
      palm_buf   <- buffer(palm_pt, width=buffer)   # Buffer point
      
      img_crop    <- mask(VI_mask, palm_buf)
      
      checkmax    <- mean(maxValue(VI_mask)) # if 0 --> do not make zonal stats
      
      if(checkmax > 0){
        i <- i + 1
        
        img_crop  <- mask(VI_mask, palm_buf)
        pixels    <- as.matrix(img_crop)
        npixels   <- length(pixels) - sum(is.na(pixels)) 
        b_mean    <- mean(pixels, na.rm=T)
        b_sd      <- sd(pixels, na.rm=T)
        newrow    <- t(as.data.frame(c(PR_FI,
                                       palm_name,
                                       substr(names(VI_mask),5,8),
                                       buffer, 
                                       brks[1],
                                       brks[2],
                                       brk_vals[1],
                                       brk_vals[2], 
                                       npixels, 
                                       b_mean, 
                                       b_sd)))
        if(BRK == 1){
          result_brk  <- newrow
        } else {
          result_brk  <- rbind(result_brk, newrow)
        }
      }
      
      if(RESULT_NR == 0){
        result     <- as.data.frame(result_brk)
        RESULT_NR  <- RESULT_NR + 1
        
      } else {
        result     <- rbind(result, result_brk)
      }
    }
    
    names(result) <- c("PR_FI", 
                       "palm_id",
                       "VI",
                       "Buffer_m", 
                       "Break1", 
                       "Break2", 
                       "Break_val1", 
                       "Break_val2", 
                       "pix_nr", 
                       "pix_mean", 
                       "pix_sd")
    rownames(result) <- NULL
    
    numcols <- c(4:11)
    result[numcols] <- sapply(result[numcols],as.numeric)
    
    return(result)
  }
}

