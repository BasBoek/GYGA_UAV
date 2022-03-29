library(rgeos) 

# Create percentile breaks
brks_all   <- c(0.00, 0.01, 0.02, 
               0.05, 0.10, 0.15,
               0.20, 0.35, 0.50,
               0.65, 0.80, 0.85,
               0.90, 0.95, 0.98,
               0.99, 1.00)

# brks_all   <- c(0.00, 0.01)
# BRK <- 1
# VI <- NDVI
# 
# plot(NDRE)
# plot(NDVI)
# plot(mask_ras)

zonal_stats_field <- function(PR_FI, VI, avg_pos, points, buffer, brks){
  
  RESULT_NR <- 0
  

  for(BRK in 1:(length(brks_all)-1)){
    i <- 0

    mask_ras <- round( avg_pos  / cellStats(avg_pos,  "max", na.rm=T) * 10000)

    brk1 <- brks_all[BRK]
    brk2 <- brks_all[BRK+1]
    print(brk2)

    brks      <- c(brk1, brk2)
    brk_vals  <- as.vector(quantile(mask_ras, brks))

    mask_ras[mask_ras <- mask_ras < brk_vals[1] | mask_ras >= brk_vals[2] ] <- NA
    mask_ras[mask_ras > 0] <- 1
    mask_ras@file@datanotation    <- NDVI@file@datanotation

    VI_mask   <- mask(VI, mask_ras)

    # Buffer polygons
    polygons  <- buffer(points, width=buffer, dissolve=F)   # Buffer points

    checkmax  <- mean(maxValue(VI)) # if 0 --> do not make zonal stats

    if(checkmax > 0){
      i <- i + 1

      img_crop  <- mask(VI_mask, polygons)
      pixels    <- as.matrix(img_crop)
      npixels   <- length(pixels) - sum(is.na(pixels)) 
      b_mean    <- mean(pixels, na.rm=T)
      b_sd      <- sd(pixels, na.rm=T)
      newrow    <- t(as.data.frame(c(PR_FI, 
                                      substr(names(VI),5,8),
                                      buffer, 
                                      brks[1],
                                      brks[2],
                                      brk_vals[1],
                                      brk_vals[2], 
                                      npixels, 
                                      b_mean, 
                                      b_sd)))
      if(i == 1){
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
  
  numcols <- c(3:10)
  result[numcols] <- sapply(result[numcols],as.numeric)
  
  return(result)
}


