# Bastiaen Boekelo, March 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Visualize / inspect how shadowiness and shininess influence the VI response.


##########################################
#### CHECK WD's & ADAPT IF NECESSARY #####
##########################################


rm(list=ls())  # Clean script <- clean mind

# Set Script and data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

library(ggplot2)
source("functions/aggregator.R")
source("functions/make_gradient.R")
source("functions/global_variables.R")

pdf(paste0(wd, "3_Output/02_Shadowiness_Shininess_Fields.pdf"), width=15, height=10)

for(VI in VIs){
  for(buf_dist in buffers){
    print(paste(VI,buf_dist))
    # Load and data
    
    for(COMBI in 1:length(combis)){
      filenames  <- list.files(paste0(wd, "2_Intermediate/"), pattern = "02_ZonalStats_")
      filenames  <- grep(VI, filenames, value = T)
      filenames  <- grep(paste0(buf_dist,"m"), filenames, value = T)
      zonals     <- do.call(rbind, lapply(paste0(wd, "2_Intermediate/", filenames), read.csv))
    }
    
    zonals$weight       <- zonals$pix_nr / max(zonals$pix_nr)
    zonals$pix_cov      <- zonals$pix_sd / zonals$pix_mean
    zonals$barcenter    <- zonals$Break2 - (zonals$Break2 - zonals$Break1)/2
    zonals$barwidth     <- zonals$Break2 - zonals$Break1
    
    # Visualize 
    g1 <- ggplot(zonals) +
      geom_line( aes(x=barcenter, y=pix_mean), size=0.8, colour="#001100") +
      geom_point( aes(x=barcenter, y=pix_mean), size=1.6, colour="#001100") +
      geom_errorbar( aes(x=barcenter, ymin=pix_mean-pix_sd, ymax=pix_mean+pix_sd), width = zonals$barwidth-zonals$barwidth/2, size=0.6, colour="#001100", alpha=0.9) +
      theme_bw() + xlim(0, 1) +
      facet_wrap(~PR_FI) +
      ylab(VI) + xlab("Shadowiness - Shininess") +
      ggtitle(paste0(VI, " - ", buf_dist, " m"))
    plot(g1)
    
    g2 <- ggplot(zonals) +
      geom_line( aes(x=barcenter, y=pix_mean, group=PR_FI,  col=PR_FI), size=0.8, alpha=0.5) +
      geom_point( aes(x=barcenter, y=pix_mean, group=PR_FI,  col=PR_FI), size=1.6) +
      theme_bw() + xlim(0, 1) +
      ylab(paste(VI, "(per field)")) + xlab("Shadowiness - Shininess") +
      ggtitle(paste0(VI, " - ", buf_dist, " m"))
    plot(g2)
    
    agg_avg      <- aggregator(zonals, "pix_mean", "barcenter", "mean")
    agg_sd       <- aggregator(zonals, "pix_mean", "barcenter", "sd")
    agg          <- cbind(agg_avg[,2:3], agg_sd$pix_mean)
    names(agg)   <- c("barcenter", "pix_mean", "pix_sd")
    agg$pix_cov  <- agg$pix_sd / agg$pix_mean
    agg$cov_norm <- 1000 - round(agg$pix_cov / max(agg$pix_cov) * 1000)
    
    
    g3 <- ggplot(agg) +
      geom_line( aes(x=barcenter, y=pix_mean), size=1, colour="#001100") +
      geom_point( aes(x=barcenter, y=pix_mean), size=1.6, colour="#001100") +
      geom_errorbar( aes(x=barcenter, ymin=pix_mean-pix_sd, ymax=pix_mean+pix_sd), size=0.6, colour="#001100", alpha=0.9) +
      theme_bw() + xlim(0, 1) +
      ylab(paste("Mean of average field", VI)) + xlab("Shadowiness - Shininess") +
      ggtitle(paste0(VI, " - ", buf_dist, " m"))
    plot(g3)
    
    g4 <- ggplot(agg) +
      geom_line( aes(x=barcenter, y=pix_cov), size=1, colour="#001100") +
      geom_point( aes(x=barcenter, y=pix_cov), size=1.6, colour="#001100") +
      theme_bw() + xlim(0, 1) +
      ylab(paste("CoV per field -", VI)) + xlab("Shadowiness - Shininess") +
      ggtitle(paste0(VI, " - ", buf_dist, " m"))
    plot(g4)
    
    # Define breaks NDVI and NDRE for vertical colouring
    
    g5 <- ggplot(agg) +
      geom_line( aes(x=barcenter, y=cov_norm), size=1, colour="#001100") +
      geom_point( aes(x=barcenter, y=cov_norm), size=1.6, colour="#001100") +
      theme_bw() + xlim(0, 1) + 
      ylab(paste("Weight")) + xlab("Shadowiness - Shininess") +
      ggtitle(paste0(VI, " - ", buf_dist, " m"))
    plot(g5)
    
  }
}
dev.off()

grad1 <- make_gradient( deg = 0,  n = 500, cols = c("#111111", "white"))
grad2 <- make_gradient( deg = 90,  n = 500, cols = c("green", "red" ))


# Compare results between different buffers
pdf(paste0(wd, "3_Output/02_Shadowiness_Shininess_Buffers.pdf"), width=10, height=8)

for(VI in VIs){
  
  for(COMBI in 1:length(combis)){
    filenames  <- list.files(paste0(wd, "2_Intermediate/"), pattern = "02_ZonalStats_")
    filenames  <- grep(VI, filenames, value = T)
    zonals     <- do.call(rbind, lapply(paste0(wd, "2_Intermediate/", filenames), read.csv))}
  zonals$weight       <- zonals$pix_nr / max(zonals$pix_nr)
  zonals$pix_cov      <- zonals$pix_sd / zonals$pix_mean
  zonals$barcenter    <- zonals$Break2 - (zonals$Break2 - zonals$Break1)/2
  zonals$barwidth     <- zonals$Break2 - zonals$Break1
  
  agg_avg      <- aggregator(zonals, "pix_mean", c("Buffer_m", "barcenter"), "mean")
  agg_sd       <- aggregator(zonals, "pix_mean", c("Buffer_m", "barcenter"), "sd")
  agg          <- cbind(agg_avg[,2:4], agg_sd$pix_mean)
  names(agg)   <- c("Buffer_m", "barcenter", "pix_mean", "pix_sd")
  agg$pix_cov  <- agg$pix_sd / agg$pix_mean
  agg$cov_norm <- 1000 - round(agg$pix_cov / max(agg$pix_cov) * 1000)
  
  g6 <- ggplot(agg) +
    annotation_custom(grob = grad2, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + 
    geom_line( aes(x=barcenter, y=pix_mean), size=0.8, colour="#001100") +
    geom_point(aes(x=barcenter, y=pix_mean), size=1.6, colour="#001100") +
    geom_errorbar( aes(x=barcenter, ymin=pix_mean-pix_sd, ymax=pix_mean+pix_sd),  size=0.6, colour="#001100", alpha=0.9) +
    theme_bw() + xlim(0, 1) + ylim(3000, 9500) +
    facet_wrap(~Buffer_m) +
    ylab(VI) + xlab("Shadowiness - Shininess") +
    ggtitle(VI)
  plot(g6)

  g7 <- ggplot(agg) +
    geom_line( aes(x=barcenter, y=pix_sd, group=Buffer_m, color=Buffer_m), size=0.8) +
    geom_point(aes(x=barcenter, y=pix_sd, group=Buffer_m, color=Buffer_m), size=1.6) +
    theme_bw() + xlim(0, 1) + ylim(0, 850) +
    ylab(paste0("Stdv ", VI)) + xlab("Shadowiness - Shininess") +
    ggtitle(VI)
  plot(g7)
  
  g8 <- ggplot(agg) +
    geom_line( aes(x=barcenter, y=cov_norm), size=0.8, colour="#001100") +
    geom_point(aes(x=barcenter, y=cov_norm), size=1.6, colour="#001100") +
    theme_bw() + xlim(0, 1) + ylim(0, 900) +
    facet_wrap(~Buffer_m) +
    ylab(paste("Weigh", VI)) + xlab("Shadowiness - Shininess") +
    ggtitle(VI)
  plot(g8)
  
}

dev.off()


# Idea:
# Remove first and last 2 percentiles
# Create graph (normal curve-ish) based on COV:
# Derive value at first percentile
# Derive (x and y) of maximum value
# Derive (x and y) of minimum value
# Create skewed histogram 
# Let curve stay flat after this minimum


# Another idea:
# For each percentile variable -> correlate with the ground data.
# See what correlates best
# 
# AND:
# Create formula including all
# 
# AND:
# Create correlation matrix









