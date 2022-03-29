# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Inspect correlations between created spatial variables and ground measurement


rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

# Determine if writing pdf while running script (takes a while)
pdf_writing <- "yes"

shapes <- c("0.5", "1.5", "1.5-1", "1", "1-0.5", "2.5", "2.5-1.5", "3.5", "3.5-2.5")
SHAPE  <- "0.5"

# Load libraries 
library(purrr)    # map (for strsplit)
library(ggplot2)
library(dplyr)    # select

source("functions/add_zero_to_single_characters_in_vector.R")
source("functions/sort_df.R")
source("functions/aggregator.R")
source("functions/pvalue.R")
omitter <- function(data, desiredCols) { 
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
} # https://stackoverflow.com/questions/11254524/omit-rows-containing-specific-column-of-na

for(SHAPE in shapes){
  
  # Read vegetation indices 1
  df_vi1 <- read.csv(paste0(wd, "3_Output/11_WDRVIs_a010.csv"), stringsAsFactors = F)
  df_vi1$ID       <- paste0(df_vi1$PR_FI, "_", substr(df_vi1$palm_id,1,4), "P", substr(df_vi1$palm_id,5,6))
  df_vi1[,c("RED_FI", "NIR_FI", "palm_id")]             <- NULL
  
  # Read vegetation indices 2
  df_vi2 <- read.csv(paste0(wd, "3_Output/08_spatial_variables_shape_", SHAPE, ".csv"), stringsAsFactors = F)
  df_vi2[,c("PR","FI","NR","TRT")]            <- NULL
  #df_vi2[,grepl("_cat2", names(df_vi2))]      <- NULL
  #df_vi2[,grepl("_cat3", names(df_vi2))]      <- NULL
  df_vi2[,grepl("G_cat", names(df_vi2))]      <- NULL
  # names(df_vi2[,grepl("_cat", names(df_vi2))]) <- paste0(names(df_vi2[,grepl("_cat", names(df_vi2))]), "_35m_25m")
  
  
  # Read vegetation indices 3
  df_vi3          <- read.csv(paste0(wd, "3_Output/11_band_norms.csv"), stringsAsFactors = F)
  df_vi3$ID       <- paste0(df_vi3$PR_FI, "_", substr(df_vi3$palm_id,1,4), "P", substr(df_vi3$palm_id,5,6))
  df_vi3[,c("palm_id", "PR_FI")] <- NULL
  
  
  # Read nutrient measurement data
  df_nu  <- read.csv(paste0(wd, "1_Input/Onground/sheet1_Intensive sampling 210310__intensive_sampling.csv"),      stringsAsFactors=F)
  df_nu$DOM_nut  <- df_nu$DOS
  df_nu$FI       <- substr(unlist(map(strsplit(df_nu$Field_ID, split = "-"), 2)),1,2)
  df_nu          <- add_zero(df_nu, "Sample_point") # Creates column "NR"
  df_nu$NR_ID    <- paste0("P", df_nu$NR_ID)
  df_nu$NR_nut   <- df_nu$No
  df_nu$ID       <- paste(df_nu$Site, df_nu$FI, df_nu$Treatment, df_nu$NR_ID, sep="_")
  #df_nu$ID       <- str_replace(df_nu$ID, "_P", "_")
  df_nu[,c("Site", "Farmer_ID", "DOS", "Lab_Ref", "Treatment", "PR",
           "NR_ID", "S_code", "No", "Field_ID", "Sample_point")] <- NULL
  
  
  # Read SPAD measurements
  df_sp          <- read.csv(paste0(wd, "1_Input/Onground/sheet2_Intensive sampling 210310__SPAD_Measurement_Year_2.csv"), stringsAsFactors=F)
  df_sp$DOM_spad <- df_sp$DOM
  df_sp$FI       <- substr( unlist(map(strsplit(df_sp$Field_ID, split = "-"), 2)), 1, 2)
  df_sp          <- add_zero(df_sp, "Sample_point") # Creates column "NR"
  df_sp$NR_ID    <- paste0("P", df_sp$NR_ID) 
  df_sp$ID       <- paste(df_sp$Site, df_sp$FI, df_sp$Treatment, df_sp$NR_ID, sep="_")
  #df_sp$ID       <- str_replace(df_sp$ID, "_P", "_")
  df_sp$SPAD     <- rowMeans(df_sp[,c('SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')])
  df_sp[,c("Site", "Field_ID", "Farmer", "Treatment", "DOM", "PR", "Sample_point", "FI", "NR_ID",
           'SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')] <- NULL
  
  
  # Read treatment boundary statistics
  df_trt            <- read.csv("../Data/2_Intermediate/13_Field_Treatment_stats.csv", stringsAsFactors=F)
  df_trt_avg        <- df_trt[df_trt$treatment == "mean",]
  names(df_trt_avg) <- c("PR_FI", "metric", "TRT", paste0(names(df_trt_avg)[4:ncol(df_trt_avg)], "_avg"))
  df_trt_sd         <- df_trt[df_trt$treatment == "sd",]
  names(df_trt_sd)  <- c("PR_FI", "metric", "TRT", paste0(names(df_trt_sd)[4:ncol(df_trt_sd)], "_sd"))
  df_trt            <- cbind(df_trt_avg[,1:3], df_trt_avg[4:ncol(df_trt_avg)], df_trt_sd[4:ncol(df_trt_sd)])
  rm(df_trt_avg, df_trt_sd)
  
  # Read crown data
  # crowns_CK_F4 <- read.csv("../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_CK_F4.csv", stringsAsFactors=F)
  # crowns_JB_F7 <- read.csv("../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_JB_F7.csv", stringsAsFactors=F)
  # crowns_RI_F4 <- read.csv("../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_RI_F4.csv", stringsAsFactors=F)
  # crowns_SS_F5 <- read.csv("../Data/2_Intermediate/14_Crown_Delineations_NDVI_NDRE_SS_F5.csv", stringsAsFactors=F)
  
  # crowns_CK_F4$PR_FI <- "CK_F4"
  # crowns_JB_F7$PR_FI <- "JB_F7"
  # crowns_RI_F4$PR_FI <- "RI_F4"
  # crowns_SS_F5$PR_FI <- "SS_F5"
  # 
  # crowns_CK_F4$ID <- NA
  # crowns_JB_F7$ID <- paste0(crowns_JB_F7$PR_FI, crowns_JB_F7$Treat)
  # 
  # crowns_CK_F4 <- crowns_CK_F4[,c("PR_FI", "NDRE")]
  # crowns_JB_F7 <- crowns_JB_F7[,c("PR_FI", "NDRE")]
  # crowns_RI_F4 <- crowns_RI_F4[,c("PR_FI", "NDRE")]
  # crowns_SS_F5 <- crowns_SS_F5[,c("PR_FI", "NDRE")]
  # crowns_all   <- rbind(crowns_CK_F4, crowns_JB_F7, crowns_RI_F4, crowns_SS_F5)
  
  
  
  
  
  # Combine all data
  m1             <- merge(df_sp, df_nu, by="ID", all.x=T, all.y=T)
  m1[,c("PR", "FI")] <- NULL
  m2             <- merge(m1, df_vi1, by="ID", all.x=F, all.y=F)
  m3             <- merge(m2, df_vi2, by="ID")
  m4             <- merge(m3, df_vi3, by="ID")
  # md           <- setdiff(m3,m2)
  vars_all       <- m4
  vars_all$TRT   <- substr(vars_all$ID, 7,9)
  vars_all$PR    <- substr(vars_all$PR_FI, 1,2)
  vars_all$palm_id <- NULL
  vars_all       <- merge(vars_all, df_trt, by=c("PR_FI", "TRT"), all.x=T, all.y=F)
  vars_all$metric<- NULL
  vars_all       <- sort_df(vars_all, c("ID", "PR", "FI", "PR_FI", "TRT", "Farmer", "NR", "NR_nut", "DOM_nut", "DOM_spad"))
  
  # Add variable for NPK
  vars_all$N_norm <- (vars_all$N - mean(vars_all$N, na.rm=T)) / sd(vars_all$N, na.rm=T)
  vars_all$P_norm <- (vars_all$P - mean(vars_all$P, na.rm=T)) / sd(vars_all$P, na.rm=T)
  vars_all$K_norm <- (vars_all$K - mean(vars_all$K, na.rm=T)) / sd(vars_all$K, na.rm=T)
  vars_all$N_def  <- 0
  vars_all$P_def  <- 0
  vars_all$K_def  <- 0
  vars_all$N_def[vars_all$N_norm < - 1.5] <- 1
  vars_all$P_def[vars_all$P_norm < - 1.5] <- 1
  vars_all$K_def[vars_all$K_norm < - 1.5] <- 1
  vars_all$NPK_nm  <- vars_all$N_norm + vars_all$P_norm + vars_all$K_norm 
  vars_all$NPK_nm  <- (vars_all$NPK_nm - min(vars_all$NPK_nm, na.rm=T) )^ 2
  vars_all$NPK_def <- vars_all$N_def  + vars_all$P_def  + vars_all$K_def 
  vars_all$weight  <- ( (1 + vars_all$N_def) * vars_all$LA_F )^2
  vars_all$weight  <- ( (1 + vars_all$NPK_nm) * vars_all$LA_F )^2
  
  # Outlier removal
  vars_all <- vars_all[vars_all$PR_FI != "CK_F3", ]     # Weirdly high reflections -> low VI
  vars_all <- vars_all[vars_all$PR_FI != "CK_F5", ]     # Weirdly high reflections
  #vars <- vars[vars$PR != "RI", ] 
  vars_all <- vars_all[vars_all$ID != "JB_F6_REF_P06", ] # Other tree interferring
  vars_all <- vars_all[vars_all$ID != "CK_F8_BMP_P15", ] # Other tree interferring
  vars_all <- vars_all[vars_all$ID != "CK_F5_BMP_P14", ] # Partly blurred and clipped
  vars_all <- vars_all[vars_all$ID != "RI_F7_BMP_P08", ] # Weirdly high LA_F (18.67, while rest between 3.55 and 14.62)
  vars_all <- vars_all[vars_all$ID != "RI_F7_REF_P20", ] # Overhanging tree
  vars_all <- vars_all[vars_all$ID != "RI_F4_BMP_P09", ] # Palm tree seems absent..
  vars_all <- vars_all[vars_all$ID != "SS_F5_REF_P05", ] # Extreme high reflection compared to rest of trees
  vars_all <- vars_all[vars_all$ID != "JB_F3_REF_P07", ] # Extreme high backscatter on soil
  # names(vars_all)
  names(vars_all)
  
  
  # SELECTION - TESTING
  ##########################################
  
  #inspect <- vars[,c("ID", "VI", "VOI")]
  #vars    <- vars[vars$PR_FI == "SS_F7",]
  
  # lo_fit <- loess(VI ~ VOI, data=vars, span=1) 
  # lm_fit <- lm(VI    ~ VOI, data=vars) 
  # plot(vars$VOI, vars$VI, ylim=c(min(vars$VI, na.rm=T), max(vars$VI, na.rm=T)))
  # points(vars$VOI, predict(lo_fit), ylim=c(min(vars$VI, na.rm=T), max(vars$VI, na.rm=T)), pch=19)
  # 
  # r2_lo  <-  1 - sum(lo_fit$residuals^2) / sum((vars$VI - mean(vars$VI))^2 ) 
  # r2_lm  <-  1 - sum(lm_fit$residuals^2) / sum((vars$VI - mean(vars$VI))^2 ) 
  # r2_dif <-  r2_lo - r2_lm
  
  # VI       <- "NDRE_cat1"
  # VOI      <- "N"
  # vars     <- vars_all
  # vars$VI  <- vars[,VI]
  # vars$VOI <- vars[,VOI]
  # vars     <- omitter(vars, c("VI", "VOI"))
  # 
  # r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=LA_F), alpha=0.6, size=1/sqrt(vars$NPK_def+0.8)*3, na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI, "weighted")) +
  #   geom_point(aes(x=VI, y=VOI, color=LA_F), alpha=0.6, size=1/sqrt(vars$NPK_def+0.8)*3, na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI, weight=weight), formula=y~x, color="black", method="lm", se = T) 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=NPK_def), alpha=0.6, size=sqrt(vars$LA_F), na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=as.factor(NPK_def), shape=TRT), alpha=0.6, size=3, na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) 
  # xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) + xlab(VOI) + ylab(VI)
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, shape=TRT), alpha=0.6, size=3, na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) 
  # xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) + xlab(VOI) + ylab(VI)
  # 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=P), alpha=0.6, size=1/sqrt(vars$NPK_def+1)*3, na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=PR_FI, shape = TRT), alpha=0.5, size=2, na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) 
  # 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=log(NRED), size = log(NNIR)), alpha=0.5, na.rm=T) +
  #   scale_color_gradient2(midpoint=mean(log(vars$NRED)), low="#EEEEAA", mid="orange", high="red", space ="Lab" )
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=log(NNIR)), size=4.5, alpha=0.5, na.rm=T) +
  #   scale_color_gradient2(midpoint=mean(log(vars$NNIR)), low="lightblue", mid="purple", high="darkred", space ="Lab" )
  # 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, size=log(GREEN_PA), color=PR_FI), alpha=0.4, na.rm=T) 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, color=log(GREEN_PA)), size=4.5, alpha=0.4, na.rm=T) +
  #   scale_color_gradient2(midpoint=mean(log(vars$GREEN_PA)), low="brown", mid="orange", high="darkgreen", space ="Lab" )
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, size=log(GREEN_PA)), alpha=0.25, na.rm=T) +
  #   xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  #   facet_wrap(~ PR_FI)
  # 
  # 
  # 
  #geom_smooth(aes(x=VI, y=VOI, color=PR_FI), formula=y~x, method=lm, se = F) +
  # geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method=lm, se = T) 
  # 
  # ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  #   geom_point(aes(x=VI, y=VOI, alpha=0.2, color = TRT), size=2, na.rm=T) +
  #   geom_smooth(aes(x=VI, y=VOI, color=TRT), formula=y~x, method=lm, se = T) + 
  #   xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  #   # geom_mark_ellipse(aes(x=VI, y=VOI), geom="polygon",level=0.95, alpha=0.1) +
  #   facet_wrap(~ PR_FI)
  # 
  # SELECTION - ALL TO PDF
  ##########################################
  
  VOIs        <- names(vars_all[9:19])
  VIs         <- names(vars_all[20:ncol(vars_all)])
  # VOIs        <- names(vars_all[9:11])
  # VIs         <- names(vars_all[20:22])
  print(VOIs)
  print(VIs)
  
  # PDF 1
  
  if(pdf_writing == "yes"){
    
    
    pdf(paste0("../Data/3_Output/15_Cor_VOI_VI_scatter_", SHAPE, ".pdf"))
    for(VOI in VOIs){
      print(VOI)
      for(VI in VIs){
        vars     <- vars_all
        vars$VI  <- vars[,VI]
        vars$VOI <- vars[,VOI]
        #vars <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
        
        vars     <- omitter(vars, c("VI", "VOI"))
        r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
        
        gg <- ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI, "  r2 =", r2)) +
          geom_point(aes(x=VI, y=VOI, color=PR_FI, alpha=0.2, shape = TRT), size=2, na.rm=T) +
          geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method=lm, se = T) +
          xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) + xlab(VOI) + ylab(VI)
        
        plot(gg)
        
      }
    }
    dev.off()
    
    
    
    # PDF 2
    
    
    pdf("../Data/3_Output/15_Cor_VOI_VI_scatter_GRN_all_", SHAPE, ".pdf")
    VOIs        <- names(vars_all[9:19])
    VIs         <- names(vars_all[20:ncol(vars_all)])

    for(VOI in VOIs){
      print(VOI)
      for(VI in VIs){
        vars     <- vars_all
        vars$VI  <- vars[,VI]
        vars$VOI <- vars[,VOI]
        #vars <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
        
        vars     <- omitter(vars, c("VI", "VOI"))
        r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
        
        gg2 <- ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI, "  r2 =", r2)) +
          geom_point(aes(x=VI, y=VOI, color=log(GREEN_PA)), size=3, alpha=0.4, na.rm=T) +
          scale_color_gradient2(midpoint=mean(log(vars$GREEN_PA)), low="brown", mid="orange", high="darkgreen", space ="Lab" ) +
          xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) + xlab(VOI) + ylab(VI)
        
        
        plot(gg2)
        
      }
    }
    dev.off()
    
    
    # PDF 3
    pdf("../Data/3_Output/15_Cor_VOI_VI_scatter_GRN__PR_FI_", SHAPE, ".pdf")
    VOIs        <- names(vars_all[9:19])
    VIs         <- names(vars_all[20:ncol(vars_all)])

    for(VOI in VOIs){
      print(VOI)
      for(VI in VIs){
        vars     <- vars_all
        vars$VI  <- vars[,VI]
        vars$VOI <- vars[,VOI]
        #vars <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
        
        vars     <- omitter(vars, c("VI", "VOI"))
        r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
        
        gg3 <- ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI, "  r2 =", r2)) +
          geom_point(aes(x=VI, y=VOI, color=log(GREEN_PA)), alpha=0.35, na.rm=T) +
          scale_color_gradient2(midpoint=mean(log(vars$GREEN_PA)), low="brown", mid="orange", high="darkgreen", space ="Lab" ) +
          xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) + xlab(VOI) + ylab(VI) +
          facet_wrap(~ PR_FI)
        
        plot(gg3)
        
      }
    }
    dev.off()
    
  }
  
  ## AGGREGATION ON FIELD LEVEL
  
  names(vars_all)
  
  
  
  VI       <- "WDRVI"
  VOI      <- "P"
  vars     <- vars_all
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  
  
  vars_agg        <- aggregator(vars, c(VI, VOI), c("PR_FI", "TRT"), "mean")
  vars_agg$VI     <- vars_agg[,VI]
  vars_agg$VOI    <- vars_agg[,VOI]
  vars_agg_sd     <- aggregator(vars, c(VI, VOI), c("PR_FI", "TRT"), "sd")
  vars_agg_sd$VI  <- vars_agg_sd[,VI]
  vars_agg_sd$VOI <- vars_agg_sd[,VOI]
  
  VOIs        <- names(vars_all[ 9:19])
  VIs         <- names(vars_all[51:ncol(vars_all)])

  VOI <- "Ash"
  VI  <- "WDRVI_avg"
  if(pdf_writing == "yes"){
    
    pdf(paste0("../Data/3_Output/15_field_vars_", SHAPE, ".pdf"))
    
    for(VOI in VOIs){
      print(VOI)
      for(VI in VIs){
        vars     <- vars_all
        
        vars$VI  <- vars[,VI]
        vars$VOI <- vars[,VOI]
        vars     <- omitter(vars, c("VI", "VOI"))
        
        vars_agg        <- aggregator(vars, c(VI, VOI), c("PR_FI", "TRT"), "mean")
        vars_agg$VI     <- vars_agg[,VI]
        vars_agg$VOI    <- vars_agg[,VOI]
        
        vars_agg_sd     <- aggregator(vars, c(VI, VOI), c("PR_FI", "TRT"), "sd")
        # vars_agg_sd$VI  <- vars_agg_sd[,VI]
        vars_agg_sd$VOI <- vars_agg_sd[,VOI]
        
        r2 <- round(summary(lm(VI ~ VOI, vars_agg))$r.squared,2)
        
        gg4 <- ggplot(vars_agg) + ggtitle(paste(VOI, "-", VI, "  r2 =", r2)) + theme_bw() +
          geom_errorbar(aes(x = vars_agg[,"VI"],
                            y = vars_agg[,"VOI"],
                            ymin=vars_agg[,"VOI"]-vars_agg_sd[,"VOI"], 
                            ymax=vars_agg[,"VOI"]+vars_agg_sd[,"VOI"])) +
          xlim(min(vars_agg[,"VI"]), max(vars_agg[,"VI"])) + ylim(min(vars_agg[,"VOI"]), max(vars_agg[,"VOI"])) + xlab(VOI) + ylab(VI)
        
        geom_point(aes(x=VI, y=VOI, color=PR_FI, shape = TRT), size=4, na.rm=T)
        
        plot(gg4)
        
        
        
      }
    }
    dev.off()
  }
  
}




 # theme_bw() +
  # geom_smooth(aes(x=VOI, y=VI), size=1.5, linetype = "longdash", color="gray", formula=y~x, method=lm, se = F) +
  #geom_point(aes(x=N, y=GWWDRVI, color=PR_FI), na.rm=T)
# geom_text(aes(x=VOI, y=VI, label=NR, hjust = -0.2, vjust=-0.2), size=3, alpha=0.5) +
# xlab(VOI) +
# geom_smooth(aes(x=VOI, y=VI, color=TRT), formula=y~x, method=lm, se = F) +
# scale_color_manual(values = c("darkgreen", "orange")) +
# theme(legend.position="none", plot.title = element_text(size = 10, face = "bold"))


# options(warn=0) # with warnings
# options(warn=-1)

# Select variables and create models










  
stats <- data.frame(PR_FI=character(),
                 VI=character(),
                 VOI=character(),
                 r2_all=double(),
                 r2_out=double(),
                 p_all=double(),
                 p_out=double(),
                 outliers=character(),
                 stringsAsFactors=FALSE)

FIs     <- unique(vars$PR_FI)
VIs     <- c("WDRVI", "GWWDRVI")
VOIs    <- c("SPAD", "Ash", "N", "P", "K", "Mg", "Ca", "B", "LA_F", "TGF", "LA_Palm")
nr_rows <- length(FIs) * length(VIs) * length(VOIs)

stats[1:nr_rows,] <- NA

i <- 0
for(FI in FIs){
  for(VI in VIs){
    for(VOI in VOIs){
      
      i <- i + 1  # Row number to fill
      
      vars_FI     <- vars[vars$PR_FI == FI, ]
      vars_FI$VI  <- vars_FI[, VI]
      vars_FI$VOI <- vars_FI[, VOI]
      vars_FI     <- vars_FI[, c("ID", "TRT", "NR", "VOI", "VI")]

      stats$PR_FI[i] <- FI
      stats$VI[i]    <- VI
      stats$VOI[i]   <- VOI
      
      if(sum(!is.na(vars_FI[,"VOI"])) > 0){
        vars_lm     <- lm(VI ~ VOI, vars_FI)
        
        vars_FI$upper   <- as.vector(vars_lm$fitted.values + 2.5 * sd(vars_lm$residuals))
        vars_FI$lower   <- as.vector(vars_lm$fitted.values - 2.5 * sd(vars_lm$residuals))
        vars_FI$outlier <- "yes"
        vars_FI$outlier[vars_FI$VI < vars_FI$upper & vars_FI$VI > vars_FI$lower] <- "no"
        
        vars_FI_out     <- vars_FI[vars_FI$outlier == "no",]
        vars_lm_out     <- lm(VI ~ VOI, vars_FI_out)
        
        stats$r2_all[i] <- summary(vars_lm)$r.squared
        stats$p_all[i]  <- lmp(vars_lm)
        stats$r2_out[i] <- summary(vars_lm_out)$r.squared
        stats$p_out[i]  <- lmp(vars_lm_out)
        
        if(nrow(vars_FI_out) > 0){
          stats$outliers[i] <- paste(vars_FI[vars_FI$outlier == "yes", "NR"], collapse=" ")
        }
      } 
    }
  }
}

stats    <- stats[complete.cases(stats),]
outliers <- as.data.frame(table(unlist(strsplit(stats$outliers, " "))))

stats$index     <- unlist(lapply(strsplit(stats$VI, "___"), `[[`, 1))
stats$shape     <- unlist(lapply(strsplit(stats$VI, "___"), `[[`, 2))
stats$weightcut <- unlist(lapply(strsplit(stats$VI, "___"), `[[`, 3))

stats_agg       <- aggregator(stats, c("r2_all", "p_all", "r2_out", "p_out"), c("VOI", "index", "shape", "weightcut"), "mean")

write.csv(stats,     paste0(wd, "2_Intermediate/05_statistics_per_field.csv"),     row.names = F)
write.csv(stats_agg, paste0(wd, "2_Intermediate/05_statistics_per_parameter.csv"), row.names = F)


# g1 <- ggplot(vars_FI) +
#   ggtitle(paste0(FI, " - With outliers\n\n", VI), ) +
#   theme_bw() +
#   geom_smooth(aes(x=VOI, y=VI), size=1.5, linetype = "longdash", color="gray", formula=y~x, method=lm, se = F) +
#   geom_point(aes(x=VOI, y=VI, color=TRT, size=factor(outlier, ordered=T))) +
#   geom_text(aes(x=VOI, y=VI, label=NR, hjust = -0.2, vjust=-0.2), size=3, alpha=0.5) +
#   xlab(VOI) +
#   geom_smooth(aes(x=VOI, y=VI, color=TRT), formula=y~x, method=lm, se = F) +
#   scale_color_manual(values = c("darkgreen", "orange")) +
#   theme(legend.position="none", plot.title = element_text(size = 10, face = "bold"))
# 
# g2 <- ggplot(vars_FI_out) +
#   ggtitle(paste0(FI, " - Without outliers\n\n", VI), ) +
#   theme_bw() +
#   geom_smooth(aes(x=VOI, y=VI), size=1.5, linetype = "longdash",color="gray", formula=y~x, method=lm, se = F) +
#   geom_point(aes(x=VOI, y=VI, color=TRT, size=factor(outlier, ordered=T))) +
#   geom_text(aes(x=VOI, y=VI, label=NR, hjust = -0.2, vjust=-0.2), size=3, alpha=0.5) +
#   xlab(VOI) +
#   geom_smooth(aes(x=VOI, y=VI, color=TRT), formula=y~x, method=lm, se = F) +
#   scale_color_manual(values = c("darkgreen", "orange")) +
#   theme(legend.position="none", plot.title = element_text(size = 10, face = "bold"))
# 
# ggarrange(g1, g2, ncol = 2, nrow = 1)


# 1) Save r2 for all combinations
# 2) Determine which VI is best on average
# 3) Create above graphs for best VI
# 4) Share results

# 5) Only do circles, calculate donuts
  # 0.0m - 0.5m
  # 0.5m - 1.0m
  # 1.0m - 1.5m
  # 1.5m - 2.0m
  # 2.0m - 2.5m
  # 2.5m - 3.0m
  # 3.0m - 3.5m

# 6) Create 

# 6) Try some other method for weighing pixels:
  # Weighing differently (so determine what weight for what range for all)
  # Weighing based on different thing besides shadows and shines (greenness)




# 
# 
# 
# test1 <- lm(NDVI___shape_1__0.5___weight_yes_cut10 ~ TRT + SPAD + N + P + K + B + Ca + Ash + LA_Palm + PR_FI, data = vars)
# all1  <- ols_step_all_possible(test1)
# #best1 <- lmer(NDVI___shape_1__0.5___weight_yes_cut10 ~ TRT + SPAD + LA_Palm + (1 | PR_FI), data = vars)
# 
# test2 <- lm(NDVI___shape_1__0.5___weight_yes_cut30 ~ TRT + SPAD + N + P + K + B + Ca + Ash + LA_Palm + PR_FI, data = vars)
# all2  <- ols_step_all_possible(test2)
# 
# test3 <- lm(NDRE___shape_1__0.5___weight_yes_cut10 ~ TRT + SPAD + N + P + K + B + Ca + Ash + LA_Palm + PR_FI, data = vars)
# all3  <- ols_step_all_possible(test3)
# best3_lm  <- lm(NDRE___shape_1__0.5___weight_yes_cut10 ~ TRT + SPAD + PR_FI, data = vars)
# best3_lm  <- lm(NDRE___shape_1__0.5___weight_yes_cut10 ~ TRT + SPAD + K + PR_FI, data = vars)
# summary(best3_lm)
# best3_mlm <- lmer(NDRE___shape_1__0.5___weight_yes_cut10 ~ TRT + SPAD + (1 | PR_FI), data = vars)
# 
# 
# 
# test4 <- lm(NDRE___shape_1__0.5___weight_yes_cut30 ~ TRT + SPAD + N + P + K + B + Ca + Ash + LA_Palm + PR_FI, data = vars)
# all4  <- ols_step_all_possible(test4)
# 
# names(VI)
# 
# test <- lm(N ~ TRT + SPAD + P + K + B + Ca + Ash + LA_Palm + NDRE___shape_1__0.5___weight_yes_cut30 + PR_FI, data = vars)
# summary(test)
# 
# test <- vars[vars$PR_FI == "JB_F6",]
# plot(test$N, test$NDRE___shape_1.5__0.5___weight_no_cut10, col=as.factor(test$TRT))
# lm_test <- lm(NDRE___shape_1.5__0.5___weight_no_cut10 ~ N, data = test)
# test$NDRE___shape_1.5__0.5___weight_no_cut10 - 
# 
# sd(lm_test$residuals)
# 
# abline(lm_test, as.factor(test$TRT))
# 
# 
# 
# # plot(test$N[test$TRT=="BMP"], test$NDRE___shape_1__0.5___weight_no_cut30[test$TRT=="BMP"])
# # points(test$N[test$TRT=="TRT"], test$NDRE___shape_1__0.5___weight_no_cut30[test$TRT=="TRT"])
# 
# unique(vars$PR_FI)

