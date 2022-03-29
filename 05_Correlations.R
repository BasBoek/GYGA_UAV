# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Inspect correlations between created spatial variables and ground measurement


rm(list=ls()) 

##########################################
#### CHECK WD's & ADAPT IF NECESSARY #####
##########################################


# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

# Load libraries 
library(purrr) # map (during strsplit)
library(corrplot)
library(lme4)
library(olsrr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)

source("functions/add_zero_to_single_characters_in_vector.R")
source("functions/sort_df.R")
source("functions/aggregator.R")

# Get p-value from lm (https://stackoverflow.com/questions/5587676/pull-out-p-values-and-r-squared-from-a-linear-regression)
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}


vi_files <- c(
  "3_Output/04_spatial_variables_NDVI_cut10.csv",
  "3_Output/04_spatial_variables_NDVI_cut30.csv",
  "3_Output/04_spatial_variables_NDRE_cut10.csv",
  "3_Output/04_spatial_variables_NDRE_cut30.csv"
)

# Read input data
descr_cols           <- read.csv(paste0(wd, "3_Output/04_spatial_variables_NDVI_cut10.csv"), stringsAsFactors = F)[,1:5]

df_NDVI_cut10        <- read.csv(paste0(wd, "3_Output/04_spatial_variables_NDVI_cut10.csv"), stringsAsFactors = F)
names(df_NDVI_cut10) <- paste0(names(df_NDVI_cut10), "_cut10")
df_NDVI_cut30        <- read.csv(paste0(wd, "3_Output/04_spatial_variables_NDVI_cut30.csv"), stringsAsFactors = F)
names(df_NDVI_cut30) <- paste0(names(df_NDVI_cut30), "_cut30")
df_NDRE_cut10        <- read.csv(paste0(wd, "3_Output/04_spatial_variables_NDRE_cut10.csv"), stringsAsFactors = F)
names(df_NDRE_cut10) <- paste0(names(df_NDRE_cut10), "_cut10")
df_NDRE_cut30        <- read.csv(paste0(wd, "3_Output/04_spatial_variables_NDRE_cut30.csv"), stringsAsFactors = F)
names(df_NDRE_cut30) <- paste0(names(df_NDRE_cut30), "_cut30")

df_vi <- cbind(descr_cols, df_NDVI_cut10[6:29], df_NDVI_cut30[6:29], df_NDRE_cut10[6:29], df_NDRE_cut30[6:29])
df_nu <- read.csv(paste0(wd, "1_Input/Onground/sheet1_Intensive sampling 210310__intensive_sampling.csv"),      stringsAsFactors=F)
df_sp <- read.csv(paste0(wd, "1_Input/Onground/sheet2_Intensive sampling 210310__SPAD_Measurement_Year_2.csv"), stringsAsFactors=F)

# Preprocess data before merge
df_nu$DOM_nut  <- df_nu$DOS
df_nu$FI       <- substr(unlist(map(strsplit(df_nu$Field_ID, split = "-"), 2)),1,2)
df_nu          <- add_zero(df_nu, "Sample_point") # Creates column "NR"
df_nu$NR_ID    <- paste0("P", df_nu$NR_ID)
df_nu$NR_nut   <- df_nu$No
df_nu$ID       <- paste(df_nu$Site, df_nu$FI, df_nu$Treatment, df_nu$NR_ID, sep="_")
df_nu[,c("Farmer", "Farmer_ID", "DOS", "Lab_Ref", "Treatment", "PR",
         "NR_ID", "FI", "Site", "S_code", "No", "Field_ID", "Sample_point")] <- NULL

df_sp$DOM_spad <- df_sp$DOM
df_sp$FI       <- substr( unlist(map(strsplit(df_sp$Field_ID, split = "-"), 2)), 1, 2)
df_sp          <- add_zero(df_sp, "Sample_point") # Creates column "NR"
df_sp$NR_ID    <- paste0("P", df_sp$NR_ID) 
df_sp$ID       <- paste(df_sp$Site, df_sp$FI, df_sp$Treatment, df_sp$NR_ID, sep="_")
df_sp$SPAD     <- rowMeans(df_sp[,c('SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')])
df_sp[,c("Site", "Field_ID", "Farmer", "Treatment", "DOM", "PR", "Sample_point", "FI", "NR_ID",
         'SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')] <- NULL

# Combine all data
vars        <- merge(df_sp, df_nu, all.x=T, all.y=T)
vars        <- merge(vars, df_vi, all.x=T, all.y=T)
vars$PR_FI  <- paste0(vars$PR, "_", vars$FI)
vars        <- vars[vars$PR_FI != "NA_NA",] # Excluding fields with no spatial info
vars$NR     <- paste0(substr(vars$TRT, 1,1), substr(vars$NR,2,3))
vars$NR     <- paste0(vars$PR_FI, "_", vars$NR)
vars        <- sort_df(vars, c("ID", "PR", "FI", "PR_FI", "TRT", "NR", "NR_nut", "DOM_nut", "DOM_spad"))

# Create correlation plots
var_cor        <- vars
var_cor        <- var_cor[complete.cases(var_cor), ]      # Remove rows containing NA
names(var_cor) <- paste0("v_", as.character(1:ncol(var_cor)))
corrplot(cor(var_cor[30:ncol(var_cor)]), method="circle")
corrplot(cor(var_cor[10:44]),            method="circle") # VI1 (NDRE)
corrplot(cor(var_cor[c(10:20, 69:92)]),  method="circle") # VI2 (NDVI)
corrplot(cor(var_cor[c(10:20, 93:116)]), method="circle") # VI2 (NDVI)


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
VIs     <- names(select(vars, contains("weight_no")))
VOIs    <- c("SPAD", "Ash", "N", "P", "K", "Mg", "Ca", "B", "LA_F", "TGF", "LA_Palm")
nr_rows <- length(FIs) * length(VIs) * length(VOIs)

stats[1:nr_rows,] <- NA

FI  <- "CK_F1"
VI  <- "NDRE___shape_1.5__0.5___weight_no_cut10"
VOI <- "SPAD"

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



