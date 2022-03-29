# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Combine and clean data


############################
#### SET-UP ENVIRONMENT ####
############################

rm(list=ls()) 

# Set script and data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd     <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

SHAPE  <- "3.5"   # Defines what delineation will be chosen

# Load libraries 
library(purrr)    # map (for strsplit)
library(dplyr)

# Functions
source("functions/add_zero_to_single_characters_in_vector.R")
source("functions/sort_df.R")
source("functions/aggregator.R")
source("functions/omitter.R")

add_class <- function(df, val_col, class_col, vals){ # Add class according to MLSY Book (page 341)
  df[, class_col] <- NA
  df[  df[,val_col] <  vals[1], class_col]                               <- 1
  df[  df[,val_col] >= vals[1]  &  df[,val_col] < vals[2]  , class_col]  <- 1
  df[  df[,val_col] >= vals[2]  &  df[,val_col] < vals[3]  , class_col]  <- 0
  df[  df[,val_col] >= vals[3]  &  df[,val_col] < vals[4]  , class_col]  <- 0
  df[  df[,val_col] >= vals[4], class_col]                               <- 0
  return(df)
}

###########################
####     READ DATA     ####
###########################

# VI1: WEIGHING WDRVI
# df_vi1          <- read.csv(paste0(wd, "3_Output/11_WDRVIs_a010.csv"), stringsAsFactors = F)
df_vi1          <- read.csv(paste0(wd, "3_Output/11_WDRVIs_a010__3.5m.csv"), stringsAsFactors = F)
df_vi1$ID       <- paste0(df_vi1$PR_FI, "_", substr(df_vi1$palm_id,1,4), "P", substr(df_vi1$palm_id,5,6))
df_vi1[,c("RED_FI", "NIR_FI", "palm_id", "PR_FI")]   <- NULL


# VI2: PIXEL REMOVAL & PIXEL REMOVAL WITH CUT-OFF VALUES
df_vi2 <- read.csv(paste0(wd, "3_Output/08_spatial_variables_shape_", SHAPE, ".csv"), stringsAsFactors = F)
df_vi2[,c("PR","FI","NR","TRT")]            <- NULL
df_vi2[,grepl("G_cat", names(df_vi2))]      <- NULL


# VI3: NRED, NGRN & NNIR
df_vi3          <- read.csv(paste0(wd, "3_Output/11_band_norms.csv"), stringsAsFactors = F)
df_vi3$ID       <- paste0(df_vi3$PR_FI, "_", substr(df_vi3$palm_id,1,4), "P", substr(df_vi3$palm_id,5,6))
df_vi3[,c("palm_id", "PR_FI")] <- NULL


# VI4: WEIGHING NDRE
df_vi4          <- read.csv(paste0(wd, "3_Output/11_NDREs.csv"), stringsAsFactors = F)
df_vi4$ID       <- paste0(df_vi4$PR_FI, "_", substr(df_vi4$palm_id,1,4), "P", substr(df_vi4$palm_id,5,6))
df_vi4[,c("RED_FI", "NIR_FI", "palm_id")]   <- NULL


# VI5: 
shapes <- c("0.5", "1", "1.5", "2.5", "3.5", "1-0.5", "1.5-1", "2.5-1.5", "3.5-2.5")
for(i in 1:length(shapes)){
  
  SHAPE <- shapes[i]
  
  temp          <- read.csv(paste0(wd, "3_Output/08_spatial_variables_shape_", SHAPE, ".csv"), stringsAsFactors = F)
  names(temp)   <- gsub("_", paste0("_shp_", SHAPE), names(temp) )
  
  collies       <- temp %>% select(-contains("_"))
  shape         <- temp %>% select(contains(c('WDRV_', 'NDRE_')))
  shape         <- shape  %>% select(contains('cat1'))
  names(shape)  <- gsub("cat1", "", names(shape))
  
  if(i == 1){
    df_vi5  <- cbind(collies, shape)
  } else {
    df_vi5  <- cbind(df_vi5, shape)
  }
}
rm(temp, shape, collies)
df_vi5[,c("PR", "FI", "NR", "TRT")] <- NULL

# Pure nutrient concentrations and performance variables
df_na <- read.csv(paste0(wd, "1_Input/Onground/Year 2 VG and Nutrient content HS211012_pure_BB.csv"), fileEncoding="UTF-8-BOM", stringsAsFactors = F)

   # No calculation of entire kg in palm because:
       # Lower green fronds probably have little relation with the vegetation index signal anymore
       # The number of green fronds is partly determined by pruning practices
       # Hence introducing this variable confounds palm performance with management practice, how to conclude anything?

df_na$Fr_Tot_kg   <- 0.0781 * df_na$PCS + 0.3948
df_na$Fr_Lf_kg    <- 0.0305 * df_na$PCS + 0.1195
df_na$Fr_Ra_kg    <- 0.0327 * df_na$PCS + 0.0707
df_na$Fr_RP_kg    <- 0.0476 * df_na$PCS + 0.2753

# Nutrient amounts in leaflets
df_na$N_Lf_Fr_kg  <- df_na$Fr_Lf_kg * df_na$N_reg_Lf
df_na$P_Lf_Fr_kg  <- df_na$Fr_Lf_kg * df_na$P_reg_Lf
df_na$K_Lf_Fr_kg  <- df_na$Fr_Lf_kg * df_na$K_reg_Lf
df_na$Mg_Lf_Fr_kg <- df_na$Fr_Lf_kg * df_na$Mg_reg_Lf
df_na$Ca_Lf_Fr_kg <- df_na$Fr_Lf_kg * df_na$Ca_reg_Lf

# Nutrient amounts in rachis
df_na$N_Ra_Fr_kg  <- df_na$Fr_Ra_kg * df_na$N_reg_Ra
df_na$P_Ra_Fr_kg  <- df_na$Fr_Ra_kg * df_na$P_reg_Ra
df_na$K_Ra_Fr_kg  <- df_na$Fr_Ra_kg * df_na$K_reg_Ra
df_na$Mg_Ra_Fr_kg <- df_na$Fr_Ra_kg * df_na$Mg_reg_Ra
df_na$Ca_Ra_Fr_kg <- df_na$Fr_Ra_kg * df_na$Ca_reg_Ra

# Nutrient amounts in rachis + petiole
df_na$N_RP_Fr_kg  <- df_na$Fr_RP_kg * df_na$N_reg_Ra
df_na$P_RP_Fr_kg  <- df_na$Fr_RP_kg * df_na$P_reg_Ra
df_na$K_RP_Fr_kg  <- df_na$Fr_RP_kg * df_na$K_reg_Ra
df_na$Mg_RP_Fr_kg <- df_na$Fr_RP_kg * df_na$Mg_reg_Ra
df_na$Ca_RP_Fr_kg <- df_na$Fr_RP_kg * df_na$Ca_reg_Ra

# Nutrient amounts in entire frond
df_na$N_Tot_Fr_kg  <- df_na$N_Lf_Fr_kg  + df_na$N_RP_Fr_kg
df_na$P_Tot_Fr_kg  <- df_na$P_Lf_Fr_kg  + df_na$P_RP_Fr_kg
df_na$K_Tot_Fr_kg  <- df_na$K_Lf_Fr_kg  + df_na$K_RP_Fr_kg
df_na$Mg_Tot_Fr_kg <- df_na$Mg_Lf_Fr_kg + df_na$Mg_RP_Fr_kg
df_na$Ca_Tot_Fr_kg <- df_na$Ca_Lf_Fr_kg + df_na$Ca_RP_Fr_kg

# Production - total dry weight
df_na$FrProd_FrYr            <- df_na$FrProd_obs / df_na$Time_Period
df_na$Pa_Tot_kg_prod   <- df_na$Fr_Tot_kg  * df_na$FrProd_FrYr
df_na$Pa_Lf_kg_prod    <- df_na$Fr_Lf_kg   * df_na$FrProd_FrYr
df_na$Pa_Ra_kg_prod    <- df_na$Fr_Ra_kg   * df_na$FrProd_FrYr
df_na$Pa_RP_kg_prod    <- df_na$Fr_RP_kg   * df_na$FrProd_FrYr

# Production - nutrient amounts in all produced fronds, assuming 1) same concentration in all fronds and 2) same PCS (hence DW's) of all fronds.
df_na$N_Pa_Tot_kg_prod  <- df_na$N_Tot_Fr_kg  * df_na$FrProd_FrYr
df_na$P_Pa_Tot_kg_prod  <- df_na$P_Tot_Fr_kg  * df_na$FrProd_FrYr
df_na$K_Pa_Tot_kg_prod  <- df_na$K_Tot_Fr_kg  * df_na$FrProd_FrYr
df_na$Mg_Pa_Tot_kg_prod <- df_na$Mg_Tot_Fr_kg * df_na$FrProd_FrYr
df_na$Ca_Pa_Tot_kg_prod <- df_na$Ca_Tot_Fr_kg * df_na$FrProd_FrYr

# Production - idem, but only leaflets
df_na$N_Pa_Lf_kg_prod  <- df_na$N_Lf_Fr_kg  * df_na$FrProd_FrYr
df_na$P_Pa_Lf_kg_prod  <- df_na$P_Lf_Fr_kg  * df_na$FrProd_FrYr
df_na$K_Pa_Lf_kg_prod  <- df_na$K_Lf_Fr_kg  * df_na$FrProd_FrYr
df_na$Mg_Pa_Lf_kg_prod <- df_na$Mg_Lf_Fr_kg * df_na$FrProd_FrYr
df_na$Ca_Pa_Lf_kg_prod <- df_na$Ca_Lf_Fr_kg * df_na$FrProd_FrYr

# Make ready for merge
df_na$PR_FI    <- gsub("-", "_", substr( df_na$Field_ID,1, 5))
df_na          <- add_zero(df_na, "Sample_Point")
df_na$ID       <- paste0(df_na$PR_FI, "_", df_na$Treatment, "_P", df_na$NR_ID)
df_na          <- sort_df(df_na, "ID")
df_na[,c("Site", "Field_ID", "Farmer_Name", "Treatment", "Period", "Time_Period", "Sample_Point", "PR_FI")] <- NULL


# Area measurements
df_ar          <- read.csv(paste0(wd, "1_Input/Onground/Areastatement 211020.csv"), stringsAsFactors=F)
df_ar$PR_FI    <- substr(df_ar$Field_ID, 1, 5)
df_ar$PR_FI    <- gsub("-", "_", df_ar$PR_FI)
df_ar$TRT      <- df_ar$Treatment
df_ar$A_ha     <- df_ar$Ha
df_ar$nr_palms <- df_ar$Palms
df_ar          <- df_ar[,c("PR_FI", "TRT", "YOP", "A_ha", "nr_palms")]


# Nutrient concentrations individual trees
df_nu          <- read.csv(paste0(wd, "1_Input/Onground/sheet1_Intensive sampling 210310__intensive_sampling.csv"),      stringsAsFactors=F)
df_nu$DOM_nut  <- df_nu$DOS
df_nu$FI       <- substr(unlist(map(strsplit(df_nu$Field_ID, split = "-"), 2)),1,2)
df_nu          <- add_zero(df_nu, "Sample_point") # Creates column "NR"
df_nu$NR_ID    <- paste0("P", df_nu$NR_ID)
df_nu$NR_nut   <- df_nu$No
df_nu$ID       <- paste(df_nu$Site, df_nu$FI, df_nu$Treatment, df_nu$NR_ID, sep="_")
df_nu[,c("Site", "Farmer_ID", "DOS", "Lab_Ref", "Treatment", "PR",
         "NR_ID", "S_code", "No", "Field_ID", "Sample_point")] <- NULL
df_nu[,c("FI")] <- NULL 

df_nu <- add_class(df_nu, "N",  "N_class",  c(2.3,   2.4,   2.81,  3.00 ) ) # Adding nutrient classes to individual measurement data
df_nu <- add_class(df_nu, "P",  "P_class",  c(0.140, 0.150, 0.181, 0.250) ) # Adding nutrient classes to individual measurement data
df_nu <- add_class(df_nu, "K",  "K_class",  c(0.75,  0.90,  1.21,  1.60 ) ) # Adding nutrient classes to individual measurement data
df_nu <- add_class(df_nu, "Mg", "Mg_class", c(0.20,  0.25,  0.41,  0.70 ) ) # Adding nutrient classes to individual measurement data
df_nu <- add_class(df_nu, "Ca", "Ca_class", c(0.25,  0.50,  0.76,  1.00 ) ) # Adding nutrient classes to individual measurement data
df_nu <- add_class(df_nu, "B",  "B_class",  c(8,     15,    26,    40   ) ) # Adding nutrient classes to individual measurement data
df_nu$def_tot <- df_nu$N_class + df_nu$P_class + df_nu$K_class + df_nu$Mg_class + df_nu$Ca_class + df_nu$B_class


# SPAD measurements
df_sp          <- read.csv(paste0(wd, "1_Input/Onground/sheet2_Intensive sampling 210310__SPAD_Measurement_Year_2.csv"), stringsAsFactors=F)
df_sp$DOM_spad <- df_sp$DOM
df_sp$FI       <- substr( unlist(map(strsplit(df_sp$Field_ID, split = "-"), 2)), 1, 2)
df_sp          <- add_zero(df_sp, "Sample_point") # Creates column "NR"
df_sp$NR_ID    <- paste0("P", df_sp$NR_ID) 
df_sp$ID       <- paste(df_sp$Site, df_sp$FI, df_sp$Treatment, df_sp$NR_ID, sep="_")
df_sp$SPAD     <- rowMeans(df_sp[,c('SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')])
df_sp[,c("Site", "Field_ID", "Farmer", "Treatment", "DOM", "PR", "Sample_point", "FI", "NR_ID",
         'SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')] <- NULL


# Treatment boundary statistics
df_trt            <- read.csv("../Data/2_Intermediate/13_Field_Treatment_stats.csv", stringsAsFactors=F)
df_trt_avg        <- df_trt[df_trt$treatment == "mean",]
names(df_trt_avg) <- c("PR_FI", "metric", "TRT", paste0(names(df_trt_avg)[4:ncol(df_trt_avg)], "_avg"))
df_trt_sd         <- df_trt[df_trt$treatment == "sd",]
names(df_trt_sd)  <- c("PR_FI", "metric", "TRT", paste0(names(df_trt_sd)[4:ncol(df_trt_sd)], "_sd"))
df_trt            <- cbind(df_trt_avg[,1:3], df_trt_avg[4:ncol(df_trt_avg)], df_trt_sd[4:ncol(df_trt_sd)])
df_trt$metric     <- NULL

rm(df_trt_avg, df_trt_sd)


###########################
####    MERGE DATA     ####
###########################

# Merge preprocessed data (1)
m1             <- merge(df_sp, df_nu, by="ID", all.x=T, all.y=T) # Include all fields in merge
m2             <- merge(m1, df_na,    by="ID", all.x=T, all.y=F) # Include all fields in merge
m3             <- merge(m2, df_vi1,   by="ID", all.x=F, all.y=F) # Exclude fields in which no imagery was present 
m4             <- merge(m3, df_vi2,   by="ID")
m5             <- merge(m4, df_vi3,   by="ID")
m6             <- merge(m5, df_vi4,   by="ID")
m7             <- merge(m6, df_vi5,   by="ID")

# Add relevant and remove irrelevant variables to merge with last two datasets
vars_all         <- m7
vars_all$TRT     <- substr(vars_all$ID, 7,9)
vars_all$PR      <- substr(vars_all$PR_FI, 1,2)
vars_all$palm_id <- NULL

# Merge preprocessed data (2)
vars_all         <- merge(vars_all, df_trt, by=c("PR_FI", "TRT"), all.x=T, all.y=F)
vars_all         <- merge(vars_all, df_ar, by=c("PR_FI", "TRT") )


# Renaming some variables
vars_all$LA_F             <- vars_all$Leaf_Area_m2
vars_all$Leaf_Area_m2     <- NULL
vars_all$RL               <- vars_all$Rachis_length_cm
vars_all$Rachis_length_cm <- NULL

# Sort dataframe
vars_all         <- sort_df(vars_all, c("ID", "PR", "FI", "PR_FI", "TRT", "Farmer", "NR", "NR_nut", "DOM_nut", "DOM_spad", "YOP", "A_ha", "nr_palms", "RL", "LA_F"))

###########################
####  OUTLIER REMOVAL  ####
###########################

vars_all <- vars_all[vars_all$PR_FI != "RI_F1", ]         # Only REF, no BMP
vars_all <- vars_all[vars_all$PR_FI != "CK_F3", ]         # Weirdly high reflections -> low VI
vars_all <- vars_all[vars_all$PR_FI != "CK_F5", ]         # Weirdly high reflections (but already lost because no nutrients / SPAD measurements)
vars_all <- vars_all[vars_all$ID    != "JB_F6_REF_P06", ] # Other tree interferring
vars_all <- vars_all[vars_all$ID    != "CK_F8_BMP_P15", ] # Other tree interferring
vars_all <- vars_all[vars_all$ID    != "CK_F5_BMP_P14", ] # Partly blurred and clipped
vars_all <- vars_all[vars_all$ID    != "RI_F7_BMP_P08", ] # Weirdly high LA_F (18.67, while rest between 3.55 and 14.62)
vars_all <- vars_all[vars_all$ID    != "RI_F7_REF_P20", ] # Overhanging tree
vars_all <- vars_all[vars_all$ID    != "RI_F4_BMP_P09", ] # Palm tree seems absent..
vars_all <- vars_all[vars_all$ID    != "SS_F5_REF_P05", ] # Extreme high reflection compared to rest of trees
vars_all <- vars_all[vars_all$ID    != "JB_F3_REF_P07", ] # Extreme high backscatter on soil

#################################################
#######     RENAME SPATIAL VARIABLES     ########
#################################################

# cat1: Nothing removed
# cat2: Shines removed
# cat3: Shadows removed
# cat4: Shines and shadows removed

collies             <- names(vars_all)
collies             <- gsub("WDRV_", "WDRVI_", collies)
names(vars_all)     <- collies
vars_all$WDRVI      <- vars_all$WDRVI_cat1   # Change WDRVI to UnWeighed WDRVI
vars_all$NDRE       <- vars_all$NDRE_cat1    # Create NDRE_cat1 to UnWeighed NDRE
vars_all$MSAVI      <- vars_all$MSAV_cat1   # Change WDRVI to UnWeighed WDRVI
vars_all$NDVI       <- vars_all$NDVI_cat1    # Create NDRE_cat1 to UnWeighed NDRE

vars_all$WDRVI_LNSY <- vars_all$WDRVI_cat2
vars_all$WDRVI_LYSN <- vars_all$WDRVI_cat3
vars_all$WDRVI_LNSN <- vars_all$WDRVI_cat4

vars_all$NDRE_LNSY  <- vars_all$NDRE_cat2
vars_all$NDRE_LYSN  <- vars_all$NDRE_cat3
vars_all$NDRE_LNSN  <- vars_all$NDRE_cat4

vars_all$NDVI_LNSY <- vars_all$NDVI_cat2
vars_all$NDVI_LYSN <- vars_all$NDVI_cat3
vars_all$NDVI_LNSN <- vars_all$NDVI_cat4

vars_all$MSAVI_LNSY  <- vars_all$MSAVI_cat2
vars_all$MSAVI_LYSN  <- vars_all$MSAVI_cat3
vars_all$MSAVI_LNSN  <- vars_all$MSAVI_cat4

vars_all$WDRVI_cat1 <- NULL
vars_all$WDRVI_cat2 <- NULL
vars_all$WDRVI_cat3 <- NULL
vars_all$WDRVI_cat4 <- NULL
vars_all$NDRE_cat1  <- NULL
vars_all$NDRE_cat2  <- NULL
vars_all$NDRE_cat3  <- NULL
vars_all$NDRE_cat4  <- NULL
vars_all$MSAV_cat1  <- NULL
vars_all$MSAV_cat2  <- NULL
vars_all$MSAV_cat3  <- NULL
vars_all$MSAV_cat4  <- NULL
vars_all$NDVI_cat1  <- NULL
vars_all$NDVI_cat2  <- NULL
vars_all$NDVI_cat3  <- NULL
vars_all$NDVI_cat4  <- NULL

# Shape vegetation indices name change
names(vars_all) <- gsub("-", "_H_", names(vars_all))
names(vars_all) <- gsub("_shp_", "__", names(vars_all))

###############################
#######     PALMS     #########
###############################

write.csv(vars_all, "../Data/2_Intermediate/17_Preprocessed_palms.csv", row.names = F)

###############################
####### ACROSS FIELDS #########
###############################

vars_agg_PRFI              <- aggregator(vars_all, names(vars_all[,c(9:ncol(vars_all))]), c("PR_FI"), "mean")
vars_agg_PRFI$Function     <- NULL
vars_agg_PRFI_TRT          <- aggregator(vars_all, names(vars_all[,c(9:ncol(vars_all))]), c("PR_FI", "TRT"), "mean", all = T)
vars_agg_PRFI_TRT$Function <- NULL

# Add nutrient classes for field averages (D, L, O, H, E)
vars_agg_PRFI <- add_class(vars_agg_PRFI, "N_reg_Lf",  "N_reg_Lf_class",  c(2.3,   2.4,   2.81,  3.00 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI <- add_class(vars_agg_PRFI, "P_reg_Lf",  "P_reg_Lf_class",  c(0.140, 0.150, 0.181, 0.250) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI <- add_class(vars_agg_PRFI, "K_reg_Lf",  "K_reg_Lf_class",  c(0.75,  0.90,  1.21,  1.60 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI <- add_class(vars_agg_PRFI, "Mg_reg_Lf", "Mg_reg_Lf_class", c(0.20,  0.25,  0.41,  0.70 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI <- add_class(vars_agg_PRFI, "Ca_reg_Lf", "Ca_reg_Lf_class", c(0.25,  0.50,  0.76,  1.00 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI$def_reg_tot <- vars_agg_PRFI$N_reg_Lf_class + vars_agg_PRFI$P_reg_Lf_class + vars_agg_PRFI$K_reg_Lf_class + vars_agg_PRFI$Mg_reg_Lf_class + vars_agg_PRFI$Ca_reg_Lf_class

vars_agg_PRFI_TRT <- add_class(vars_agg_PRFI_TRT, "N_reg_Lf",  "N_reg_Lf_class",  c(2.3,   2.4,   2.81,  3.00 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI_TRT <- add_class(vars_agg_PRFI_TRT, "P_reg_Lf",  "P_reg_Lf_class",  c(0.140, 0.150, 0.181, 0.250) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI_TRT <- add_class(vars_agg_PRFI_TRT, "K_reg_Lf",  "K_reg_Lf_class",  c(0.75,  0.90,  1.21,  1.60 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI_TRT <- add_class(vars_agg_PRFI_TRT, "Mg_reg_Lf", "Mg_reg_Lf_class", c(0.20,  0.25,  0.41,  0.70 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI_TRT <- add_class(vars_agg_PRFI_TRT, "Ca_reg_Lf", "Ca_reg_Lf_class", c(0.25,  0.50,  0.76,  1.00 ) ) # Adding nutrient classes to individual measurement data
vars_agg_PRFI_TRT$def_reg_tot <- vars_agg_PRFI_TRT$N_reg_Lf_class + vars_agg_PRFI_TRT$P_reg_Lf_class + vars_agg_PRFI_TRT$K_reg_Lf_class + vars_agg_PRFI_TRT$Mg_reg_Lf_class + vars_agg_PRFI_TRT$Ca_reg_Lf_class

write.csv(vars_agg_PRFI,     "../Data/2_Intermediate/17_Preprocessed_PR_FI.csv",     row.names = F)
write.csv(vars_agg_PRFI_TRT, "../Data/2_Intermediate/17_Preprocessed_PR_FI_TRT.csv", row.names = F)







