# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Inspect correlations between created spatial variables and ground measurement


rm(list=ls())  # Clean script <- clean mind

# Set Script and Data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.
wd <- "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/"

shapes <- c("0.5", "1.5", "1.5-1", "1", "1-0.5", "2.5", "2.5-1.5", "3.5", "3.5-2.5")
SHAPE  <- "3.5"

# Load libraries 
library(purrr)    # map (for strsplit)
library(ggplot2)
library(dplyr)    # select
library(gridExtra)
library(ggcorrplot)
library(maditr)

source("functions/add_zero_to_single_characters_in_vector.R")
source("functions/sort_df.R")
source("functions/aggregator.R")
source("functions/pvalue.R")

## FUNCTIONS

omitter <- function(data, desiredCols) { 
  completeVec <- complete.cases(data[, desiredCols])
  return(data[completeVec, ])
} # https://stackoverflow.com/questions/11254524/omit-rows-containing-specific-column-of-na

# Add class according to MLSY Book (page 341)
add_class <- function(df, val_col, class_col, vals){
  df[, class_col] <- NA
  df[  df[,val_col] <  vals[1], class_col]                               <- 8
  df[  df[,val_col] >= vals[1]  &  df[,val_col] < vals[2]  , class_col]  <- 4
  df[  df[,val_col] >= vals[2]  &  df[,val_col] < vals[3]  , class_col]  <- 2
  df[  df[,val_col] >= vals[3]  &  df[,val_col] < vals[4]  , class_col]  <- 1
  df[  df[,val_col] >= vals[4], class_col]                               <- 0
  return(df)
}

## READ DATA

# Read vegetation indices 1
df_vi1          <- read.csv(paste0(wd, "3_Output/11_WDRVIs_a010.csv"), stringsAsFactors = F)
df_vi1$ID       <- paste0(df_vi1$PR_FI, "_", substr(df_vi1$palm_id,1,4), "P", substr(df_vi1$palm_id,5,6))
df_vi1[,c("RED_FI", "NIR_FI", "palm_id")]   <- NULL

# Read vegetation indices 2
df_vi2 <- read.csv(paste0(wd, "3_Output/08_spatial_variables_shape_", SHAPE, ".csv"), stringsAsFactors = F)
df_vi2[,c("PR","FI","NR","TRT")]            <- NULL
df_vi2[,grepl("G_cat", names(df_vi2))]      <- NULL

# Read vegetation indices 3
df_vi3          <- read.csv(paste0(wd, "3_Output/11_band_norms.csv"), stringsAsFactors = F)
df_vi3$ID       <- paste0(df_vi3$PR_FI, "_", substr(df_vi3$palm_id,1,4), "P", substr(df_vi3$palm_id,5,6))
df_vi3[,c("palm_id", "PR_FI")] <- NULL

# Read nutrient amounts
df_na       <- read.csv(paste0(wd, "1_Input/Onground/Year 2 VG and Nutrient content HS211012.csv"), fileEncoding="UTF-8-BOM", stringsAsFactors = F)

df_na$PR_FI <- gsub("-", "_", substr( df_na$Field_ID,1, 5))
df_na       <- add_zero(df_na, "Sample_Point")
df_na$ID    <- paste0(df_na$PR_FI, "_", df_na$Treatment, "_P", df_na$NR_ID)
df_na       <- df_na[,c("ID",
                        "Frond_kg_Lf_prod", "Frond_kg_Ra_prod", "Frond_kg_RaPe_prod", "Frond_kg_Tot_prod", 
                        "N_Lf_prod", "P_Lf_prod", "K_Lf_prod", "Mg_Lf_prod", "Ca_Lf_prod", "N_RP_prod", "P_RP_prod", "K_RP_prod", "Mg_RP_prod", "Ca_RP_prod", "N_Fr_prod", "P_Fr_prod", "K_Fr_prod", "Mg_Fr_prod", "Ca_Fr_prod",
                        "N_Lf","P_Lf","K_Lf","Mg_Lf","Ca_Lf","N_RP","P_RP","K_RP","Mg_RP","Ca_RP","N_Fr","P_Fr","K_Fr","Mg_Fr","Ca_Fr", 
                        "N_reg_Lf","P_reg_Lf","K_reg_Lf","Mg_reg_Lf","Ca_reg_Lf","N_reg_Ra","P_reg_Ra","K_reg_Ra","Mg_reg_Ra","Ca_reg_Ra", 
                        "Rachis_length_cm", "PCS", "No_Leaflet", "Leaf_Area_m2", "TGF_Palm")]
df_na       <- as.data.frame(sapply(df_na, function(x) gsub(" ", "", x) ))


df_na[,2:ncol(df_na)]        <- as.data.frame(sapply(2:ncol(df_na), function(x) as.numeric(df_na[,x]) ))
df_na$N_div_Mg  <- df_na$N_reg_Lf/df_na$Mg_reg_Lf
# df_na$N_plus_Mg <- scale(df_nu$N) + scale(df_nu$Mg)
# df_na$PK        <- scale(df_nu$P) + scale(df_nu$K)
# df_na$N_div_Mg[is.infinite(df_nu$N_div_Mg)] <- NA


# Read area measurement
df_ar          <- read.csv(paste0(wd, "1_Input/Onground/Areastatement 211020.csv"), stringsAsFactors=F)
df_ar$PR_FI    <- substr(df_ar$Field_ID, 1, 5)
df_ar$PR_FI    <- gsub("-", "_", df_ar$PR_FI)
df_ar$TRT      <- df_ar$Treatment
df_ar$A_ha     <- df_ar$Ha
df_ar$nr_palms <- df_ar$Palms
df_ar          <- df_ar[,c("PR_FI", "TRT", "YOP", "A_ha", "nr_palms")]


# Read nutrient measurement data
df_nu          <- read.csv(paste0(wd, "1_Input/Onground/sheet1_Intensive sampling 210310__intensive_sampling.csv"),      stringsAsFactors=F)
df_nu$DOM_nut  <- df_nu$DOS
df_nu$FI       <- substr(unlist(map(strsplit(df_nu$Field_ID, split = "-"), 2)),1,2)
df_nu          <- add_zero(df_nu, "Sample_point") # Creates column "NR"
df_nu$NR_ID    <- paste0("P", df_nu$NR_ID)
df_nu$NR_nut   <- df_nu$No
df_nu$ID       <- paste(df_nu$Site, df_nu$FI, df_nu$Treatment, df_nu$NR_ID, sep="_")
df_nu[,c("Site", "Farmer_ID", "DOS", "Lab_Ref", "Treatment", "PR",
         "NR_ID", "S_code", "No", "Field_ID", "Sample_point")] <- NULL

df_nu <- add_class(df_nu, "N",  "N_class",  c(2.3,   2.4,   2.81,  3.00 ) )
df_nu <- add_class(df_nu, "P",  "P_class",  c(0.140, 0.150, 0.181, 0.250) )
df_nu <- add_class(df_nu, "K",  "K_class",  c(0.75,  0.90,  1.21,  1.60 ) )
df_nu <- add_class(df_nu, "Mg", "Mg_class", c(0.20,  0.25,  0.41,  0.70 ) )
df_nu <- add_class(df_nu, "Ca", "Ca_class", c(0.25,  0.50,  0.76,  1.00 ) )
df_nu <- add_class(df_nu, "B",  "B_class",  c(8,     15,    26,    40   ) )
df_nu$def_tot <- df_nu$N_class + df_nu$P_class + df_nu$K_class + df_nu$Mg_class + df_nu$Ca_class + df_nu$B_class

df_nu$def_class <- 2
df_nu$def_class[df_nu$def_tot < quantile(df_nu$def_tot, 0.40, na.rm=T)] <- 1
df_nu$def_class[df_nu$def_tot > quantile(df_nu$def_tot, 0.80, na.rm=T)] <- 3

df_nu$nut_lim <- "no_def"

trh <- 3
df_nu$nut_lim[df_nu$N_class > trh & df_nu$P_class < trh & df_nu$K_class < trh] <- "N_def"
df_nu$nut_lim[df_nu$N_class < trh & df_nu$P_class > trh & df_nu$K_class < trh] <- "P_def"
df_nu$nut_lim[df_nu$N_class < trh & df_nu$P_class < trh & df_nu$K_class > trh] <- "K_def"
df_nu$nut_lim[df_nu$N_class > trh & df_nu$P_class > trh & df_nu$K_class < trh] <- "NP_def"
df_nu$nut_lim[df_nu$N_class > trh & df_nu$P_class < trh & df_nu$K_class > trh] <- "NK_def"
df_nu$nut_lim[df_nu$N_class < trh & df_nu$P_class > trh & df_nu$K_class > trh] <- "PK_def"
df_nu$nut_lim[df_nu$N_class > trh & df_nu$P_class > trh & df_nu$K_class > trh] <- "NPK_def"



# temp <- df_nu[,c("N_class", "P_class", "K_class", "Mg_class", "Ca_class", "B_class")] 
# temp[temp == 8] <- "D"
# temp[temp == 4] <- "L"
# temp[temp == 2] <- "O"
# temp[temp == 1] <- "H"
# temp[temp == 0] <- "E"
# df_nu[,c("N_class", "P_class", "K_class", "Mg_class", "Ca_class", "B_class")]  <- temp
# 
# neworder <- c("D", "L", "O", "H", "E")
# temp     <- df_nu
# temp    <- arrange(transform(temp, N_class = factor(N_class, levels=neworder)), temp)




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

# Combine all data
m1             <- merge(df_sp, df_nu, by="ID", all.x=T, all.y=T)
m1[,c("PR", "FI")] <- NULL
m2             <- merge(m1, df_na, "ID", all.x=T, all.y=F)
m3             <- merge(m2, df_vi1, by="ID", all.x=F, all.y=F)
m4             <- merge(m3, df_vi2, by="ID")
m5             <- merge(m4, df_vi3, by="ID")

# md           <- setdiff(m3,m2)
vars_all       <- m5
vars_all$TRT   <- substr(vars_all$ID, 7,9)
vars_all$PR    <- substr(vars_all$PR_FI, 1,2)
vars_all$palm_id <- NULL
vars_all       <- merge(vars_all, df_trt, by=c("PR_FI", "TRT"), all.x=T, all.y=F)
vars_all       <- merge(vars_all, df_ar, by=c("PR_FI", "TRT") )

vars_all$metric<- NULL
vars_all       <- sort_df(vars_all, c("ID", "PR", "FI", "PR_FI", "TRT", "Farmer", "NR", "NR_nut", "DOM_nut", "DOM_spad", "YOP", "A_ha", "nr_palms"))


# Outlier removal
vars_all <- vars_all[vars_all$PR_FI != "RI_F1", ]      # Only REF, no BMP
vars_all <- vars_all[vars_all$PR_FI != "CK_F3", ]      # Weirdly high reflections -> low VI
vars_all <- vars_all[vars_all$PR_FI != "CK_F3", ]      # Weirdly high reflections -> low VI
vars_all <- vars_all[vars_all$PR_FI != "CK_F5", ]      # Weirdly high reflections (but already lost because no nutrients / SPAD measurements)
vars_all <- vars_all[vars_all$ID != "JB_F6_REF_P06", ] # Other tree interferring
vars_all <- vars_all[vars_all$ID != "CK_F8_BMP_P15", ] # Other tree interferring
vars_all <- vars_all[vars_all$ID != "CK_F5_BMP_P14", ] # Partly blurred and clipped
vars_all <- vars_all[vars_all$ID != "RI_F7_BMP_P08", ] # Weirdly high LA_F (18.67, while rest between 3.55 and 14.62)
vars_all <- vars_all[vars_all$ID != "RI_F7_REF_P20", ] # Overhanging tree
vars_all <- vars_all[vars_all$ID != "RI_F4_BMP_P09", ] # Palm tree seems absent..
vars_all <- vars_all[vars_all$ID != "SS_F5_REF_P05", ] # Extreme high reflection compared to rest of trees
vars_all <- vars_all[vars_all$ID != "JB_F3_REF_P07", ] # Extreme high backscatter on soil

vars_all$LA_F         <- vars_all$Leaf_Area_m2
vars_all$Leaf_Area_m2 <- NULL

vars_agg_PRFI_TRT <- aggregator(vars_all, names(vars_all[,c(9:ncol(vars_all))]), c("PR_FI", "TRT"), "mean", all = T)
vars_agg_PRFI     <- aggregator(vars_all, names(vars_all[,c(9:ncol(vars_all))]), c("PR_FI"), "mean", all = T)

### EXTRA ANALYSIS

# Calculate r values

sel <- c("WDRV_cat1", "N", "P", "K", "Mg", "Ca", "B", 
         "Rachis_length_cm", "LA_F", "Frond_kg_Lf_prod")

PR_FIs  <- unique(vars_all$PR_FI)

new_df <- data.frame(
  PR_FI=character(),
  var1=character(),
  var2=character(),
  r=numeric()
)

i <- 0
for(PR_FI in PR_FIs){
  for(VAR1 in sel){
    for(VAR2 in sel){
      
      df         <- vars_all[vars_all$PR_FI == PR_FI,]
      
      nas1 <- sum(is.na(df[,VAR1]))
      nas2 <- sum(is.na(df[,VAR2]))
      
      if(nas1 == 0 & nas2 == 0){
        
        i          <- i + 1
        
        lm         <- lm(df[,VAR1] ~ df[,VAR2])
        r          <- sqrt(summary(lm)$r.squared) * as.numeric(lm$coefficients[2]) / abs(lm$coefficients[2]) # direction negative or positive
        new_df[i,] <- c(PR_FI, VAR1, VAR2, r)
        
      }
    }
  }
}
new_df$r <- as.numeric(new_df$r)

r_avg    <- aggregator(new_df, "r", c("var1", "var2"), "mean")
r_avg    <- as.data.frame(dcast(r_avg, var1 ~ var2, value.var = "r"))
row.names(r_avg) <- r_avg$var1
r_avg$var1 <- NULL

r_std    <- aggregator(new_df, "r", c("var1", "var2"), "sd")
r_std    <- as.data.frame(dcast(r_std, var1 ~ var2, value.var = "r"))
row.names(r_std) <- r_avg$var1
r_std$var1 <- NULL

# Calculate p-values

vars <- vars_all

new_df <- data.frame(
  var1=character(),
  var2=character(),
  p=numeric()
)

i <- 0
for(VARS1 in sel){
  for(VARS2 in sel){
    i <- i + 1
    
    lm       <- lm(vars[,VARS1] ~ vars[,VARS2] + vars[,"PR_FI"])
    p_coef   <- summary(lm)$coefficients[2,4]

    new_df[i,] <- c(VARS1, VARS2, p_coef)
    
  }
}

p_val            <- as.data.frame(dcast(new_df, var1 ~ var2, value.var = "p"))
row.names(p_val) <- p_val$var1
p_val$var1       <- NULL
p_val            <- sapply(p_val, as.numeric)
rownames(p_val)  <- colnames(test)

test <- p_val
test <- test %>%
  slice(match(sel, var1))


# Correlation plot


ggcorrplot(r_avg,
           hc.order = F,
           lab = F,
           type = 'lower',
           tl.cex = 9,
           tl.srt = 90,
           p.mat=p_val,
           sig.level=0.05,
           insig="blank",
           pch = 1,
)




plot(vars_all$K_reg_Lf * vars_all$Ash, (vars_all$Mg_reg_Lf + vars_all$Ca_reg_Lf))
plot(vars_all$YOP,  vars_all$WDRVI_avg)
plot(vars_all$K_reg_Lf,  vars_all$K_reg_Ra)
plot(vars_all$P_reg_Lf,  vars_all$P_reg_Ra)
plot(vars_all$PCS,  vars_all$Rachis_length_cm)
plot(vars_all$PCS,  vars_all$LA_F)
plot(vars_all$PCS,  vars_all$LA_F)
plot(vars_all$WDRV_cat1,  vars_all$Mg_Lf_prod)
plot(vars_all$PCS,  vars_all$LA_Palm)
plot(vars_all$LA_F,  vars_all$LA_F)
plot(vars_all$N_reg_Lf,  vars_all$N_reg_Ra)
plot(vars_all$Mg_reg_Lf, vars_all$Mg_reg_Ra)
plot(vars_all$Ca_reg_Lf, vars_all$Ca_reg_Ra)

# Concentration vs production
plot(vars_all$N, vars_all$N_Lf_prod)


plot(vars_all$K_reg_Lf, (vars_all$Mg_reg_Lf + vars_all$Ca_reg_Lf))
plot(vars_all$WDRVI,    (vars_all$N_reg_Lf + vars_all$P_reg_Lf) * (vars_all$Mg_reg_Lf + vars_all$K_reg_Lf) )
plot(vars_all$WDRVI,     vars_all$Rachis_length_cm )
plot(vars_all$WDRVI,     vars_all$PCS)
plot(vars_all$WDRVI,     vars_all$Frond_kg_Tot_prod, col=as.factor(vars_all$PR_FI))
plot(vars_all$WDRVI,     vars_all$Rachis_length_cm * vars_all$PCS)
plot(vars_all$WDRVI,     vars_all$Rachis_length_cm * sqrt(vars_all$PCS))
plot(vars_all$NWWDRVI,   vars_all$Rachis_length_cm * sqrt(vars_all$PCS))
plot(vars_all$NDRE_cat1, vars_all$Rachis_length_cm * sqrt(vars_all$PCS))
plot(vars_all$K_reg_Lf,  vars_all$Mg_reg_Lf)

install.packages("MuMIn", dependencies = T)
library(MuMIn)
names(vars_agg_PRFI)
test <- vars_agg_PRFI_TRT[complete.cases(vars_agg_PRFI_TRT),]

blub <- lm(WDRV_cat1 ~ Mg_reg_Lf + N_reg_Lf + P_reg_Lf + Ca_reg_Lf + K_reg_Lf, data=vars_agg_PRFI_TRT)
blub <- lm(WDRV_cat1 ~ Rachis_length, data=vars_agg_PRFI_TRT)
summary(blub)
test <- dredge(blub)

plot(vars_agg_PRFI_TRT$Mg_reg_Ra, vars_agg_PRFI_TRT$WDRVI)
plot(vars_all$YOP, vars_all$Ash)
plot(vars_all$YOP, vars_all$K)
plot(vars_all$P_reg_Lf, vars_all$WDRV_cat1)
plot(vars_all$K_reg_Lf, vars_all$WDRV_cat1)

#  hist3D and ribbon3D with greyish background, rotated, rescaled,...
install.packages("plot3D")
library(plot3D)

x <- vars_all$YOP
y <- vars_all$LA_F
z <- vars_all$WDRV_cat1

# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)

# predict values on regular xy grid
grid.lines = 101
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy     <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)

# fitted points for droplines to surface
fitpoints <- predict(fit)

# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 1, 
          theta = 0, phi = 0, ticktype = "detailed",
          xlab = "x", ylab = "y", zlab = "WDRVI",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "FIT")

# scatter plot with regression plane
scatter3D(x, y, z, pch = 18, cex = 1, 
          theta = -10, phi = -30, ticktype = "detailed",
          xlab = "x", ylab = "y", zlab = "WDRVI",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "FIT")




rm(m1, m2, m3, m4, m5, trh)


# SELECTION - TESTING
##########################################

nut_lim_scatter <- function(VI, VOI, data){
  vars     <- data
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  
  ggplot(transform(vars, nut_lim=factor(nut_lim,levels=c("no_def", "N_def", "P_def", "K_def", "NP_def", "NK_def", "PK_def", "NPK_def")))) + 
    theme_bw() + xlab(VI) + ylab(VOI) +
    scale_color_manual(values=c( "#757575", "#fae105", "#0dd91b", "#09a5ed", "orange", "#6f0091", "#cc0e00", "black")) +
    geom_smooth(aes(x=VI, y=VOI ), formula=y~x, color="black", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=nut_lim), alpha=0.4, size=2.5, na.rm=T) +
    xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"]))
}

nut_lim_scatter("N", "Frond_kg_Lf_prod", vars_all)

PR_FI_scatter <- function(VI, VOI, data){
  vars     <- data
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  
  # g1 <- ggplot(vars) + 
  #   theme_bw() + xlab(VI) + ylab(VOI) +
  #   geom_smooth(aes(x=VI, y=VOI ), formula=y~x, color="black", method="lm", se = T) +
  #   geom_point(aes(x=VI, y=VOI, color = as.factor(nut_lim)), alpha=0.4, size=2.5, na.rm=T) +
  #   xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  #   facet_wrap(~PR_FI)
  
  g1 <- ggplot(vars) + 
    xlab(VI) + ylab(VOI) +
    geom_smooth(aes(x=VI, y=VOI ), formula=y~x, color="black", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color = TRT), alpha=0.4, size=2.5, na.rm=T) +
    xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
    facet_wrap(~PR_FI) + theme_bw()
  
  g1
  
  # g2 <- ggplot(vars) + 
  #   theme_bw() + xlab(VI) + ylab(VOI) +
  #   geom_smooth(aes(x=VI, y=VOI ), formula=y~x, color="black", method="lm", se = T) +
  #   geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.4, size=2.5, na.rm=T) +
  #   xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  #   facet_wrap(~nut_lim)
  # 
  # grid.arrange(g1, g2)
  
}

scatter2 <- function(VI, VOI, data){
  vars     <- data
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  
  g1 <- ggplot(vars) + 
    theme_bw() + xlab(VI) + ylab(VOI) +
    geom_smooth(aes(x=VI, y=VOI ), formula=y~x, color="black", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=nut_lim), alpha=0.4, size=2.5, na.rm=T) +
    xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) 
  
  g1
  
  
}


# VI versus VEGETATIVE PERFORMANCE
PR_FI_scatter("WDRV_cat1",   "Rachis_length_cm", vars_all)
PR_FI_scatter("WDRV_cat1",   "PCS", vars_all) 
PR_FI_scatter("WDRV_cat1",   "LA_F", vars_all) 
PR_FI_scatter("WDRV_cat1",   "Frond_kg_Lf_prod", vars_all) # Derivative of PCS

summary(lm(Rachis_length_cm ~ WDRV_cat1 + PR_FI, vars_all)) # super significant
summary(lm(PCS ~ WDRV_cat1              + PR_FI, vars_all)) # super significant
summary(lm(Frond_kg_Lf_prod ~ WDRV_cat1 + PR_FI, vars_all)) # super significant

# Conclusions:
# - All show that there is a relation between vegetative performance and WDRVI
# - Visual interpretation: some fields clearly responding to BMP treatment
# - Hence, in-field relative differences can be indication for spread of nutrients

# VEGETATIVE PERFORMANCE versus NUTRIENT CONCENTRATIONS - KG Frond Production / year
PR_FI_scatter("N",   "Frond_kg_Lf_prod", vars_all)
PR_FI_scatter("P",   "Frond_kg_Lf_prod", vars_all)
PR_FI_scatter("K",   "Frond_kg_Lf_prod", vars_all)
PR_FI_scatter("Mg",  "Frond_kg_Lf_prod", vars_all)
PR_FI_scatter("Ca",  "Frond_kg_Lf_prod", vars_all)
PR_FI_scatter("B",   "Frond_kg_Lf_prod", vars_all)

summary(lm(Frond_kg_Lf_prod ~ N  + PR_FI, vars_all)) # A bit positive
summary(lm(Frond_kg_Lf_prod ~ P  + PR_FI, vars_all)) # A bit positive
summary(lm(Frond_kg_Lf_prod ~ K  + PR_FI, vars_all)) # A bit positive
summary(lm(Frond_kg_Lf_prod ~ Mg + PR_FI, vars_all)) # ns.
summary(lm(Frond_kg_Lf_prod ~ Ca + PR_FI, vars_all)) # A bit negative
summary(lm(Frond_kg_Lf_prod ~ B  + PR_FI, vars_all)) # ns.

PR_FI_scatter("N",   "Rachis_length_cm", vars_all)
PR_FI_scatter("P",   "Rachis_length_cm", vars_all)
PR_FI_scatter("K",   "Rachis_length_cm", vars_all)
PR_FI_scatter("Mg",  "Rachis_length_cm", vars_all)
PR_FI_scatter("Ca",  "Rachis_length_cm", vars_all)
PR_FI_scatter("B",   "Rachis_length_cm", vars_all)
summary(lm(Frond_kg_Lf_prod ~ N  + PR_FI, vars_all)) #
summary(lm(Frond_kg_Lf_prod ~ P  + PR_FI, vars_all)) #
summary(lm(Frond_kg_Lf_prod ~ K  + PR_FI, vars_all)) #
summary(lm(Frond_kg_Lf_prod ~ Mg + PR_FI, vars_all)) #
summary(lm(Frond_kg_Lf_prod ~ Ca + PR_FI, vars_all)) #
summary(lm(Frond_kg_Lf_prod ~ B  + PR_FI, vars_all)) #

# VI versus NUTRIENT CONCENTRATIONS
PR_FI_scatter("N",   "WDRV_cat1", vars_all) #
PR_FI_scatter("P",   "WDRV_cat1", vars_all) # 
PR_FI_scatter("K",   "WDRV_cat1", vars_all) # 
PR_FI_scatter("Mg",  "WDRV_cat1", vars_all) # 
PR_FI_scatter("Ca",  "WDRV_cat1", vars_all) # 
PR_FI_scatter("B",   "WDRV_cat1", vars_all) # 

PR_FI_scatter("Ash", "Frond_kg_Lf_prod", vars_all) # A bit negative

summary(lm(WDRV_cat1 ~ N  + PR_FI, vars_all))
summary(lm(WDRV_cat1 ~ P  + PR_FI, vars_all))
summary(lm(WDRV_cat1 ~ K  + PR_FI, vars_all))
summary(lm(WDRV_cat1 ~ Mg + PR_FI, vars_all)) 
summary(lm(WDRV_cat1 ~ Ca + PR_FI, vars_all)) 
summary(lm(WDRV_cat1 ~ B  + PR_FI, vars_all)) 

# ACROSS FIELDS
scatter2("WDRV_cat1", "N_Lf_prod", vars_agg_PRFI)
scatter2("WDRV_cat1", "P_Lf_prod", vars_agg_PRFI)
scatter2("WDRV_cat1", "K_Lf_prod", vars_agg_PRFI)
scatter2("WDRV_cat1", "Mg_Lf_prod",vars_agg_PRFI)
scatter2("WDRV_cat1", "Ca_Lf_prod",vars_agg_PRFI)

# Conclusions
# - BMP / REF effect field dependent
# - In general 









VI       <- "WDRVI_cat1"
VOI      <- "TGF"
VI       <- "LA_Palm"
VOI      <- "def_tot"

vars     <- vars_all
vars$VI  <- vars[,VI]
vars$VOI <- vars[,VOI]
vars     <- omitter(vars, c("VI", "VOI"))

ggplot(transform(vars, nut_lim=factor(nut_lim,levels=c("no_def", "N_def", "P_def", "K_def", "NP_def", "NK_def", "PK_def", "NPK_def")))) + 
  theme_bw() + xlab(VI) + ylab(VOI) +
  scale_color_manual(values=c( "#757575", "#fae105", "#0dd91b", "#09a5ed", "orange", "#6f0091", "#cc0e00", "black")) +
  geom_smooth(aes(x=VI, y=VOI ), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=nut_lim), alpha=0.4, size=2.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ PR_FI)


r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)


ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_point(aes(x=VI, y=VOI, color=def_tot), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ PR_FI)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ def_class)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_point(aes(x=VI, y=VOI, color=K_class), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ PR_FI)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ K_class)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_point(aes(x=VI, y=VOI, color=P_class), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ PR_FI)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ P_class)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_point(aes(x=VI, y=VOI, color=N_class), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ PR_FI)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  facet_wrap(~ nut_lim)

ggplot(transform(vars, nut_lim=factor(nut_lim,levels=c("N_lim", "P_lim", "K_lim", "NP_lim", "NK_lim", "PK_lim", "NPK_lim", "no_lim")))) + 
  theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) 
facet_wrap(~ nut_lim)

# Nutrient amounts -- ALL PALMS

VI       <- "WDRV_cat1"
VOI      <- "PK"

vars     <- vars_all
vars$VI  <- vars[,VI]
vars$VOI <- vars[,VOI]
vars     <- omitter(vars, c("VI", "VOI"))
vars     <- omitter(vars, c("VI", "VOI"))
unique(vars$PR_FI)

ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=TRT), size=2, alpha=0.3, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  xlab(VI) + ylab(VOI) +
  facet_wrap(~ PR_FI)

plot(vars$N_div_Mg_Fr/ vars$N_plus_Mg_Fr, vars$B)
plot(vars$LA_F, vars$LA_Palm)

vars     <- vars[vars$PR_FI == "SS_F1",]
ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=TRT), size=2, alpha=0.3, na.rm=T) +
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  xlab(VI) + ylab(VOI) + geom_text(aes(VI, VOI, label=ID), hjust=0, vjust=-2, size=2) 
facet_wrap(~ TGF_Palm)

plot(vars$TGF, vars$TGF_Palm)

print("stop")



# Nutrient amounts -- FIELD AVERAGES
#################################################################################

VI       <- "LA_F"
VOI      <- "TGF_Palm"
vars     <- vars_agg_PRFI
vars$VI  <- vars[,VI]
vars$VOI <- vars[,VOI]
vars     <- omitter(vars, c("VI", "VOI"))
vars     <- omitter(vars, c("VI", "VOI"))

LM       <- lm(VI ~ VOI , data=vars)
r2       <- round(summary(LM)$r.squared, 2)
pval     <- round(lmp(LM), 4)

# vars$PR <- substr(vars$PR_FI, 1,2)
# test <- aggregator(vars, c("YOP", "TGF_Palm"), "PR", c("mean", "sd"))

ggplot(vars) + theme_bw() + 
  geom_smooth(aes(x=VI, y=VOI), alpha=0.05, formula=y~x, color="gray", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=YOP), size=3, alpha=0.6, na.rm=T) +
  geom_text(aes(x=VI, y=VOI, label=PR_FI), hjust=-0.1, vjust=-0.5, position = position_dodge(width=0.9),  size=2.4, alpha=0.5) +
  xlab(VI) + ylab(VOI) +
  theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 7)) +
  guides(color=guide_legend(ncol=1, bycol=TRUE))

ggplot(vars) + theme_bw() + 
  geom_smooth(aes(x=VI, y=VOI), alpha=0.05, formula=y~x, color="gray", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=LA_F, shape=TRT), size=3, alpha=0.6, na.rm=T) +
  geom_text(aes(x=VI, y=VOI, label=PR_FI), hjust=-0.1, vjust=-0.5, position = position_dodge(width=0.9),  size=2.4, alpha=0.5) +
  xlab(VI) + ylab(VOI) +
  theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 7)) +
  guides(color=guide_legend(ncol=1, bycol=TRUE))

pdf("../Data/3_Output/16_WDRVI_VGMs.pdf", width=9, height=7)
VOIs_orig <- c("PCS", "Rachis_length_cm", "LA_F")
for(i in 1:length(VOIs_orig)){
  
  VOIs     <- VOIs_orig
  
  VI       <- "WDRV_cat1"
  VOI      <- VOIs[i]
  
  vars     <- vars_agg_PRFI_TRT
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  vars     <- omitter(vars, c("VI", "VOI"))
  
  xminval  <- min(vars[,"VI"])
  xmaxval  <- max(vars[,"VI"])
  yminval  <- min(vars[,"VOI"])
  ymaxval  <- max(vars[,"VOI"])
  
  LM       <- lm(VI ~ VOI , data=vars)
  r2       <- round(summary(LM)$r.squared, 2)
  pval     <- round(lmp(LM), 4)
  
  gg <- ggplot(vars) + theme_bw() + ggtitle(paste0("r2 = ", r2)) +
    geom_smooth(aes(x=VI, y=VOI), alpha=0.1, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=PR_FI, shape=TRT), size=2, alpha=0.8, na.rm=T) +
    xlim(xminval, xmaxval) +ylim(yminval, ymaxval) +
    xlab("WDRVI") + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  print(gg)
}
dev.off()



vars         <- vars_agg_PRFI_TRT
vars$palm_ha <- vars$nr_palms   / vars$A_ha
vars$N_ha    <- vars$N_Fr_prod  * vars$palm_ha
vars$P_ha    <- vars$P_Fr_prod  * vars$palm_ha
vars$K_ha    <- vars$K_Fr_prod  * vars$palm_ha
vars$Mg_ha   <- vars$Mg_Fr_prod * vars$palm_ha
vars$Ca_ha   <- vars$Ca_Fr_prod * vars$palm_ha
# VI           <- "WDRVI_avg"
# VOI          <- "Ca_ha"
pdf("../Data/3_Output/16_Nutrients_per_hectare", width=8, height=6)
VOIs_orig    <- c("palm_ha", "N_ha", "P_ha", "K_ha", "Mg_ha", "Ca_ha")
for(i in 1:length(VOIs_orig)){
  VOIs       <- VOIs_orig
  VOI        <- VOIs[i]
  
  vars$VI    <- vars[,VI]
  vars$VOI   <- vars[,VOI]
  vars       <- omitter(vars, c("VI", "VOI"))
  vars       <- omitter(vars, c("VI", "VOI"))
  
  
  gg <- ggplot(vars) + theme_bw() + 
    geom_smooth(aes(x=VI, y=VOI), alpha=0.05, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=LA_F, shape=TRT), size=3, alpha=0.6, na.rm=T) +
    geom_text(aes(x=VI, y=VOI, label=PR_FI), hjust=-0.1, vjust=-0.5, position = position_dodge(width=0.9),  size=2.4, alpha=0.5) +
    xlab(VI) + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  print(gg)
  
}
dev.off()

vars         <- vars_agg_PRFI_TRT
vars$palm_ha <- vars$nr_palms   / vars$A_ha
VI       <- "palm_ha"
VOI      <- "LA_F"

vars$VI  <- vars[,VI]
vars$VOI <- vars[,VOI]
vars     <- omitter(vars, c("VI", "VOI"))
vars     <- omitter(vars, c("VI", "VOI"))

test     <- vars[,c("PR_FI", "YOP")]

ggplot(vars) + theme_bw() + 
  geom_smooth(aes(x=VI, y=VOI), alpha=0.05, formula=y~x, color="gray", method="lm", se = T) +
  geom_point(aes(x=VI, y=VOI, color=YOP, shape=TRT), size=3, alpha=0.6, na.rm=T) +
  geom_text(aes(x=VI, y=VOI, label=PR_FI), hjust=-0.1, vjust=-0.5, position = position_dodge(width=0.9),  size=2.4, alpha=0.5) +
  xlab(VI) + ylab(VOI) +
  theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
        legend.title = element_text(size = 10, face = "bold"),
        axis.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 7)) +
  guides(color=guide_legend(ncol=1, bycol=TRUE))
print("")





pdf("../Data/3_Output/16_Nutrient_contents.pdf", width=15, height=11)
for(i in 1:length(VOIs_orig)){
  
  VOIs     <- VOIs_orig
  
  VI       <- "WDRV_cat1"
  VOI      <- VOIs[i]
  
  names(vars)
  vars     <- vars_agg_PRFI_TRT
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  vars     <- omitter(vars, c("VI", "VOI"))
  
  xminval  <- min(vars[,"VI"])
  xmaxval  <- max(vars[,"VI"])
  yminval  <- min(vars[,"VOI"])
  ymaxval  <- max(vars[,"VOI"])
  
  LM <- lm(VI ~ VOI , data=vars)
  r2   <- round(summary(LM)$r.squared,2)
  pval <- round(lmp(LM),4)
  
  gg1 <- ggplot(vars) + theme_bw() + ggtitle(paste0("Total Content: r2 = ", r2)) +
    geom_smooth(aes(x=VI, y=VOI), alpha=0.1, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=PR_FI, shape=TRT), size=2, alpha=0.8, na.rm=T) +
    xlim(xminval, xmaxval) +ylim(yminval, ymaxval) +
    xlab("WDRVI") + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  
  vars     <- vars_agg_PRFI
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  vars     <- omitter(vars, c("VI", "VOI"))
  
  LM <- lm(VI ~ VOI , data=vars)
  r2   <- round(summary(LM)$r.squared,2)
  pval <- round(lmp(LM),4)
  
  gg2 <- ggplot(vars) + theme_bw() + ggtitle(paste0("Total Content: r2 = ", r2, ", p = ", pval)) +
    geom_smooth(aes(x=VI, y=VOI), alpha=0.1, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=PR_FI), size=2, alpha=0.8, na.rm=T) +
    xlim(xminval, xmaxval) +ylim(yminval, ymaxval) +
    xlab("WDRVI") + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  
  
  VOIs <- paste0(VOIs, "_prod")
  VOI <- VOIs[i]
  
  
  names(vars)
  vars     <- vars_agg_PRFI_TRT
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  vars     <- omitter(vars, c("VI", "VOI"))
  
  xminval  <- min(vars[,"VI"])
  xmaxval  <- max(vars[,"VI"])
  yminval  <- min(vars[,"VOI"])
  ymaxval  <- max(vars[,"VOI"])
  
  LM <- lm(VI ~ VOI , data=vars)
  r2   <- round(summary(LM)$r.squared,2)
  pval <- round(lmp(LM),4)
  
  gg3 <- ggplot(vars) + theme_bw() + ggtitle(paste0("Production: r2 = ", r2)) +
    geom_smooth(aes(x=VI, y=VOI), alpha=0.1, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=PR_FI, shape=TRT), size=2, alpha=0.8, na.rm=T) +
    xlim(xminval, xmaxval) +ylim(yminval, ymaxval) +
    xlab("WDRVI") + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  
  vars     <- vars_agg_PRFI
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  vars     <- omitter(vars, c("VI", "VOI"))
  
  LM <- lm(VI ~ VOI , data=vars)
  r2   <- round(summary(LM)$r.squared,2)
  pval <- round(lmp(LM),4)
  
  gg4 <- ggplot(vars) + theme_bw() + ggtitle(paste0("Production: r2 = ", r2, ", p = ", pval)) +
    geom_smooth(aes(x=VI, y=VOI), alpha=0.1, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=PR_FI), size=2, alpha=0.8, na.rm=T) +
    xlim(xminval, xmaxval) +ylim(yminval, ymaxval) +
    xlab("WDRVI") + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  grid.arrange(gg2, gg1, gg4, gg3, ncol=2)
  
}
dev.off()



# Nutrient concentrations -- FIELD AVERAGES: pooled samples
#################################################################################


VOIs     <- c("N_reg_Lf","P_reg_Lf","K_reg_Lf","Mg_reg_Lf","Ca_reg_Lf","N_reg_Ra","P_reg_Ra","K_reg_Ra","Mg_reg_Ra","Ca_reg_Ra")


pdf("../Data/3_Output/16_Nutrient_concentrations_pooled_samples.pdf", width=15, height=5.5)
for(i in 1:length(VOIs)){
  
  VI       <- "WDRV_cat1"
  VOI      <- VOIs[i]
  
  vars     <- vars_agg_PRFI_TRT
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  vars     <- omitter(vars, c("VI", "VOI"))
  
  xminval  <- min(vars[,"VI"])
  xmaxval  <- max(vars[,"VI"])
  yminval  <- min(vars[,"VOI"])
  ymaxval  <- max(vars[,"VOI"])
  
  LM <- lm(VI ~ VOI , data=vars)
  r2   <- round(summary(LM)$r.squared,2)
  pval <- round(lmp(LM),4)
  
  gg1 <- ggplot(vars) + theme_bw() + ggtitle(paste0("Nutrient concentrations: r2 = ", r2)) +
    geom_smooth(aes(x=VI, y=VOI), alpha=0.1, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=PR_FI, shape=TRT), size=2, alpha=0.8, na.rm=T) +
    xlim(xminval, xmaxval) +ylim(yminval, ymaxval) +
    xlab("WDRVI") + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  vars     <- vars_agg_PRFI
  vars$VI  <- vars[,VI]
  vars$VOI <- vars[,VOI]
  vars     <- omitter(vars, c("VI", "VOI"))
  vars     <- omitter(vars, c("VI", "VOI"))
  
  LM <- lm(VI ~ VOI , data=vars)
  r2   <- round(summary(LM)$r.squared,2)
  pval <- round(lmp(LM),4)
  
  gg2 <- ggplot(vars) + theme_bw() + ggtitle(paste0("Nutrient concentrations: r2 = ", r2, ", p = ", pval)) +
    geom_smooth(aes(x=VI, y=VOI), alpha=0.1, formula=y~x, color="gray", method="lm", se = T) +
    geom_point(aes(x=VI, y=VOI, color=PR_FI), size=2, alpha=0.8, na.rm=T) +
    xlim(xminval, xmaxval) +ylim(yminval, ymaxval) +
    xlab("WDRVI") + ylab(VOI) +
    theme(legend.key.size = unit(0.3, 'cm'), plot.title = element_text(size = 12, face = "bold"), 
          legend.title = element_text(size = 10, face = "bold"),
          axis.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 7)) +
    guides(color=guide_legend(ncol=1, bycol=TRUE))
  
  
  grid.arrange(gg2, gg1, ncol=2)
}
dev.off()




print("stop")





ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
  geom_point(aes(x=VI, y=VOI, alpha=0.2, color = TRT), size=2, na.rm=T) +
  geom_smooth(aes(x=VI, y=VOI, color=TRT), formula=y~x, method=lm, se = T) + 
  xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
  # geom_mark_ellipse(aes(x=VI, y=VOI), geom="polygon",level=0.95, alpha=0.1) +
  facet_wrap(~ PR_FI)

# SELECTION - ALL TO PDF
##########################################

VOIs        <- names(vars_all[9:19])
VOIs        <- VOIs[VOIs != "SPAD"]

VIs         <- c("WDRV_cat1", "WDRV_cat2", "WDRV_cat3", "WDRV_cat4", "GWWDRVI", "NWWDRVI")
#VIs         <- c("NDRE_cat1","NDRE_cat2","NDRE_cat3","NDRE_cat4", "NDREGC_cat1", "NDREGC_cat2", "NDREGC_cat3", "NDREGC_cat4")
VI_max      <- max(vars[,VIs])
VI_min      <- min(vars[,VIs])
names(vars_all)
VI <- "WDRV_cat1"
VOI <- "N"
# PDF 1

# if(pdf_writing == "yes"){


pdf(paste0("../Data/3_Output/16_Cor_VOI_VI_def_classes_WDRVI_", SHAPE, ".pdf"), width=15, height=5)
for(VOI in VOIs){
  print(VOI)
  for(VI in VIs){
    vars     <- vars_all
    vars$VI  <- vars[,VI]
    vars$VOI <- vars[,VOI]
    #vars <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
    
    vars     <- omitter(vars, c("VI", "VOI"))
    r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
    
    gg <- ggplot(vars) + theme_bw() + ggtitle(paste(VOI, "-", VI)) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ def_class)
    
    plot(gg)
    
  }
}
dev.off()

VIs         <- c("WDRV_cat1")
pdf(paste0("../Data/3_Output/16_Cor_VOI_VI_def_class_per_nut_WDRVI_", SHAPE, ".pdf"), width=15, height=5)
for(VOI in VOIs){
  print(VOI)
  for(VI in VIs){
    vars     <- vars_all
    vars$VI  <- vars[,VI]
    vars$VOI <- vars[,VOI]
    #vars <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
    
    vars     <- omitter(vars, c("VI", "VOI"))
    r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
    
    defN  <- ggplot(vars) + theme_bw() + ggtitle(paste0(VOI, " - ", VI, " (Deficit class: N)")) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ N_class, ncol=5)
    plot(defN)
    
    defP  <- ggplot(vars) + theme_bw() + ggtitle(paste0(VOI, " - ", VI, " (Deficit class: P)")) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ P_class, ncol=5)
    plot(defP)
    
    defK  <- ggplot(vars) + theme_bw() + ggtitle(paste0(VOI, " - ", VI, " (Deficit class: K)")) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ K_class, ncol=5)
    plot(defK)
    
    defCa  <- ggplot(vars) + theme_bw() + ggtitle(paste0(VOI, " - ", VI, " (Deficit class: Ca)")) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ Ca_class, ncol=5)
    plot(defCa)
    
    defMg  <- ggplot(vars) + theme_bw() + ggtitle(paste0(VOI, " - ", VI, " (Deficit class: Mg)")) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ Mg_class, ncol=5)
    plot(defMg)
    
    defB  <- ggplot(vars) + theme_bw() + ggtitle(paste0(VOI, " - ", VI, " (Deficit class: B)")) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ B_class, ncol=5)
    plot(defB)
    
  }
}
dev.off()

VIs         <- c("WDRV_cat1")
pdf(paste0("../Data/3_Output/16_Cor_VOI_VI_def_lims_WDRVI_", SHAPE, ".pdf"), width=12, height=8)
for(VOI in VOIs){
  print(VOI)
  for(VI in VIs){
    vars     <- vars_all
    vars$VI  <- vars[,VI]
    vars$VOI <- vars[,VOI]
    #vars <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
    
    vars     <- omitter(vars, c("VI", "VOI"))
    r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
    
    gg <- ggplot(transform(vars, nut_lim=factor(nut_lim,levels=c("N_lim", "P_lim", "K_lim", "NP_lim", "NK_lim", "PK_lim", "NPK_lim", "no_lim")))) + 
      theme_bw() + ggtitle(paste(VOI, "-", VI)) +
      geom_smooth(aes(x=VI, y=VOI), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=PR_FI), alpha=0.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) +
      facet_wrap(~ nut_lim)
    plot(gg)
    
    vars_agg    <- aggregator(vars, VI, "nut_lim", "mean")
    vars_agg$sd <- aggregator(vars, VI, "nut_lim", "sd")$WDRV_cat1
    
    ggplot(vars_agg, aes(x=nut_lim, y=WDRV_cat1)) + 
      theme_bw() + ggtitle(paste(VOI, "-", VI)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient2(position="bottom" , low = "blue", mid = "blue", high = "red", midpoint = 200) +
      geom_errorbar(aes(ymin=WDRV_cat1-sd, ymax=WDRV_cat1+sd), width=.2, position=position_dodge(.9)) 
    
    
    
  }
}
dev.off()

VIs         <- c("LA_Palm", "LA_F", "TGF")
VOIs        <- c("N", "P", "K")
pdf(paste0("../Data/3_Output/16_Cor_LAIs_nut_defs.pdf"), width=12, height=8)
for(VOI in VOIs){
  print(VOI)
  for(VI in VIs){
    vars     <- vars_all
    vars$VI  <- vars[,VI]
    vars$VOI <- vars[,VOI]
    #vars <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
    
    vars     <- omitter(vars, c("VI", "VOI"))
    r2       <- round(summary(lm(VI ~ VOI, vars))$r.squared,2)
    pval     <- round(lmp(lm(VI ~ VOI, vars)),5)
    
    gg <- ggplot(transform(vars, nut_lim=factor(nut_lim,levels=c("no_def", "N_def", "P_def", "K_def", "NP_def", "NK_def", "PK_def", "NPK_def")))) + 
      theme_bw() + xlab(VI) + ylab(VOI) + ggtitle(paste0("r2 = ", r2, ", p-value = ", pval)) + 
      scale_color_manual(values=c( "#757575", "#fae105", "#0dd91b", "#09a5ed", "orange", "#6f0091", "#cc0e00", "black")) +
      geom_smooth(aes(x=VI, y=VOI ), formula=y~x, color="black", method="lm", se = T) +
      geom_point(aes(x=VI, y=VOI, color=nut_lim), alpha=0.4, size=2.5, na.rm=T) +
      xlim(min(vars[,"VI"]), max(vars[,"VI"])) + ylim(min(vars[,"VOI"]), max(vars[,"VOI"])) 
    
    plot(gg)
    
  }
}
dev.off()

##########################
## Aggregated per field ##
##########################

VOIs        <- names(vars_all[9:19])
VOIs        <- VOIs[VOIs != "SPAD"]
vars_agg        <- aggregator(vars, VOIs, c("PR_FI"), "mean")
vars_agg_sd     <- aggregator(vars, VOIs, c("PR_FI"), "sd")

vars_agg <- add_class(vars_agg, "N",  "N_class",  c(2.3,   2.4,   2.81,  3.00 ) )
vars_agg <- add_class(vars_agg, "P",  "P_class",  c(0.140, 0.150, 0.181, 0.250) )
vars_agg <- add_class(vars_agg, "K",  "K_class",  c(0.75,  0.90,  1.21,  1.60 ) )

vars_agg$nut_lim <- "no_def"
trh <- 3
vars_agg$nut_lim[vars_agg$N_class > trh & vars_agg$P_class < trh & vars_agg$K_class < trh] <- "N_def"
vars_agg$nut_lim[vars_agg$N_class < trh & vars_agg$P_class > trh & vars_agg$K_class < trh] <- "P_def"
vars_agg$nut_lim[vars_agg$N_class < trh & vars_agg$P_class < trh & vars_agg$K_class > trh] <- "K_def"
vars_agg$nut_lim[vars_agg$N_class > trh & vars_agg$P_class > trh & vars_agg$K_class < trh] <- "NP_def"
vars_agg$nut_lim[vars_agg$N_class > trh & vars_agg$P_class < trh & vars_agg$K_class > trh] <- "NK_def"
vars_agg$nut_lim[vars_agg$N_class < trh & vars_agg$P_class > trh & vars_agg$K_class > trh] <- "PK_def"
vars_agg$nut_lim[vars_agg$N_class > trh & vars_agg$P_class > trh & vars_agg$K_class > trh] <- "NPK_def"

VIs         <- c("LA_Palm", "LA_F", "TGF")
VOIs        <- c("N", "P", "K")



pdf(paste0("../Data/3_Output/16_Cor_LAIs_nut_defs_field_average.pdf"), width=12, height=8)
for(VOI in VOIs){
  print(VOI)
  for(VI in VIs){
    
    vars_agg$VI     <- vars_agg[,VI]
    vars_agg$VOI    <- vars_agg[,VOI]
    vars_agg_sd$VI  <- vars_agg_sd[,VI]
    vars_agg_sd$VOI <- vars_agg_sd[,VOI]
    
    r2       <- round(summary(lm(VI ~ VOI, vars_agg))$r.squared,2)
    pval     <- round(lmp(lm(VI ~ VOI, vars_agg)),5)
    
    gg <- ggplot() + theme_bw() + xlab(VI) + ylab(VOI) + ggtitle(paste0("r2 = ", r2, " p-value = ", pval)) + 
      geom_errorbar(aes(x = vars_agg[,VI], y = vars_agg[,VOI],
                        ymin=vars_agg[,VOI]-vars_agg_sd[,VOI],
                        ymax=vars_agg[,VOI]+vars_agg_sd[,VOI]), color="#777777") +
      geom_errorbarh(aes(x = vars_agg[,VI], y = vars_agg[,VOI],
                         xmin=vars_agg[, VI]-vars_agg_sd[, VI],
                         xmax=vars_agg[, VI]+vars_agg_sd[, VI]), color="#777777") +
      geom_smooth(data=vars_agg, aes(x=VI, y=VOI), method="lm", color="black") +
      geom_point(data=transform(vars_agg, nut_lim=factor(nut_lim,levels=c("no_def", "N_def", "P_def", "K_def", "NP_def", "NK_def", "PK_def", "NPK_def"))), 
                 mapping=aes(x=VI, y=VOI, color=nut_lim), size=3, alpha=0.7) +
      scale_color_manual(values=c( "#888888", "#09a5ed", "#6f0091", "#cc0e00", "black"))
    
    plot(gg)
    
  }
}
dev.off()






















































