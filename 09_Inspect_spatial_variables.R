# Bastiaen Boekelo, April 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Inspect correlations between created spatial variables and ground measurement


rm(list=ls())  # Clean script <- clean mind

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
  "../Data/3_Output/08_spatial_variables_shape_1-0.5.csv"
)

# Read input data

df_vi <- read.csv(vi_files[1], stringsAsFactors = F)
df_nu <- read.csv(paste0(wd, "1_Input/Onground/sheet1_Intensive sampling 210310__intensive_sampling.csv"),      stringsAsFactors=F)
df_sp <- read.csv(paste0(wd, "1_Input/Onground/sheet2_Intensive sampling 210310__SPAD_Measurement_Year_2.csv"), stringsAsFactors=F)

# Preprocess data before merge
#df_vi$ID       <- paste0(substr(df_vi$ID,1,10), "P", substr(df_vi$ID,11,12))

df_nu$DOM_nut  <- df_nu$DOS
df_nu$FI       <- substr(unlist(map(strsplit(df_nu$Field_ID, split = "-"), 2)),1,2)
df_nu          <- add_zero(df_nu, "Sample_point") # Creates column "NR"
df_nu$NR_ID    <- paste0("P", df_nu$NR_ID)
df_nu$NR_nut   <- df_nu$No
df_nu$ID       <- paste(df_nu$Site, df_nu$FI, df_nu$Treatment, df_nu$NR_ID, sep="_")
#df_nu$ID       <- str_replace(df_nu$ID, "_P", "_")
df_nu[,c("Site", "Farmer_ID", "DOS", "Lab_Ref", "Treatment", "PR",
         "NR_ID", "S_code", "No", "Field_ID", "Sample_point")] <- NULL

df_sp$DOM_spad <- df_sp$DOM
df_sp$FI       <- substr( unlist(map(strsplit(df_sp$Field_ID, split = "-"), 2)), 1, 2)
df_sp          <- add_zero(df_sp, "Sample_point") # Creates column "NR"
df_sp$NR_ID    <- paste0("P", df_sp$NR_ID) 
df_sp$ID       <- paste(df_sp$Site, df_sp$FI, df_sp$Treatment, df_sp$NR_ID, sep="_")
#df_sp$ID       <- str_replace(df_sp$ID, "_P", "_")
df_sp$SPAD     <- rowMeans(df_sp[,c('SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')])
df_sp[,c("Site", "Field_ID", "Farmer", "Treatment", "DOM", "PR", "Sample_point", "FI", "NR_ID",
         'SPAD1', 'SPAD2', 'SPAD3', 'SPAD4', 'SPAD5', 'SPAD6')] <- NULL

# Combine all data
m1       <- merge(df_sp, df_nu, by="ID", all.x=T, all.y=T)
m1[,c("PR", "FI")] <- NULL
m2       <- merge(m1, df_vi,    by="ID", all.x=F, all.y=F)
m2$PR_FI <- paste0(m2$PR, "_", m2$FI)
# md       <- setdiff(m3,m2)
vars       <- m2
vars$NR    <- paste0(substr(vars$TRT, 1,1), substr(vars$NR,2,3))
vars$NR    <- paste0(vars$PR_FI, "_", vars$NR)
vars       <- sort_df(vars, c("ID", "PR", "FI", "PR_FI", "TRT", "Farmer", "NR", "NR_nut", "DOM_nut", "DOM_spad"))

unique(vars$PR_FI)

test <- vars[vars$PR_FI == "CK_F1" | vars$PR_FI == "CK_F7",]
ggplot(test, aes(NDVI_cat4, N)) + 
  geom_point(aes(colour=TRT), size=2) + 
  geom_smooth(method = "lm", se=F) + 
  labs(x="NDRE", y="N concentration (% DM)",title="CK-F1&F7")

lm_test <- lm()

# # Create correlation plots
# var_cor        <- vars
# var_cor        <- var_cor[complete.cases(var_cor), ]      # Remove rows containing NA
# names(var_cor) <- paste0("v_", as.character(1:ncol(var_cor)))
# corrplot(cor(var_cor[10:ncol(var_cor)]), method="circle")
# corrplot(cor(var_cor[30:ncol(var_cor)]), method="circle")
# corrplot(cor(var_cor[10:44]),            method="circle") # VI1 (NDRE)
# corrplot(cor(var_cor[c(10:20, 69:92)]),  method="circle") # VI2 (NDVI)
# corrplot(cor(var_cor[c(10:20, 93:116)]), method="circle") # VI2 (NDVI)
# 

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
VOIs    <- c("SPAD", "Ash", "N", "P", "K", "Mg", "Ca", "B", "LA_F", "TGF", "LA_Palm")
VIs     <- names(vars)[21:length(names(vars))]
nr_rows <- length(FIs) * length(VOIs) * length(VIs)

stats[1:nr_rows,] <- NA

FI  <- "CK_F1"
VI  <- "NDRE_cat4"
VOI <- "SPAD"

i <- 0
for(FI in FIs){
  for(VI in VIs){
    for(VOI in VOIs){
      
      vars_FI     <- vars[vars$PR_FI == FI, ]
      vars_FI$VI  <- vars_FI[, VI]
      vars_FI$VOI <- vars_FI[, VOI]
      vars_FI     <- vars_FI[, c("ID", "TRT", "NR", "VOI", "VI")]

      stats$PR_FI[i] <- FI
      stats$VI[i]    <- VI
      stats$VOI[i]   <- VOI
      
      if(sum(!is.na(vars_FI[,"VOI"])) > 36 & sum(!is.na(vars_FI[,"VI"])) > 36){
        i <- i + 1  # Row number to fill
        
        vars_FI <- vars_FI[complete.cases(vars_FI), ]
        vars_lm <- lm(VI ~ VOI, vars_FI, na.action=na.omit)
        
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

# write.csv(stats,     paste0(wd, "2_Intermediate/05_statistics_per_field.csv"),     row.names = F)
# write.csv(stats_agg, paste0(wd, "2_Intermediate/05_statistics_per_parameter.csv"), row.names = F)


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

