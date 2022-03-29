# Bastiaen Boekelo, November 2021
# Nebraska project - UAV and Oil Palm Nutrients
# Goal: Create correlation plots


############################
#### SET-UP ENVIRONMENT ####
############################

rm(list=ls()) 

# Set script and data wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Works if you have a recent version of RStudio.

# Functions
library(ggcorrplot)
library(maditr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(ggrepel) # prevents overlapping point labels
library(viridis)
source("functions/aggregator.R")
source("functions/sort_df.R")

###########################
####     READ DATA     ####
###########################

df_ind     <- read.csv("../Data/2_Intermediate/17_Preprocessed_palms.csv",     stringsAsFactors = F)
df_fi      <- read.csv("../Data/2_Intermediate/17_Preprocessed_PR_FI.csv",     stringsAsFactors = F)
df_fi_trt  <- read.csv("../Data/2_Intermediate/17_Preprocessed_PR_FI_TRT.csv", stringsAsFactors = F)

# Rename variable
df_ind$Fr_Yr               <- df_ind$FrProd_FrYr
df_ind$FrProd_FrYr         <- NULL
df_fi$Fr_Yr                <- df_fi$FrProd_FrYr
df_fi$FrProd_FrYr          <- NULL
df_fi_trt$Fr_Yr            <- df_fi_trt$FrProd_FrYr
df_fi_trt$FrProd_FrYr      <- NULL

# Rename variable
df_ind$DW_prod             <- df_ind$Pa_Tot_kg_prod
df_ind$Pa_Tot_kg_prod      <- NULL
df_fi$DW_prod              <- df_fi$Pa_Tot_kg_prod
df_fi$Pa_Tot_kg_prod       <- NULL
df_fi_trt$DW_prod          <- df_fi_trt$Pa_Tot_kg_prod
df_fi_trt$Pa_Tot_kg_prod   <- NULL


############################################################
####    CALCULATE IN_FIELD CORRELATIONS AND P-VALUES    ####
############################################################


r_calc_ind <- function(data, VOIs, VIs, fields){
  
  new_df   <- data.frame(
    PR_FI=character(),
    var1=character(),
    var2=character(),
    r=numeric()
  )
  
  i <- 0
  for(PR_FI in fields){
    for(VAR1 in VIs){
      for(VAR2 in VOIs){
        
        df   <- data[data$PR_FI == PR_FI,]
        
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
  new_df$r         <- as.numeric(new_df$r)
  
  r_avg            <- aggregator(new_df, "r", c("var1", "var2"), "mean")
  r_avg            <- as.data.frame(dcast(r_avg, var1 ~ var2, value.var = "r"))
  r_avg            <- r_avg %>% slice(match(VIs, var1))
  
  row.names(r_avg) <- r_avg$var1
  r_avg$var1       <- NULL
  
  r_avg            <- sort_df(r_avg, VOIs)
  
  return(r_avg)
}

# Calculate p-values
p_calc_ind <- function(data, VOIs, VIs){

  new_df <- data.frame(
    var1=character(),
    var2=character(),
    p=numeric()
  )
  
  i <- 0
  for(VAR1 in VIs){
    for(VAR2 in VOIs){
      i <- i + 1
      
      lm       <- lm(data[,VAR2] ~ data[,VAR1] + data[,"PR_FI"])
      p_coef   <- summary(lm)$coefficients[2,4]
      
      new_df[i,] <- c(VAR1, VAR2, p_coef)
      
    }
  }
  
  p_val            <- as.data.frame(dcast(new_df, var1 ~ var2, value.var = "p"))
  p_val            <- p_val %>% slice(match(VIs, var1))
  
  rowies           <- p_val$var1
  rownames(p_val)  <- rowies
  p_val$var1       <- NULL
  p_val            <- sort_df(p_val, VOIs)
  p_val            <- sapply(p_val, as.numeric)
  rownames(p_val)  <- rowies
  
  return(p_val)
  
}

# Calculate correlation plots
VOIs_inf_pc    <- c("N", "P", "K", "Mg", "Ca", "B")
VOIs_inf_kg    <- c("N_Tot_Fr_kg",  "P_Tot_Fr_kg",  "K_Tot_Fr_kg",  "Mg_Tot_Fr_kg",  "Ca_Tot_Fr_kg",  "Fr_Tot_kg") # Leaflets
# VOIs_inf_kg    <- c("N_Fr_Tot_kg", "P_Fr_Tot_kg", "K_Fr_Tot_kg", "Mg_Fr_Tot_kg", "Ca_Fr_Tot_kg", "Fr_Tot_kg") # Rachis + Petiole + Leaflets
VOIs_inf_prod  <- c("N_Pa_Tot_kg_prod",  "P_Pa_Tot_kg_prod",  "K_Pa_Tot_kg_prod",  "Mg_Pa_Tot_kg_prod",  "Ca_Pa_Tot_kg_prod", "DW_prod")
VOIs_inf_pf    <- c("RL", "LA_F", "DW_prod")

#VIs0     <- c("NDVI", "NDRE", "WDRVI", "MSAVI", "NRED", "NGRN", "NNIR")
VIs0     <- c("NDVI", "NDRE", "WDRVI", "MSAVI")
VIs1     <- c("WDRVI", "NWWDRVI", "GWWDRVI", "WDRVI_LNSY", "WDRVI_LYSN", "WDRVI_LNSN")
VIs2     <- c("NDRE",  "NWNDRE",   "GWNDRE",  "NDRE_LNSY",  "NDRE_LYSN",  "NDRE_LNSN")
VIs3     <- names(df_ind)[grepl("WDRVI__", names(df_ind))]
VIs4     <- names(df_ind)[grepl("NDRE__", names(df_ind))]

PR_FIs   <- unique(df_ind$PR_FI)

create_corplot_PR_FI <- function(df_ind, VOIs, VIs, PR_FIs, TYPE="full"){
  
  r_avg    <- r_calc_ind(df_ind, VOIs, VIs, PR_FIs)
  p_val    <- p_calc_ind(df_ind, VOIs, VIs)
  
  gg <- ggcorrplot(r_avg,
             hc.order = F,
             lab = T,
             type=TYPE,
             lab_size = 3,
             outline.color = "#777777",
             colors = c("blue", "white", "red"),
             tl.cex = 9,
             tl.srt = 90,
             p.mat=p_val,
             sig.level=0.05,
             pch = 19, #or 21
             pch.cex = 9,
             pch.col = "#999999",
             legend.title = "Pearson's r"
  )
  plot(gg)
  return(gg)
}

#pdf("../Data/3_Output/18_Correlation_plots.pdf")
VIs_all_inf   <- create_corplot_PR_FI(df_ind, c(VOIs_inf_pc, VOIs_inf_kg, VOIs_inf_pf), VIs0, PR_FIs)
# VIs_all_inf_kg   <- create_corplot_PR_FI(df_ind, c(VOIs_inf_kg, VOIs_inf_pf), VIs0, PR_FIs)
# VIs_all_inf_prod <- create_corplot_PR_FI(df_ind, VOIs_inf_prod, VIs0, PR_FIs)

WDRVI_all     <- create_corplot_PR_FI(df_ind, c(VOIs_inf_pc, VOIs_inf_kg, VOIs_inf_pf), c(VIs1, VIs3), PR_FIs)
NDRE_all      <- create_corplot_PR_FI(df_ind, c(VOIs_inf_pc, VOIs_inf_kg, VOIs_inf_pf), c(VIs2, VIs4), PR_FIs)

#dev.off()

####################################################
#### EXAMPLES BASED ON RESULTS - CONCENTRATIONS  ###
####################################################

sel <- df_ind[is.na(df_ind$N) == F,]

# KEY FIGURE 1 - Concentrations
ggplot(sel) + 
  geom_smooth(aes(x=NDRE/10000, y=N ), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=NDRE/10000, y=N, color = as.character(def_tot)), alpha=0.8, size=1.75, na.rm=T) +
  labs(color='Limited (#)') + ylab("N (%)") + xlab("NDRE") +
  scale_color_viridis(direction = -1, begin=0.08, discrete = T, end=1) +
  facet_wrap(~PR_FI) + 
  theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))

# KEY FIGURE 1 - Amounts
ggplot(sel) + 
  geom_smooth(aes(x=NDRE/10000, y=N_Tot_Fr_kg ), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=NDRE/10000, y=N_Tot_Fr_kg, color = as.character(def_tot)), alpha=0.8, size=1.75, na.rm=T) +
  labs(color='Limited (#)') + ylab("N (kg / 17th frond)") + xlab("NDRE") +
  scale_color_viridis(direction = -1, begin=0.08, discrete = T, end=1) +
  facet_wrap(~PR_FI) + 
  theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))

# # Other figures
# ggplot(sel) +
#   geom_smooth(aes(x=NDRE/10000, y=N ), formula=y~x, color="black", method="lm", se = T) +
#   geom_point(aes(x=NDRE/10000, y=N, color = TRT), alpha=0.5, size=2, na.rm=T) +
#   scale_color_manual(values=c("#0088FF", "#FF9922")) +
#   labs(color='Treatment') + ylab("N (%)") + xlab("NDRE") +
#   facet_wrap(~PR_FI) +
#   theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))
# 
# ggplot(sel) +
#   geom_smooth(aes(x=NDRE/10000, y=N ), formula=y~x, color="black", method="lm", se = T) +
#   geom_point(aes(x=NDRE/10000, y=N, shape=TRT, color = as.character(def_tot)), alpha=0.8, size=1.75, na.rm=T) +
#   labs(color='Limited (#)', shape="Treatment") + ylab("N (%)") + xlab("NDRE") +
#   scale_color_viridis(direction = -1, begin=0.08, discrete = T, end=1) +
#   facet_wrap(~PR_FI) +
#   theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))
# 
# ggplot(sel) +
#   geom_smooth(aes(x=NDRE/10000, y=P ), formula=y~x, color="black", method="lm", se = T) +
#   geom_point(aes(x=NDRE/10000, y=P, shape=TRT, color = as.character(def_tot)), alpha=0.8, size=1.75, na.rm=T) +
#   labs(color='Limited (#)', shape="Treatment") + ylab("P (%)") + xlab("NDRE") +
#   scale_color_viridis(direction = -1, begin=0.08, discrete = T, end=1) +
#   facet_wrap(~PR_FI) +
#   theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))
# 
# ggplot(sel) +
#   geom_smooth(aes(x=NDRE/10000, y=K ), formula=y~x, color="black", method="lm", se = T) +
#   geom_point(aes(x=NDRE/10000, y=K, shape=TRT, color = as.character(def_tot)), alpha=0.8, size=1.75, na.rm=T) +
#   labs(color='Limited (#)', shape="Treatment") + ylab("K (%)") + xlab("NDRE") +
#   scale_color_viridis(direction = -1, begin=0.08, discrete = T, end=1) +
#   facet_wrap(~PR_FI) +
#   theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))
# 
# ggplot(sel) +
#   geom_smooth(aes(x=NDRE/10000, y=LA_F ), formula=y~x, color="black", method="lm", se = T) +
#   geom_point(aes(x=NDRE/10000, y=LA_F, shape=TRT, color = as.character(def_tot)), alpha=0.8, size=2, na.rm=T) +
#   labs(color='Limited (#)', shape="Treatment") + ylab(bquote("Leaf area / frond  ("~m^2~")")) + xlab("NDRE") +
#   scale_color_viridis(direction = -1, begin=0.08, discrete = T, end=1) +
#   facet_wrap(~PR_FI) + #, nrow = 4,  ncol = 3, scales = "fixed") +
#   theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))

ggplot(sel) + 
  geom_smooth(aes(x=def_tot, y=N ), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=def_tot, y=N, color = NDRE/10000), alpha=0.8, size=1.75, na.rm=T) +
  labs(color='NDRE') + ylab("N (%)") + xlab("Limited (#)") +
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  facet_wrap(~PR_FI) + 
  theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))


############################################################
#### CALCULATE ACROSS FIELD CORRELATIONS AND P-VALUES   ####
############################################################

# Calculate correlation plots
VOIs_afi_pc  <- c("N_reg_Lf", "P_reg_Lf", "K_reg_Lf", "Mg_reg_Lf", "Ca_reg_Lf")


VIs  <- VIs1
data <- df_fi
p_r_calc_fi <- function(data, VOIs, VIs){
  
  df <- data
  
  df_p <- data.frame(
    var1=character(),
    var2=character(),
    p=numeric()
  )
  
  df_r <- data.frame(
    var1=character(),
    var2=character(),
    r=numeric()
  )

  i <- 0
  for(VAR1 in VIs){
    for(VAR2 in VOIs){
      i <- i + 1

      lm       <- lm(df[,VAR2] ~ df[,VAR1])
      p_coef   <- summary(lm)$coefficients[2,4]
      r        <- sqrt(summary(lm)$r.squared) * as.numeric(lm$coefficients[2]) / abs(lm$coefficients[2]) # direction negative or positive
      
      df_p[i,] <- c(VAR1, VAR2, p_coef)
      df_r[i,] <- c(VAR1, VAR2, r)
      
    }
  }
  
  # Reorder r values
  df_r$r           <- as.numeric(as.character(df_r$r))
  r_all            <- as.data.frame(dcast(df_r, var1 ~ var2, value.var = "r"))
  r_all            <- r_all %>% slice(match(VIs, var1))
  row.names(r_all) <- r_all$var1
  r_all$var1       <- NULL
  r_all            <- sort_df(r_all, VOIs)

  # Reorder p values
  p_val            <- df_p
  p_val$p          <- as.numeric(as.character(p_val$p))
  p_val            <- as.data.frame(dcast(df_p, var1 ~ var2, value.var = "p"))
  p_val            <- p_val %>% slice(match(VIs, var1))
  
  rowies           <- p_val$var1
  rownames(p_val)  <- rowies
  p_val$var1       <- NULL
  
  p_val            <- sort_df(p_val, VOIs)
  p_val            <- sapply(p_val, as.numeric)
  rownames(p_val)  <- rowies
  
  return(list(r_all, p_val))
  
}

# Correlation plot

create_corplot_across <- function(df, VOIs, VIs, TYPE="full"){
  
  r_all <- p_r_calc_fi(df, VOIs, VIs)[[1]] # When CK_F1 is gone, all is insignificamt
  p_all <- p_r_calc_fi(df, VOIs, VIs)[[2]] # When CK_F1 is gone, all is insignificant
  
  colnames(r_all) <- gsub("_reg_Lf", " (%)", colnames(r_all))
  colnames(p_all) <- gsub("_reg_Lf", " (%)", colnames(p_all))
  rownames(r_all) <- gsub("_reg_Lf", " (%)", rownames(r_all))
  rownames(p_all) <- gsub("_reg_Lf", " (%)", rownames(p_all))
  
  colnames(r_all) <- gsub("_Tot_Fr_kg", " (kg)", colnames(r_all))
  colnames(p_all) <- gsub("_Tot_Fr_kg", " (kg)", colnames(p_all))
  rownames(r_all) <- gsub("_Tot_Fr_kg", " (kg)", rownames(r_all))
  rownames(p_all) <- gsub("_Tot_Fr_kg", " (kg)", rownames(p_all))
  
  
  test <- p_all
  test[test <   0.05]   <- 3
  test[test < 0.1]      <- 2
  test <- as.vector(as.matrix(test))
  
  gg <- ggcorrplot(r_all,
             hc.order = F,
             lab = T,
             type = TYPE,
             #lab_size = test,
             lab_size = 3,
             outline.color = "#777777",
             colors = c("blue", "white", "red"),
             tl.cex = 9,
             tl.srt = 90,
             p.mat=p_all,
             sig.level=0.05,
             pch = 19, #or 21
             pch.cex = 9,
             pch.col = "#999999",
             legend.title = "Pearson's r"
  )
  plot(gg)
  return(gg)
}


VIs_all_acr      <- create_corplot_across(df_fi, c(VOIs_afi_pc, VOIs_inf_kg, VOIs_inf_pf), VIs0)

# pix_man_WDRVI <- create_corplot_across(df_fi, VOIs_afi, VIs1)
# pix_man_NDRE  <- create_corplot_across(df_fi, VOIs_afi, VIs2)
# pix_man_both  <- create_corplot_across(df_fi, VOIs_afi, c(VIs1, VIs2))
# shapes_WDRVI  <- create_corplot_across(df_fi, VOIs_afi, VIs3)
# shapes_NDRE   <- create_corplot_across(df_fi, VOIs_afi, VIs4)
# shapes_both   <- create_corplot_across(df_fi, VOIs_afi, c(VIs3, VIs4))
WDRVI_all        <- create_corplot_across(df_fi, c(VOIs_afi_pc, VOIs_inf_kg, VOIs_inf_pf), c(VIs1, VIs3))
NDRE_all         <- create_corplot_across(df_fi, c(VOIs_afi_pc, VOIs_inf_kg, VOIs_inf_pf), c(VIs2, VIs4))

# # Without CK_F1
# WDRVI_all_19     <- create_corplot_across(df_fi[df_fi$PR_FI != "CK_F1",], c(VOIs_afi_pc, VOIs_afi_kg), c(VIs1, VIs3))
# NDRE_all_19      <- create_corplot_across(df_fi[df_fi$PR_FI != "CK_F1",], c(VOIs_afi_pc, VOIs_afi_kg), c(VIs2, VIs4))


create_corplot_across(df_fi,
   c("N_reg_Lf",  "P_reg_Lf", "K_reg_Lf", "N_Tot_Fr_kg", "P_Tot_Fr_kg", "K_Tot_Fr_kg", "RL", "LA_F", "DW_prod"),
   c("N_reg_Lf",  "P_reg_Lf", "K_reg_Lf", "N_Tot_Fr_kg", "P_Tot_Fr_kg", "K_Tot_Fr_kg", "RL", "LA_F", "DW_prod"),
   TYPE="lower")

############################################################
#### EXAMPLES BASED ON RESULTS   ####
############################################################


# WDRVI ~ Leaf area / frond
ggplot(df_fi) + 
  geom_smooth(aes(x=WDRVI/1000, y=LA_F), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=WDRVI/1000, y=LA_F, color=P_reg_Lf), alpha=0.9, size=2, na.rm=T) +
  labs(color="P concentration (%)") + 
  ylab(bquote("Leaf area / frond  ("~m^2~")")) + xlab("WDRVI") +
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  geom_text(aes(x=WDRVI/1000, y=LA_F, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 
# LA_F correlates with WDRVI (and higher LA_F seems to have higher P concentrations)

# WDRVI ~ P (%)
ggplot(df_fi) + 
  geom_smooth(aes(x=WDRVI/1000, y=P_reg_Lf), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=WDRVI/1000, y=P_reg_Lf, color=LA_F), alpha=0.9, size=2, na.rm=T) +
  ylab("P concentration (%)") + xlab("WDRVI") +
  labs(color=bquote("Leaf area / frond  ("~m^2~")")) + 
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  geom_text(aes(x=WDRVI/1000, y=P_reg_Lf, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 
# P concentration also correlates with WDRVI (and seems to relate with the size of leafs)

# WDRVI ~ P (%)
ggplot(df_fi) + 
  geom_smooth(aes(x=WDRVI/1000, y=P_reg_Lf), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=WDRVI/1000, y=P_reg_Lf, color=Lf_prod), alpha=0.9, size=2, na.rm=T) +
  ylab("P concentration (%)") + xlab("WDRVI") +
  labs(color="Lf prod. (kg / year)") + 
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  geom_text(aes(x=WDRVI/1000, y=P_reg_Lf, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 
# P concentration also correlates with WDRVI (and seems to relate with leaf production as well)

# Leaf area / frond ~ P (%)   [ZONDER CK_F1]
ggplot(df_fi[df_fi$PR_FI != "CK_F1",]) + 
  geom_smooth(aes(x=LA_F, y=P_reg_Lf), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=LA_F, y=P_reg_Lf, color=WDRVI), alpha=0.9, size=2, na.rm=T) +
  ylab("P concentration (%)") + xlab(bquote("Leaf area / frond  ("~m^2~")")) +
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  geom_text(aes(x=LA_F, y=P_reg_Lf, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 
# Size of leafs have a relation with the P concentration

# Leaflet production ~ P (%)
ggplot(df_fi) + 
  geom_smooth(aes(x=P_reg_Lf, y=Lf_prod), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=P_reg_Lf, y=Lf_prod, color=WDRVI), alpha=0.9, size=2, na.rm=T) +
  ylab("Leaflet production (kg / year)") + xlab("P concentration (%)") + 
  xlim(  c( min(df_fi$P_reg_Lf), max(df_fi$P_reg_Lf) )  ) +
  ylim(  c( min(df_fi$Lf_prod) - 1,  max(df_fi$Lf_prod)  )  ) +
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  geom_text(aes(x=P_reg_Lf, y=Lf_prod, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 

# Frond production ~ P (%)
ggplot(df_fi) + 
  geom_smooth(aes(x=P_reg_Lf, y=FrProd_FrYr), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=P_reg_Lf, y=FrProd_FrYr, color=WDRVI), alpha=0.9, size=3, na.rm=T) +
  ylab("Frond production (# / year)") + xlab("P concentration (%)") + 
  xlim(  c( min(df_fi$P_reg_Lf), max(df_fi$P_reg_Lf) )  ) +
  ylim(  c( min(df_fi$FrProd_FrYr) - 1,  max(df_fi$FrProd_FrYr)  )  ) +
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  geom_text(aes(x=P_reg_Lf, y=FrProd_FrYr, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 


# Leaflet production ~ P (%)   [ZONDER CK_F1]
ggplot(df_fi[df_fi$PR_FI != "CK_F1",]) + 
  geom_smooth(aes(x=P_reg_Lf, y=Lf_prod), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=P_reg_Lf, y=Lf_prod, color=WDRVI), alpha=0.9, size=2, na.rm=T) +
  ylab("Leaflet production (kg / year)") + xlab("P concentration (%)") + 
  xlim(  c( min(df_fi$P_reg_Lf), max(df_fi$P_reg_Lf) )  ) +
  ylim(  c( min(df_fi$Lf_prod) - 1,  max(df_fi$Lf_prod)  )  ) +
  scale_color_viridis(direction = 1, begin=0.08, discrete = F, end=1) +
  geom_text(aes(x=P_reg_Lf, y=Lf_prod, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 

summary(lm(Lf_prod ~ P_reg_Lf , df_fi[df_fi$PR_FI != "CK_F1",]))


# Leaflet production ~ P (%)   [ZONDER CK_F1]
ggplot(df_fi) + 
  geom_smooth(aes(x=P_reg_Lf, y=Lf_prod), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=P_reg_Lf, y=Lf_prod, color=as.character(def_reg_tot), size=WDRVI), alpha=0.9, na.rm=T) +
  ylab("Leaflet production (kg / year)") + xlab("P concentration (%)") + 
  xlim(  c( min(df_fi$P_reg_Lf), max(df_fi$P_reg_Lf) )  ) +
  ylim(  c( min(df_fi$Lf_prod) - 1,  max(df_fi$Lf_prod)  )  ) +
  scale_color_viridis(direction = -1, begin=0, end=0.98,discrete = T, ) +
  geom_text(aes(x=P_reg_Lf, y=Lf_prod, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 


# KEY FIGURE 1 - WDRVI and P relation seems to be influenced by number of limiting nutrients (Related to Fronds / year)
ggplot(df_fi) + 
  geom_smooth(aes(x=WDRVI/1000, y=P_reg_Lf), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=WDRVI/1000, y=P_reg_Lf, color=as.character(def_reg_tot), size=Fr_Yr), alpha=0.9, na.rm=T) +
  geom_text_repel(aes(x=WDRVI/1000, y=P_reg_Lf, label=PR_FI), direction = "y" , hjust=1, cex=3, alpha = 0.5, force = 0.1, force_pull = 2) +  
  ylab("P concentration (%)") + xlab("WDRVI") + 
  labs(size="# Fronds / year", color="# Limited") + 
 # xlim(  c( 0, 0.5 ) ) + ylim(c(0.13, 0.18)) +
  scale_color_viridis(direction = -1, begin=0, end=0.98,discrete = T, ) +
  theme_bw() 

summary(lm(P_reg_Lf ~ WDRVI, df_fi))
summary(lm(P_reg_Lf ~ WDRVI + def_reg_tot, df_fi))
summary(lm(WDRVI    ~ P_reg_Lf + log(def_reg_tot), df_fi))   # Everything significant
summary(lm(P_reg_Lf ~ WDRVI, df_fi[df_fi$def_reg_tot > 2,]))

ggplot(df_fi) + 
  geom_smooth(aes(x=WDRVI/1000, y=P_reg_Lf), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=WDRVI/1000, y=P_reg_Lf, color=as.character(def_reg_tot), size=P_Pa_Tot_kg_prod), alpha=0.9, na.rm=T) +
  geom_text_repel(aes(x=WDRVI/1000, y=P_reg_Lf, label=PR_FI), direction = "y" , hjust=1, cex=3, alpha = 0.5, force = 0.1, force_pull = 2) +  
  ylab("P concentration (%)") + xlab("WDRVI") + 
  labs(size="kg DW / year", color="# Limited") + 
  # xlim(  c( 0, 0.5 ) ) + ylim(c(0.13, 0.18)) +
  scale_color_viridis(direction = -1, begin=0, end=0.98,discrete = T, ) +
  theme_bw() 


# KEY FIGURE 2 - Indeed, it's all connected (And reflected in Fronds / year)
ggplot(df_fi) + 
  geom_smooth(aes(x=N_reg_Lf, y=P_reg_Lf), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=N_reg_Lf, y=P_reg_Lf, color=K_reg_Lf, size=Fr_Yr), alpha=0.9, na.rm=T) +
  geom_text_repel(aes(x=N_reg_Lf, y=P_reg_Lf, label=PR_FI), direction = "y" , hjust=1, cex=3, alpha = 0.5, force = 0.1, force_pull = 2) +  
  ylab("P concentration (%)") + xlab("N concentration (%)") + 
  labs(size="# Fronds / year", color="K concentration (%)") + 
  #xlim(  c( 0, 0.5 ) ) + ylim(c(0.13, 0.18)) +
  scale_color_viridis(direction = 1, begin=0, end=1,discrete = F, ) +
  theme_bw() 

summary(lm(df_fi$N_reg_Lf ~ df_fi$P_reg_Lf))




##################################################################
##### SHOW VARIATION CAUSES DIFFERENCE IN QUALITY OF PREDICTION #######
##################################################################

# Same function, but now averaging r's over PR_FI per vegetation index
r_calc_FIELD <- function(data, VOIs, VIs, fields){
  
  new_df   <- data.frame(
    PR_FI=character(),
    var1=character(),
    var2=character(),
    r=numeric()
  )
  
  i <- 0
  for(PR_FI in fields){
    for(VAR1 in VIs){
      for(VAR2 in VOIs){
        
        df   <- data[data$PR_FI == PR_FI,]
        
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
  new_df$r         <- as.numeric(new_df$r)
  

  return(new_df)
}

# Add in-field nutrient deficiency variation - concentrations
######################################################################

# Calculate variance per field
temp             <- aggregator(df_ind, "def_tot", c("PR_FI"), "sd")
temp             <- na.omit(temp)
temp$var_def_tot <- temp$def_tot^2
temp$def_tot     <- NULL
temp$Function    <- NULL

variables        <- c("N", "P", "K", "Mg", "Ca", "B", "RL", "LA_F", "DW_prod")

# Calculate Pearson's r for every combination
test1            <- r_calc_FIELD(df_ind, variables, VIs0, unique(temp$PR_FI))
test2            <- test1[test1$var1=="NDRE",]
test3            <- merge(test2, temp, by="PR_FI")
test4            <- transform(test3, var2=factor(var2, levels=variables ) )
rm(temp, test1, test2, test3)

ggplot(test4) + 
  geom_smooth(aes(x=var_def_tot, y=r ), formula=y~x, color="black", method="lm", se = T) +
  geom_point(aes(x=var_def_tot, y=r), alpha=0.8, size=1.75, na.rm=T) +
  labs(color='Limited (#)', shape="Treatment") + ylab("Pearon's r") + xlab("# Limited nutrients - Variance") +
  facet_wrap(~var2) +
  geom_text_repel(aes(x=var_def_tot, y=r, label=PR_FI), direction = "y" , hjust=1, cex=3, alpha = 0.5, force = 0.1, force_pull = 2) +  
  theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))

summary(lm(r ~ var_def_tot, test4[test4$var2 == "N",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "P",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "K",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "Mg",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "Ca",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "B",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "RL",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "LA_F",]))$r.squared
summary(lm(r ~ var_def_tot, test4[test4$var2 == "DW_prod",]))$r.squared

summary(lm(r ~ var_def_tot, test4[test4$var2 == "N",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "P",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "K",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "Mg",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "Ca",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "B",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "RL",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "LA_F",]))$coefficients[2,4]
summary(lm(r ~ var_def_tot, test4[test4$var2 == "DW_prod",]))$coefficients[2,4]


# # Relation variation nutrient deficiency and correlation quality
# 
# variables        <- c("N_Tot_Fr_kg", "P_Tot_Fr_kg", "K_Tot_Fr_kg", "Mg_Tot_Fr_kg", "Ca_Tot_Fr_kg", "RL", "LA_F", "Pa_Tot_kg_prod")
# 
# # Calculate Pearson's r for every combination
# test1            <- r_calc_FIELD(df_ind, variables, VIs0, unique(temp$PR_FI))
# test2            <- test1[test1$var1=="NDRE",]
# test3            <- merge(test2, temp, by="PR_FI")
# test4            <- transform(test3, var2=factor(var2, levels=variables ) )
# rm(temp, test1, test2, test3)
# 
# # Not useful
# ggplot(test4) + 
#   geom_smooth(aes(x=var_def_tot, y=r ), formula=y~x, color="black", method="lm", se = T) +
#   geom_point(aes(x=var_def_tot, y=r), alpha=0.8, size=1.75, na.rm=T) +
#   labs(color='Limited (#)', shape="Treatment") + ylab("Pearon's r") + xlab("# Limited nutrients - Variance") +
#   facet_wrap(~var2) +
#   geom_text_repel(aes(x=var_def_tot, y=r, label=PR_FI), direction = "y" , hjust=1, cex=3, alpha = 0.5, force = 0.1, force_pull = 2) +  
#   theme_bw() + theme(strip.background = element_rect(fill="#FFFFDD"))




# Nutrient limitations table
######################################################################

df_ind$nut_class <- NA
df_ind$nut_class[df_ind$N_class == 1 & df_ind$P_class == 1 & df_ind$K_class == 1] <- "NPK_lim"
df_ind$nut_class[df_ind$N_class == 1 & df_ind$P_class == 1 & df_ind$K_class == 0] <- "NP_lim"
df_ind$nut_class[df_ind$N_class == 1 & df_ind$P_class == 0 & df_ind$K_class == 1] <- "NK_lim"
df_ind$nut_class[df_ind$N_class == 0 & df_ind$P_class == 1 & df_ind$K_class == 1] <- "PK_lim"
df_ind$nut_class[df_ind$N_class == 0 & df_ind$P_class == 0 & df_ind$K_class == 1] <- "K_lim"
df_ind$nut_class[df_ind$N_class == 0 & df_ind$P_class == 1 & df_ind$K_class == 0] <- "P_lim"
df_ind$nut_class[df_ind$N_class == 1 & df_ind$P_class == 0 & df_ind$K_class == 0] <- "N_lim"
df_ind$nut_class[df_ind$N_class == 0 & df_ind$P_class == 0 & df_ind$K_class == 0] <- "no_lim"
df_ind$counts <- 1

nutlim               <- aggregator(df_ind, "counts", c("nut_class", "PR_FI"), "count")
nutlim$Function      <- NULL
nutlim               <- as.data.frame(dcast(nutlim, PR_FI ~ nut_class, value.var = "counts"))
nutlim               <- sort_df(nutlim, c("PR_FI", "no_lim", "N_lim", "P_lim", "K_lim", "NP_lim", "NK_lim", "PK_lim", "NPK_lim"))
is.na(nutlim)        <- 0

nuts                 <- df_ind[is.na(df_ind$N) == F & df_ind$nut_class != "N_lim" ,] # removed N_lim since only 1 point
nuts$classes_ordered <- factor(nuts$nut_class, levels = c("no_lim", "N_lim", "P_lim", "K_lim", "NP_lim", "NK_lim", "PK_lim", "NPK_lim"))

# Nutrient limitations boxplot 
ggplot(nuts, aes(x=classes_ordered, y=NDRE/10000)) +
  geom_boxplot(outlier.colour = NA) + # outliers already visualized in jitter..
  geom_jitter(width=0.25, alpha=0.7, aes(colour = PR_FI)) +
  xlab("Limitation categories") + ylab("NDRE") + guides(colour=guide_legend(title="Fields")) +
  theme_bw()

# Nutrient limitations boxplot (to show variation)
ggplot(nuts, aes(x=classes_ordered, y=NDRE/10000)) +
  geom_boxplot(outlier.colour = NA) + # outliers already visualized in jitter..
  geom_jitter(width=0.25, alpha=0.7, aes(colour = TRT)) + 
  xlab("Limitation categories") + ylab("NDRE") + guides(colour=guide_legend(title="Treatment")) +
  theme_bw()


# WDRVI / NDRE discrepancy
VI_difs <- aggregator(df_ind, c("NDRE", "WDRVI", "NIR_avg", "RE_avg", "TGF_Palm", "LA_F", "RL", "DW_prod"), "PR_FI", "mean")
VI_difs <- na.omit(VI_difs)

ggplot(VI_difs) + 
  geom_smooth(aes(x=WDRVI/1000, y=NDRE/10000), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=WDRVI/1000, y=NDRE/10000, size=RE_avg/10000, color=LA_F), alpha=0.9, na.rm=T) +
  ylab("NDRE") + xlab("WDRVI") + labs(size="Red Edge", color="Leaf area / Frond (m2)") +
  scale_color_viridis(direction = 1, begin=0, end=0.98,discrete = F) +
  geom_text(aes(x=WDRVI/1000, y=NDRE/10000, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 

ggplot(VI_difs) + 
  geom_smooth(aes(x=LA_F, y=NDRE/10000), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=LA_F, y=NDRE/10000, size=RE_avg/10000, color=DW_prod), alpha=0.9, na.rm=T) +
  ylab("NDRE") + xlab("Leaf area / frond (m2)") + labs(size="Red Edge", color="Dw prod.") +
  scale_color_viridis(direction = 1, begin=0, end=0.98,discrete = F) +
  geom_text(aes(x=LA_F, y=NDRE/10000, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 

ggplot(VI_difs) + 
  geom_smooth(aes(x=LA_F, y=WDRVI/10000), formula=y~x, color="black", alpha=0.15, method="lm", se = T) +
  geom_point(aes(x=LA_F, y=WDRVI/10000, size=RE_avg/10000, color=DW_prod), alpha=0.9, na.rm=T) +
  ylab("WDRVI") + xlab("Leaf area / frond (m2)") + labs(size="Red Edge", color="Dw prod.") +
  scale_color_viridis(direction = 1, begin=0, end=0.98,discrete = F) +
  geom_text(aes(x=LA_F, y=WDRVI/10000, label=PR_FI), hjust=0, vjust=-1, cex=3) + 
  theme_bw() 


# RI_F6 and SS_F7 causing this


#############################################################################################
# Inspect relationships between concentrations and absolute values

cov_N  <- sd(df_ind$N,  na.rm=T)/mean(df_ind$N,  na.rm=T)
cov_P  <- sd(df_ind$P,  na.rm=T)/mean(df_ind$P,  na.rm=T)
cov_K  <- sd(df_ind$K,  na.rm=T)/mean(df_ind$K,  na.rm=T)
cov_Mg <- sd(df_ind$Mg, na.rm=T)/mean(df_ind$Mg, na.rm=T)
cov_Ca <- sd(df_ind$Ca, na.rm=T)/mean(df_ind$Ca, na.rm=T)
cov_Fr <- sd(df_ind$Frond_kg_Tot, na.rm=T)/mean(df_ind$Frond_kg_Tot, na.rm=T)

sel <- df_ind[is.na(df_ind$N) == F,]
sel$N_Fr <- sel$N * sel$Frond_kg_Lf


pdf("../Data/3_Output/18_Relations_concentrations_biomass.pdf", width=12, height=8)

ggplot(sel) + 
  geom_point(aes(x=N, y=N_Fr, color=PR_FI), size = 2, alpha=0.5, na.rm=T) +
  xlab("N") + ylab("N_Fr") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=P, y=P_Fr, color=PR_FI), size = 2, alpha=0.5, na.rm=T) +
  xlab("P") + ylab("P_Fr") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=K, y=K_Fr, color=PR_FI), size = 2, alpha=0.5, na.rm=T) +
  xlab("K") + ylab("K_Fr") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=Mg, y=Mg_Fr, color=PR_FI), size = 2, alpha=0.5, na.rm=T) +
  xlab("Mg") + ylab("Mg_Fr") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=Ca, y=Ca_Fr, color=PR_FI), size = 2, alpha=0.5, na.rm=T) +
  xlab("Ca") + ylab("Ca_Fr") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)

ggplot(sel) + 
  geom_smooth(aes(x=N_Fr, y=Frond_kg_Tot, color=PR_FI), formula=y~x, alpha=0.15, method="lm", se = F) +
  geom_point(aes(x=N_Fr, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("N_Fr") + ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_smooth(aes(x=P_Fr, y=Frond_kg_Tot, color=PR_FI), formula=y~x, alpha=0.15, method="lm", se = F) +
  geom_point(aes(x=P_Fr, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("P_Fr") + ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_smooth(aes(x=K_Fr, y=Frond_kg_Tot, color=PR_FI), formula=y~x, alpha=0.15, method="lm", se = F) +
  geom_point(aes(x=K_Fr, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("K_Fr") + ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_smooth(aes(x=Mg_Fr, y=Frond_kg_Tot, color=PR_FI), formula=y~x, alpha=0.15, method="lm", se = F) +
  geom_point(aes(x=Mg_Fr, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("Mg_Fr") + ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_smooth(aes(x=Ca_Fr, y=Frond_kg_Tot, color=PR_FI), formula=y~x, alpha=0.15, method="lm", se = F) +
  geom_point(aes(x=Ca_Fr, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("Ca_Fr") + ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)

ggplot(sel) + 
  geom_point(aes(x=N, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("N") +  ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=P, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("P") +  ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=K, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("K") +  ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=Mg, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("Mg") +  ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)
ggplot(sel) + 
  geom_point(aes(x=Ca, y=Frond_kg_Tot, color=PR_FI), size=2, alpha=0.5, na.rm=T) +
  xlab("Ca") +  ylab("Frond_kg_Tot") + labs(color="Field") + theme_bw() + facet_wrap(~TRT)

dev.off()





