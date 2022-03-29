# Determine suitable image dependent threshold for soil cover
# Bastiaen Boekelo, March 2022

library(raster)
library(rgdal)

ras_raw <- raster("C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/1_Input/UAV/RI/3_Raster_derivations_RI/NDRE_RI/NDRE_Micasense_UNL_RI_F2.tif")
shp     <- readOGR("C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/1_Input/UAV/RI/5_Treatment_Boundary_RI", "RI_F2-Polygone_Ahmad_Kusairi")
ras     <- mask(ras_raw, shp)
plot(ras)

break_min <- as.numeric(quantile(ras, 0.01))
break_max <- as.numeric(quantile(ras, 0.99))
total     <- cellStats(ras > -9999, stat = "sum")

nr_breaks <- 20
breaks    <- round(seq(break_min, break_max, (break_max - break_min)/nr_breaks ), 5)[1:nr_breaks]

for(i in 1:length(breaks[1:(nr_breaks)])){

  BRK        <- breaks[i]
  ras_break  <-  ras < BRK
  cover      <-  round(cellStats(ras_break, stat = "sum") / total * 100, 1)
  
  if(i == 1){
    covers <- cover
  } else {
    covers <- c(covers, cover)
  }

  print(paste0(BRK, ": ", cover, "%"))
  
}

# Define 2nd derivative
deriv           <- function(x, y) {diff(y) / diff(x)}
middle_pts      <- function(x)    {x[-1] - diff(x) / 2}
second_d        <- deriv(middle_pts(breaks), deriv(breaks, covers))
smooth_2nd_d    <- loess(second_d ~ midpts, data.frame(second_d = second_d, midpts = middle_pts(middle_pts(breaks))), model = T)

plot(breaks, covers, xlab = "NDRE_breakpoint", ylab="% covered")
plot(middle_pts(middle_pts(breaks)), deriv(middle_pts(breaks), deriv(breaks, covers)), ylab = "2nd derivative", xlab="NDRE_breakpoint")

derivs          <- unlist(unlist(smooth_2nd_d)[2:(nr_breaks-1)])

# Determine breakpoints & cut-off raster
highest_slope   <- max(derivs)
break_point     <- match(highest_slope,derivs) + 1 

break_value1    <- breaks[break_point - 3]
break_value2    <- breaks[break_point]
break_value3    <- breaks[break_point + 3]

canopy_cover1   <- ras > break_value1
canopy_cover2   <- ras > break_value2
canopy_cover3   <- ras > break_value3

writeRaster(canopy_cover1, "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/0_Temp/soilcover_RI_F2_v1.tif", overwrite=T)
writeRaster(canopy_cover2, "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/0_Temp/soilcover_RI_F2_v2.tif", overwrite=T)
writeRaster(canopy_cover3, "C:/Users/boeke015/OneDrive - Wageningen University & Research/UAV_Palm/Data/0_Temp/soilcover_RI_F2_v3.tif", overwrite=T)
