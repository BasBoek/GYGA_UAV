# Bastiaen Boekelo, March 2021
# Goal1: Reduce size of all data to clip it only at palm tree level
# Goal2: Store new reduced data in a consistent way for further processing

print("Preparing...")
import os
import re
import arcpy

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension('Spatial')

# Set Geoprocessing environments
arcpy.env.newPrecision = "SINGLE"
arcpy.env.autoCommit = "1000"
arcpy.env.XYResolution = ""
arcpy.env.processingServerUser = ""
arcpy.env.XYDomain = ""
arcpy.env.processingServerPassword = ""
arcpy.env.cartographicPartitions = ""
arcpy.env.terrainMemoryUsage = "false"
arcpy.env.MTolerance = ""
arcpy.env.compression = "LZ77"
arcpy.env.coincidentPoints = "MEAN"
arcpy.env.randomGenerator = "0 ACM599"
arcpy.env.outputCoordinateSystem = ""
arcpy.env.rasterStatistics = "STATISTICS 1 1"
arcpy.env.ZDomain = ""
arcpy.env.transferDomains = "false"
arcpy.env.maintainAttachments = "true"
arcpy.env.resamplingMethod = "NEAREST"
arcpy.env.snapRaster = ""
arcpy.env.projectCompare = "NONE"
arcpy.env.cartographicCoordinateSystem = ""
arcpy.env.configKeyword = ""
arcpy.env.outputZFlag = "Same As Input"
arcpy.env.qualifiedFieldNames = "true"
arcpy.env.tileSize = "128 128"
arcpy.env.parallelProcessingFactor = ""
arcpy.env.pyramid = "PYRAMIDS -1 NEAREST DEFAULT 75 NO_SKIP"
arcpy.env.referenceScale = ""
arcpy.env.processingServer = ""
arcpy.env.extent = "DEFAULT"
arcpy.env.XYTolerance = ""
arcpy.env.tinSaveVersion = "CURRENT"
arcpy.env.nodata = "NONE"
arcpy.env.MDomain = ""
arcpy.env.spatialGrid1 = "0"
arcpy.env.cellSize = "MAXOF"
arcpy.env.outputZValue = ""
arcpy.env.outputMFlag = "Same As Input"
arcpy.env.geographicTransformations = ""
arcpy.env.spatialGrid2 = "0"
arcpy.env.ZResolution = ""
arcpy.env.mask = ""
arcpy.env.spatialGrid3 = "0"
arcpy.env.maintainSpatialIndex = "false"
arcpy.env.MResolution = ""
arcpy.env.derivedPrecision = "HIGHEST"
arcpy.env.ZTolerance = ""


########
## SET WORKING DIRECTORY - MANUALLY
########

wd = "C:\\Users\\boeke015\\OneDrive - Wageningen University & Research\\UAV_Palm\\Data\\"

########
## ANALYSIS
########

print("\nGO!\n")

provinces = ["CK", "JB", "RI", "SS"]
#provinces = ["CK"]
for PROV in provinces:

    ## LIST FILES
    orthos = [f for f in os.listdir(wd + "1_Input\\UAV\\" + PROV + "\\2_Orthophoto_" + PROV + "\\") if re.match(r'.*\.tif$', f)]
    NDVIs  = [f for f in os.listdir(wd + "1_Input\\UAV\\" + PROV + "\\3_Raster_derivations_" + PROV + "\\NDVI_" + PROV + "\\") if re.match(r'.*\.tif$', f)]
    NDREs  = [f for f in os.listdir(wd + "1_Input\\UAV\\" + PROV + "\\3_Raster_derivations_" + PROV + "\\NDRE_" + PROV + "\\") if re.match(r'.*\.tif$', f)]
    points = [f for f in os.listdir(wd + "1_Input\\UAV\\" + PROV + "\\4_Palm_Points_" + PROV + "\\") if re.match(r'.*\.shp$', f)]

    fields = range(1,9)

    for FIELD in fields:

        FIELDNAME = "F" + str(FIELD)
        
        ## INPUT - OUTPUT - VARIABLES

        # Only perform action if a field is present in all the four input folders
        if any(FIELDNAME in s for s in orthos) & any(FIELDNAME in s for s in NDVIs) & any(FIELDNAME in s for s in NDREs) & any(FIELDNAME in s for s in points):

            print(PROV + "_" + FIELDNAME)

            name_ortho = [s for s in orthos if FIELDNAME in s][0]
            name_NDVI  = [s for s in NDVIs  if FIELDNAME in s][0]
            name_NDRE  = [s for s in NDREs  if FIELDNAME in s][0]
            name_pnt   = [s for s in points if FIELDNAME in s][0]

            # Inputs
            input_shp        = wd + "1_Input\\UAV\\" + PROV + "\\4_Palm_Points_" + PROV + "\\" + name_pnt
            input_ortho_tif  = wd + "1_Input\\UAV\\" + PROV + "\\2_Orthophoto_" + PROV + "\\" + name_ortho
            input_NDVI_tif   = wd + "1_Input\\UAV\\" + PROV + "\\3_Raster_derivations_" + PROV + "\\NDVI_" + PROV + "\\" + name_NDVI
            input_NDRE_tif   = wd + "1_Input\\UAV\\" + PROV + "\\3_Raster_derivations_" + PROV + "\\NDRE_" + PROV + "\\" + name_NDRE

            # Variables
            crs_in           = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],VERTCS['unknown',VDATUM['unknown'],PARAMETER['Vertical_Shift',0.0],PARAMETER['Direction',1.0],UNIT['Meter',1.0]]"
            crs_out          = "PROJCS['UTM Zone 49, Southern Hemisphere',GEOGCS['Geographic Coordinate System',DATUM['WGS84',SPHEROID['WGS84',6378137.0,298.257223560493]],PRIMEM['Greenwich',0.0],UNIT['degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',10000000.0],PARAMETER['central_meridian',111.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]"
            dist_integer     = 4
            buf_dist         = str(dist_integer)

            # Temps
            buffered_shp   = wd + "0_Temp\\01_Points_buf_" + buf_dist + "m_"+ PROV + "_" + FIELDNAME +  ".shp"
            mask_NDVI_tif  = wd + "0_Temp\\01_NDVI_"  + PROV + "_" + FIELDNAME + "_Mask_float.tif"
            mask_NDRE_tif  = wd + "0_Temp\\01_NDRE_"  + PROV + "_" + FIELDNAME + "_Mask_float.tif"
            NDVI_10000_tif = wd + "0_Temp\\01_NDVI_"  + PROV + "_" + FIELDNAME + "_10000_float.tif"
            NDRE_10000_tif = wd + "0_Temp\\01_NDRE_"  + PROV + "_" + FIELDNAME + "_10000_float.tif"

            # Intermediates
            output_ortho_tif = wd + "2_Intermediate\\01_Ortho_"  + PROV + "_" + FIELDNAME + "_Mask.tif"
            output_NDVI_tif  = wd + "2_Intermediate\\01_NDVI_"   + PROV + "_" + FIELDNAME + "_Mask.tif"
            output_NDRE_tif  = wd + "2_Intermediate\\01_NDRE_"   + PROV + "_" + FIELDNAME + "_Mask.tif"
            reprojected_shp  = wd + "2_Intermediate\\01_Points_" + PROV + "_" + FIELDNAME + ".shp"
            

            # Reproject palm points
            arcpy.Project_management(input_shp, reprojected_shp, crs_in, "", crs_out, "NO_PRESERVE_SHAPE", "", "NO_VERTICAL")

            # Buffer palm points
            arcpy.Buffer_analysis(reprojected_shp, buffered_shp, str(buf_dist) + " Meters", "FULL", "ROUND", "ALL", "", "PLANAR")

            # Mask orthomosaic
            arcpy.gp.ExtractByMask_sa(input_ortho_tif, buffered_shp, output_ortho_tif)

            # Mask NDVI
            arcpy.gp.ExtractByMask_sa(input_NDVI_tif, buffered_shp, mask_NDVI_tif)
            arcpy.gp.Times_sa(mask_NDVI_tif, "10000", NDVI_10000_tif)
            arcpy.gp.Int_sa(NDVI_10000_tif, output_NDVI_tif)

            # Mask NDRE
            arcpy.gp.ExtractByMask_sa(input_NDRE_tif, buffered_shp, mask_NDRE_tif)
            arcpy.gp.Times_sa(mask_NDRE_tif, "10000", NDRE_10000_tif)
            arcpy.gp.Int_sa(NDRE_10000_tif, output_NDRE_tif)























