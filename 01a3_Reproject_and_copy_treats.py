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

for PROV in provinces:

    ## LIST FILES
    treats = [f for f in os.listdir(wd + "1_Input\\UAV\\" + PROV + "\\5_Treatment_Boundary_" + PROV + "\\") if re.match(r'.*\.shp$', f)]

    fields = range(1,9)

    for FIELD in fields:

        FIELDNAME = "F" + str(FIELD)
        
        ## INPUT - OUTPUT - VARIABLES

        # Only perform action if a field is present in all the four input folders
        if any(FIELDNAME in s for s in treats):

            print(PROV + "_" + FIELDNAME)

            name_trt   = [s for s in treats if FIELDNAME in s][0]

            # Inputs
            input_shp        = wd + "1_Input\\UAV\\" + PROV + "\\5_Treatment_Boundary_" + PROV + "\\" + name_trt

            # Variables
            crs_in           = "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]],VERTCS['unknown',VDATUM['unknown'],PARAMETER['Vertical_Shift',0.0],PARAMETER['Direction',1.0],UNIT['Meter',1.0]]"
            crs_out          = "PROJCS['UTM Zone 49, Southern Hemisphere',GEOGCS['Geographic Coordinate System',DATUM['WGS84',SPHEROID['WGS84',6378137.0,298.257223560493]],PRIMEM['Greenwich',0.0],UNIT['degree',0.0174532925199433]],PROJECTION['Transverse_Mercator'],PARAMETER['false_easting',500000.0],PARAMETER['false_northing',10000000.0],PARAMETER['central_meridian',111.0],PARAMETER['scale_factor',0.9996],PARAMETER['latitude_of_origin',0.0],UNIT['Meter',1.0]]"

            # Intermediates
            reprojected_shp  = wd + "2_Intermediate\\01a3_FieldMask\\01_Treats_" + PROV + "_" + FIELDNAME + ".shp"

            # Reproject palm points
            arcpy.Project_management(input_shp, reprojected_shp, crs_in, "", crs_out, "NO_PRESERVE_SHAPE", "", "NO_VERTICAL")
























