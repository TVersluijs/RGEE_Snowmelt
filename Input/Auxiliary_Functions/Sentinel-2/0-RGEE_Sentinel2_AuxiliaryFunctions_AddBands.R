####################################################################################################################################
#
#                                       Sentinel 2 - Calculate new image bands
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a Sentinel-2 image collection to add new bands to each image.
#At this moment, four different Normalized Difference functions are depicted.

####################################################################################################################################

#NOTE ON RESAMPLING

  #SPATIAL RESOLUTIONS OF L2A-BANDS:
  # -B3 (Green)   = 10m 
  # -B4 (Red)     = 10m 
  # -B8 (NIR)     = 10m 
  # -B11 (SWIR 1) = 20m 
  
  #SPATIAL RESOLUTION OF NORMALIZED DIFFERENCE FUNCTIONS:
  # -NDSI: band B3 & B11  --> 20m
  # -NDVI: band B8 & B4   == 10m
  # -NDMI: band B8 & B11  --> 20m
  # -NDWI: band B3 & B8   == 10m
  
  #Thus for NDSI, and NDMI we resample the B3 and B8 band respectively to match the 20m resolution of the B11 band.
  #The default resampling is 'nearest-neighbour' in Google Earth Engine. Instead, we specify 'bicubic' using the
  #resample function. As both FSC measures are derived from NDSI, they also inherit a 20m spatial resolution.

####################################################################################################################################

#(I): Normalized Difference Snow Index (NDSI) function for Sentinel 2 images:
  getNDSI <- function(img){
    return(img$addBands(img$expression('(B3 - B11) / (B3 + B11)',
                                       list('B3'  = img$select('B3')$resample('bicubic')$reproject(img$select('B11')$projection()),
                                            'B11' = img$select('B11')))$rename('NDSI')$copyProperties(img, img$propertyNames())))}

#(II): Normalized Difference Vegetation Index (NDVI) function for Sentinel 2 images:
  getNDVI <- function(img){
    return(img$addBands(img$expression('(B8 - B4) / (B8 + B4)',
                                       list('B4'=img$select('B4'),
                                            'B8'=img$select('B8')))$rename('NDVI')$copyProperties(img, img$propertyNames())))}

#(III): Normalized Difference Moisture Index (NDMI) function for Sentinel 2 images:
  getNDMI <- function(img){
    return(img$addBands(img$expression('(B8 - B11) / (B8 + B11)',
                                       list('B8'=img$select('B8')$resample('bicubic')$reproject(img$select('B11')$projection()),
                                            'B11'=img$select('B11')))$rename('NDMI')$copyProperties(img, img$propertyNames())))}

#(IV): Normalized Difference Water Index (NDWI) function for Sentinel 2 images:
  getNDWI <- function(img){
    return(img$addBands(img$expression('(B3 - B8) / (B3 + B8)',
                                       list('B3'=img$select('B3'),
                                            'B8'=img$select('B8')))$rename('NDWI')$copyProperties(img, img$propertyNames())))}

#(V): Fractional Snow Cover (FSC) function for Sentinel 2 images:
  get_FSC <- function(img){
    
    #Gascoin (2020) Remote Sensing:
    x <- img$expression('2.65 * ndsi - 1.42', list('ndsi'=img$select('NDSI')))
    FSC_Gascoin2020 <- img$expression('0.5 * ((exp((x)) - exp(-(x))) / (exp((x)) + exp(-(x)))) + 0.5', list('x'=x))$rename('FSC_Gascoin2020')
    
    #Aalstad (2020) Remote Sensing of Environment
    FSC_Aalstad2020 <- img$expression('(1.45 * ndsi) - 0.01', list('ndsi'=img$select('NDSI')))$rename('FSC_Aalstad2020')
    FSC_Aalstad2020 <- FSC_Aalstad2020$where(FSC_Aalstad2020$lt(0), ee$Image(0))
    FSC_Aalstad2020 <- FSC_Aalstad2020$where(FSC_Aalstad2020$gt(1), ee$Image(1))
    
    #Return both measures as separate bands
    return(img$addBands(FSC_Gascoin2020)$
               addBands(FSC_Aalstad2020)$
               copyProperties(img, img$propertyNames()))
    
    }
  

######################################################################################################################################