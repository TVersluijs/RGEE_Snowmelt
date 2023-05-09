####################################################################################################################################
#
#                                       Sentinel 2 - Calculate new image bands
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a Sentinel-2 image collection to add new bands to each image.
#At this moment, four different Normalized Difference functions are depicted.

####################################################################################################################################

#(I): Normalized Difference Snow Index (NDSI) function for Sentinel 2 images:
  getNDSI <- function(img){
    return(img$addBands(img$expression('(B3 - B11) / (B3 + B11)',
                                       list('B11' = img$select('B11'),
                                            'B3'  = img$select('B3')))$rename('NDSI')$copyProperties(img, img$propertyNames())))}

#(II): Normalized Difference Vegetation Index (NDVI) function for Sentinel 2 images:
  getNDVI <- function(img){
    return(img$addBands(img$expression('(B8 - B4) / (B8 + B4)',
                                       list('B8'=img$select('B8'),
                                            'B4'=img$select('B4')))$rename('NDVI')$copyProperties(img, img$propertyNames())))}

#(III): Normalized Difference Moisture Index (NDMI) function for Sentinel 2 images:
  getNDMI <- function(img){
    return(img$addBands(img$expression('(B8 - B11) / (B8 + B11)',
                                       list('B8'=img$select('B8'),
                                            'B11'=img$select('B11')))$rename('NDMI')$copyProperties(img, img$propertyNames())))}

#(IV): Normalized Difference Water Index (NDWI) function for Sentinel 2 images:
  getNDWI <- function(img){
    return(img$addBands(img$expression('(B3 - B8) / (B3 + B8)',
                                       list('B3'=img$select('B3'),
                                            'B8'=img$select('B8')))$rename('NDWI')$copyProperties(img, img$propertyNames())))}


######################################################################################################################################