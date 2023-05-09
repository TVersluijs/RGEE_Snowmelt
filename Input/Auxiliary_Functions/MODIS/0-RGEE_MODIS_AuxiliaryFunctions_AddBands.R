####################################################################################################################################
#
#                                         MODIS - Calculate new image bands
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a MODIS image collection to add new bands to each image.
#At this moment, four different Normalized Difference functions are depicted.

####################################################################################################################################

#(I): Normalized Difference Snow Index (NDSI) function for MODIS images:
  getNDSI <- function(img){
    return(img$addBands(img$expression('(B4 - B6) / (B4 + B6)',
                                       list('B4'=img$select('sur_refl_b04'),
                                            'B6'=img$select('sur_refl_b06')))$rename('NDSI')$copyProperties(img, img$propertyNames())))}

#(II): Normalized Difference Vegetation Index (NDVI) function for MODIS images:
  getNDVI <- function(img){
    return(img$addBands(img$expression('(B2 - B1) / (B2 + B1)',
                                       list('B2'=img$select('sur_refl_b02'),
                                            'B1'=img$select('sur_refl_b01')))$rename('NDVI')$copyProperties(img, img$propertyNames())))}

#(III): Normalized Difference Moisture Index (NDMI) function for MODIS images:
  getNDMI <- function(img){
    return(img$addBands(img$expression('(B2 - B6) / (B2 + B6)',
                                       list('B2'=img$select('sur_refl_b02'),
                                            'B6'=img$select('sur_refl_b06'))) $rename('NDMI')$copyProperties(img, img$propertyNames())))}

#(IV): Normalized Difference Water Index (NDWI) function for MODIS images:
  getNDWI <- function(img){
    return(img$addBands(img$expression('(B4 - B2) / (B4 + B2)',
                                       list('B4'=img$select('sur_refl_b04'),
                                            'B2'=img$select('sur_refl_b02')))$rename('NDWI')$copyProperties(img, img$propertyNames())))}


######################################################################################################################################