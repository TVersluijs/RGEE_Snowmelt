####################################################################################################################################
#
#                                       Sentinel 2 - Snow detection and masking functions
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a Sentinel-2 image collection to detect pixels that are
#covered by snow.Snow is detected when NDSI > the specified NDSI_threshold.

####################################################################################################################################

#(I): Define a function that extracts which pixels are marked as snow
 computeSnow = function(img){
    
    #(A): Select the NDSI band
     ndsi <- img$select('NDSI')
    
    #(B): Create a binary layer using logical operations.
     snow <- ndsi$gt(NDSI_threshold)$rename('SNOW')
    
    #(C): Return the binary snow parameter
     return(img$addBands(snow))
     }
 
#(II): Calculate fraction of pixels within region of interest that are marked as snow (excluding cloud-masked pixels):
 AddSnowFraction = function(img) {
    
    #(A): Extract band containing the calculated cloud property
     snow <- img$select('SNOW')
    
    #(B): Calculate the fraction of pixels within aoi that are marked as snow
     SnowFraction = snow$reduceRegion(
       reducer = ee$Reducer$mean(),
       geometry = aoi_Shapefile,
       scale = resolution,
       maxPixels=100000000
       )$get('SNOW')
    
      #Note that the Reducer function ee$Reducer$mean() results in a fraction of pixels covered in snow as the band
      #'SNOW' contains a value of 1 for pixels marked as snow and 0 for those marked as non-snow. In addition, all
      #reducer functions do not include pixels where the mask is equal to 0. We have specified as cloud/water mask
      #previously to each image in the image collection, thus pixels containing clouds (with a mask equal to 0) are 
      #not used in the reduction (i.e. in calculating the mean, i.e. fraction of snow covered pixels)
    
    #(C): Add SnowFraction as an image property to each image
     return(img$set('SnowFraction', SnowFraction))
     }
 
#(III): If there are empty images within the collection, then the property SnowFraction is missing for those images. We manually 
#have to add an empty SnowFraction property to those images. We therefore set the SnowFraction value to a no data value 
#of -9999 for all images where the SnowFraction value is NULL.  
 AddNULLSNOW = function(img){
   
   #(A): Set the SnowFraction value to a no data value of -9999 if equal to NULL
    SnowFraction <- ee$List(list(img$get('SnowFraction'), -9999))$reduce(ee$Reducer$firstNonNull())
    
   #(B): Return the updated SnowFraction property  
    return(img$set("SnowFraction", SnowFraction))
    } 

 ###################################################################################################################################