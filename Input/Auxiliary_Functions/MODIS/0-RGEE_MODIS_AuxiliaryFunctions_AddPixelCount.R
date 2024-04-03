#######################################################################################################################################
#
#                                           MODIS - Pixel counting functions
#
#######################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a MODIS image collection to count the number of unmasked pixels.

#######################################################################################################################################

#Count the number of unmasked pixels within aoi_Shapefile
AddPixelCount <- function(img){
  
  #Count the number of unmasked pixels within aoi_Shapefile
  pixels_unmasked <- img$select('NDSI')$reduceRegion(
    reducer = ee$Reducer$count(),
    geometry = aoi_Shapefile,
    scale = resolution,
    maxPixels = 100000000
    )$get('NDSI')
  
  #Count the total number of pixels within aoi_Shapefile
  pixels_total <- img$select('NDSI')$unmask()$reduceRegion(
    reducer = ee$Reducer$count(),
    geometry = aoi_Shapefile,
    scale = resolution,
    maxPixels = 100000000
    )$get('NDSI')
  
  #Add pixel counts within aoi_Shapefile as image properties to each image
  return(img$
           set('unmasked', pixels_unmasked)$
           set('total', pixels_total))
  
}

#If there are empty images within the collection, then the pixel_count properties are missing for those images. We manually 
#have to add empty pixel count properties to those images. We therefore set the pixels_unmasked and pixels_total value to a
#no data value of -9999 for all images where these values were NULL.  
AddNULLPixelCount = function(img){
  
  #(A): Set the pixels_unmasked value to a no data value of -9999 if equal to NULL
  pixels_unmasked <- ee$List(list(img$get('unmasked'), -9999))$reduce(ee$Reducer$firstNonNull())
  
  #(B): Set the pixels_total value to a no data value of -9999 if equal to NULL  
  pixels_total <- ee$List(list(img$get('total'), -9999))$reduce(ee$Reducer$firstNonNull())
  
  #(C): Return the updated clouds fraction properties
  return(img$
           set("unmasked", pixels_unmasked)$
           set("total", pixels_total))
  
}