####################################################################################################################################
#
#                                       Sentinel 2 - Water detection and masking functions
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a Sentinel-2 image collection to detect and mask pixels that are
#covered by permanent (i.e. year-round) waterbodies.

####################################################################################################################################

#(I.a): (Option 1) Detect water pixels using the ESA Worldcover map with a 10 meter resolution
 compute_Water_ESA <- function(ESA_index=80){
   
   #Load ESA Worldcover map
   ESA_WorldCover <- ee$ImageCollection("ESA/WorldCover/v100")$first()
   
   #Select all pixels that correspond to permanent waterbodies (index 80)
   ESA_water <- ESA_WorldCover$select('Map')$eq(ESA_index)
   
   #Return water mask
   return(ESA_water)
   
 }


#(I.b): (Option 2) Detect water pixels using a manual function based on NDWI, NDSI and NIR bands
 compute_Water_Manual <- function(img_col=s2_col, start_date_NDWI=start_date_NDWI, end_date_NDWI=end_date_NDWI, NDWI_threshold=NDWI_threshold,
                                  start_date_NDSI=start_date_NDSI, end_date_NDSI=end_date_NDSI, NDSI_threshold=NDSI_threshold, NIR_threshold=NIR_threshold){
   
   #Water can in general be detected by looking at NDWI values (>0 indicates water). This method works well in non build-up areas (Du 2016, Xu 2006)
   #However, high NDWI values per se do not indicate water as snow can also have high NDWI values. To properly detect waterbodies, we thus have to mask
   #snowcovered areas first. Snow can be differentiated with water by looking at NIR reflectance values as water has NIR values < 0.15, while
   #snow has NIR values > 0.15. Thus, to mask waterbodies we follow the following steps:
   
   #(A): select all images between start_date_NDWI and end_date_NDWI with a cloudcover percentage < 5% (s2_NDWI_mask)
   #(B): Calculate a median NDWI image for all images in s2_NDWI_mask
   #(C): select all pixels with median_NDWI > 0.0 (these are WATER or SNOW pixels)
   #(D): select all images between start_date_NDSI and end_date_NDSI with a cloudcover percentage < 5% (s2_NDSI_mask)
   #(E): Calculate median NDSI and NIR (B8) images for all images in s2_NDSI_mask
   #(F): Select all pixels with a median_NDSI > 0.42 and a median_B8 > 1000 (these are SNOW pixels)
   #(G): Subtract the snow-covered pixels from the water-and-snow pixels in step 13 to retain water-pixels only
   #(H): Return these waterpixels as binary output to be used for masking water pixels
   
   #(A): select all images between start_date_NDWI and end_date_NDWI with a cloudcover percentage < 5% (s2_NDWI_mask)
   s2_NDWI_mask <- img_col$
     #Filter dates for which a median NDWI image will be generaged
     filterDate(start_date_NDWI, end_date_NDWI)$
     #Filter selection to only contain images with less than 5% cloudcover (CloudFraction=0.05)
     filter(ee$Filter$lt('CloudFraction', 0.05))$
     #Apply cloudmask for individual pixels
     map(Add_CloudMask)
   
   #(B): Calculate a median NDWI image for all images in s2_NDWI_mask
   s2_NDWI_mask <- s2_NDWI_mask$select("B4", "B3", "B2", "NDWI")$reduce(ee$Reducer$median())
   ndwi_median <-  s2_NDWI_mask$select("NDWI_median")
   
   # #Plot median NDWI image (for debugging)
   # ndwiViz <- list(min=NDWI_threshold, max = 1, palette = c('ffffff','0524ff','0dffff','000000'))
   # Map$addLayer(s2_NDWI_mask, list(bands=c("B4_median", "B3_median", "B2_median"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR Median')+
   # Map$addLayer(s2_NDWI_mask$select("NDWI_median"), ndwiViz, 'NDWI_median')
   
   #(C): select all pixels with median_NDWI > NDWI_threshold (these are WATER or SNOW pixels)
   ndwi_mask <- ndwi_median$gt(NDWI_threshold)
   Map$addLayer(ndwi_mask, list(min=0, max=1, palette = c('ffffff', '000000')), 'NDWI_mask')
   
   #(D): select all images between start_date_NDSI and end_date_NDSI with a cloudcover percentage < 5% (s2_NDSI_mask)
   s2_NDSI_mask <- img_col$
     #Filter dates for which a median NDWI image will be generaged
     filterDate(start_date_NDSI, end_date_NDSI)$
     #Filter selection to only contain images with less than 5% cloudcover (CloudFraction=0.05)
     filter(ee$Filter$lt('CloudFraction', 0.05))$
     #Apply cloudmask for individual pixels
     map(Add_CloudMask)
   
   #(E): Calculate median NDSI and NIR (B8) images for all images in s2_NDSI_mask
   s2_NDSI_mask <- s2_NDSI_mask$select("B4", "B3", "B2", "NDSI", "B8")$reduce(ee$Reducer$median())
   ndsi_median <-  s2_NDSI_mask$select("NDSI_median")
   B8_median <-  s2_NDSI_mask$select("B8_median")
   
   # #Plot median NDSI image (for debugging)
   # ndsiViz <- list(min=0.42, max = 1, palette = c('ffffff','0524ff','0dffff','000000'))
   # Map$addLayer(s2_NDSI_mask, list(bands=c("B4_median", "B3_median", "B2_median"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR Median')+
   # Map$addLayer(s2_NDSI_mask$select("NDSI_median"), ndsiViz, 'NDSI_median')
   
   # #Plot Median NIR (B8) image (for debugging)
   # B8Viz <- list(min=NIR_threshold, max = NIR_threshold+1, palette = c('000000', 'ffffff'))
   # Map$addLayer(s2_NDSI_mask, list(bands=c("B4_median", "B3_median", "B2_median"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR Median')+
   # Map$addLayer(s2_NDSI_mask$select("B8_median"), B8Viz, 'B8_median')
   
   #(F): Select all pixels with a median_NDSI > NDSI_threshold and a median_B8 > NIR_threshold (these are SNOW pixels)
   ndsi_mask <- ndsi_median$gt(NDSI_threshold)
   B8_mask <- B8_median$gt(NIR_threshold)
   snow_mask <- ndsi_mask$And(B8_mask)$rename('snow_mask')
   
   #(G): Subtract the snow-covered pixels from the water-and-snow pixels (NDWI>0) to retain water-pixels only
   water_mask <- ndwi_mask$And(snow_mask$lt(1))
   
   # #Plot stepwise selection of waterbodies in median image (for debugging)
   # Map$addLayer(s2_NDWI_mask, list(bands=c("B4_median", "B3_median", "B2_median"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR Median')+
   # Map$addLayer(ndwi_mask,list(min=0, max = 1, palette = c('ffffff', '000000')), 'NDWI > 0')+
   # Map$addLayer(snow_mask,list(min=0, max = 1, palette = c('ffffff', '000000')), 'SNOW')+
   # Map$addLayer(water_mask,list(min=0, max = 1, palette = c('ffffff', '000000')), 'Water')
   
   #(H): Return water mask
   return(water_mask)
    
 }
 
 
#(II): Mask permanent water-pixels from all images in the image collection
 Add_WaterMask = function(img){
   
   #Return the binary water parameter
   return(img$updateMask(water_mask$lt(1)))
   
 }
 

 

#######################################################################################################################################