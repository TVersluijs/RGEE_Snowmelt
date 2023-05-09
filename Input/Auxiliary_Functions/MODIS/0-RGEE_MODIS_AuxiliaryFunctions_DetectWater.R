####################################################################################################################################
#
#                                       modis - Water detection and masking functions
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a modis image collection to detect and mask pixels that are
#covered by permanent (i.e. year-round) waterbodies.

####################################################################################################################################


#(I): Detect water pixels using the MODIS XXX map with a xx meter resolution
 compute_Water_MODIS <- function(){
   
   #Load ESA Worldcover map
   MODIS_MOD44W <- ee$ImageCollection('MODIS/006/MOD44W')$first()
   
   #Select all pixels that correspond to permanent waterbodies (index 80)
   MODIS_water <- MODIS_MOD44W$select('water_mask')
   
   #Return water mask
   return(MODIS_water)
   
 }


#(II): Mask permanent water-pixels from all images in the image collection
 Add_WaterMask_MODIS = function(img){
   
   #Return the binary water parameter
   return(img$updateMask(water_mask$lt(1)))
   
 }
 

 

#######################################################################################################################################