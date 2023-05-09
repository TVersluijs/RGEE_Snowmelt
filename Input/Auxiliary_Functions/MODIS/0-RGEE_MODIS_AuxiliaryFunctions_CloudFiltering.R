####################################################################################################################################
#
#                                       MODIS - Cloud detection and masking functions
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a MODIS image collection to detect and mask pixels that are
#covered by clouds. This includes cirrus clouds and clouds as detected according to different MODIS algorithms (MOD35 and PGE11).

####################################################################################################################################

#(I): Detect which pixels are considered to be cloudy
 computeClouds <- function(img){
   
   #(A): Select quality band
    qaBand <- img$select("state_1km")
   
   #(B): Detect cirrus clouds (bits 8-9): only average and high cirrus will be filtered; gt(1)
    bit_start=8
    bit_end=9
    summed_bits <- ee$Number(sum(2^seq(bit_start, bit_end, 1)))
    cirrus <- qaBand$bitwiseAnd(summed_bits)$gt(1)$rename('cirrus') 
   
   #(C): MOD35 (bits 0-1), Filter cloudy pixels; eq(1)
    bit_start=0
    bit_end=1
    summed_bits <- ee$Number(sum(2^seq(bit_start, bit_end, 1)))
    MOD35_1 <- qaBand$bitwiseAnd(summed_bits)$eq(1)
    MOD35_clouds <- MOD35_1$rename("MOD35_clouds")
    MOD35_clouds <- MOD35_clouds$Or(cirrus) #Add cirrus to Mod35 cloud detection
   
   #(D): PGE11 (bit 10), Internal cloud algorithm flag: clouds will be filtered; gt(0)
    bit_start=10
    bit_end=10
    summed_bits <- ee$Number(sum(2^seq(bit_start, bit_end, 1)))
    PGE11_clouds <- qaBand$bitwiseAnd(summed_bits)$gt(0)$rename('PGE11_clouds')
    PGE11_clouds <- PGE11_clouds$Or(cirrus) #Add cirrus to PGE11 cloud detection
   
   #(E): Combine all algorithms into a single binary 'Combined_clouds' parameter
    Combined_clouds <- MOD35_clouds$Or(PGE11_clouds)$rename('Combined_clouds')
   
   #(F): Create an inverse cloud band for all three options
    MOD35_clouds_inv <- MOD35_clouds$lt(1)$rename('MOD35_clouds_inv')
    PGE11_clouds_inv <- PGE11_clouds$lt(1)$rename('PGE11_clouds_inv')
    Combined_clouds_inv <- Combined_clouds$lt(1)$rename('Combined_clouds_inv')
    
   #(G): Return the binary cloud parameter
    return(img$
             addBands(MOD35_clouds$gt(0))$
             addBands(PGE11_clouds$gt(0))$
             addBands(Combined_clouds$gt(0))$
             addBands(MOD35_clouds_inv)$
             addBands(PGE11_clouds_inv)$
             addBands(Combined_clouds_inv))
    }

#(II): Calculate fraction of pixels within region of interest that are marked as clouds:
 AddCloudFraction = function(img){
   
   #(A): Extract band containing the calculated cloud property
    Combined_clouds <- img$select('Combined_clouds')
    PGE11_clouds <- img$select('PGE11_clouds')
    MOD35_clouds <- img$select('MOD35_clouds')
   
   #(B): Calculate the fraction of pixels within aoi_Shapefile that are marked as either PGE11 or MOD35 clouds
    Combined_clouds_Fraction = Combined_clouds$reduceRegion(
      reducer = ee$Reducer$mean(),
      geometry = aoi_Shapefile,
      scale = resolution_qaband #Quality band has a resolution of 1km
      )$get('Combined_clouds')
   
   #(C): Calculate the fraction of pixels within aoi_Shapefile that are marked as PGE35 clouds
    PGE11_clouds_Fraction = PGE11_clouds$reduceRegion(
      reducer = ee$Reducer$mean(),
      geometry = aoi_Shapefile,
      scale = resolution_qaband #Quality band has a resolution of 1km
      )$get('PGE11_clouds')
   
   #(D): Calculate the fraction of pixels within aoi_Shapefile that are marked as MOD35 clouds
    MOD35_clouds_Fraction = MOD35_clouds$reduceRegion(
      reducer = ee$Reducer$mean(),
      geometry = aoi_Shapefile,
      scale = resolution_qaband #Quality band has a resolution of 1km
      )$get('MOD35_clouds')
   
   #(E): Add CloudFraction within aoi_Shapefile as an image property to each image
    return(img$
             set('Combined_clouds_Fraction', Combined_clouds_Fraction)$
             set('PGE11_clouds_Fraction', PGE11_clouds_Fraction)$
             set('MOD35_clouds_Fraction', MOD35_clouds_Fraction))
    }

#(III): If there are empty images within the collection, then the property clouds_Fraction is missing for those images. We manually 
#have to add an empty clouds_Fraction property to those images. We therefore set the clouds_Fraction value to a no data value 
#of -9999 for all images where the clouds_Fraction value is NULL.  
 AddNULLCloudFraction = function(img){
   
   #(A): Set the MOD35_clouds_Fraction value to a no data value of -9999 if equal to NULL
    MOD35_clouds_Fraction <- ee$List(list(img$get('MOD35_clouds_Fraction'), -9999))$reduce(ee$Reducer$firstNonNull())
   
   #(B): Set the PGE11_clouds_Fraction value to a no data value of -9999 if equal to NULL  
    PGE11_clouds_Fraction <- ee$List(list(img$get('PGE11_clouds_Fraction'), -9999))$reduce(ee$Reducer$firstNonNull())
   
   #(C): Set the Combined_clouds_Fraction value to a no data value of -9999 if equal to NULL   
   Combined_clouds_Fraction <- ee$List(list(img$get('Combined_clouds_Fraction'), -9999))$reduce(ee$Reducer$firstNonNull())
   
   #(D): Return the updated clouds fraction properties
   return(img$
             set("MOD35_clouds_Fraction", MOD35_clouds_Fraction)$
             set("PGE11_clouds_Fraction", PGE11_clouds_Fraction)$
             set("Combined_clouds_Fraction", Combined_clouds_Fraction))
   }

#(IV): Mask all cloudy pixels in all images in the image collection for the selected cloud masking algorithm
 AddCloudMask = function(img){
   
   #(A): Extract cloud information
    cloudmask_MODIS <- img$select(cloud_algorithm)
   
   #(B): Mask cloudy pixels with value zero
    return(img$updateMask(cloudmask_MODIS$lt(1)))
    }

 
#######################################################################################################################################