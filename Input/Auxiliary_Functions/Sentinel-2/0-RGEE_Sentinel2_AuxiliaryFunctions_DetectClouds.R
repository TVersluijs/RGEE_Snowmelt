####################################################################################################################################
#
#                                       Sentinel 2 - Cloud detection and masking functions
#
####################################################################################################################################

#DESCRIPTION: This script contains functions that can be mapped over a Sentinel-2 image collection to detect and mask pixels that are
#covered by clouds. This includes cirrus clouds and opaque clouds as specified in the QA60 information band.

####################################################################################################################################

#Add the 'S2_CLOUD_PROBABILITY' dataset to our sentinel-2 image collection as a new image property:
 Add_CloudProb_Property <- function(s2_col, s2_cloudless_col){

   #Join cloud probability information from s2_cloudless_col to s2_col
  
     #Use an equals filter to specify how the collections match (here based on 'system:index').
      Filter_system_index <- ee$Filter$equals(
        leftField= 'system:index',
        rightField= 'system:index')
  
     #Define the join to store matches from the secondary collection as a named property of the 
     #features of the primary collection. For each element in the primary collection the saveFirst 
     #join saves only the first element in the secondary collection that matches the conditions 
     #specified in the Filter_system_index. The parameter matchKey specifies the new property
     #name used to save the match in the primary collection. Unmatched elements in the primary 
     #collection are dropped unless the parameter outer is set to TRUE. 
      save_first_join <- ee$Join$saveFirst(matchKey='s2cloudless', outer=TRUE)
  
     #Apply the join.
      s2_col_cloud <- save_first_join$apply('primary'= s2_col, 
                                            'secondary' = s2_cloudless_col, 
                                            'condition' = Filter_system_index)
  
   #transform s2_col_cloud to image collection
    s2_col_cloud <- ee$ImageCollection(s2_col_cloud)
  
   #Return s2_col_cloud
    return(s2_col_cloud)
  
    }

#(I.a): (Option 1) Detect cloudy pixels as either containing opague clouds (bit 10) or cirrus clouds (bit 11) using the QA60 band 
 compute_Clouds_QA60 <- function(img){
  
   #(A): Select quality band
    qaBand <- img$select("QA60")
  
   #(B): Extract which pixels are not marked as opague clouds (bit 10)
    opague <- qaBand$bitwiseAnd(ee$Number(2^10))$eq(0)
  
   #(C): Extract which pixels are not marked as cirrus clouds (bit 11)
    cirrus <- qaBand$bitwiseAnd(ee$Number(2^11))$eq(0)
  
   #(D): Combine both cloud types into a single binary 'clouds' parameter (where 1 indicates clear conditions and 0 indicates clouds)
    clouds <- opague$And(cirrus)$rename('clouds')
  
   #(E): Return the binary cloud parameter
    return(img$
             addBands(clouds$lt(1))$
             addBands(clouds$gt(0)$rename('clouds_inv')))
    }

#(I.b): (Option 2) Detect cloudy pixels using s2cloudless probability layer and a preset cloud probability threshold
 compute_Clouds_cldprb <- function(img){
   
   #(A): Get s2cloudless image, subset the probability band.
   cld_prb <- ee$Image(img$get('s2cloudless'))$select('probability')$rename("clouds_prob")
   
   #(B): Extract which pixels are defined as clouds (i.e. probability > cldprb_threshold) and store as a binary 'clouds' band
   clouds <- cld_prb$gt(cldprb_threshold)$rename('clouds_unadjusted')

   #(C): Remove small cloud-shadow patches (erosion) and dilate remaining pixels by 'cld_buffer'.
   #20 m scale is for speed, and assumes clouds don't require 10 m precision.
   clouds_buffered <- clouds$
                        focalMin(radius=cld_erosion, units="meters")$
                        focalMax(radius=cld_buffer, units="meters")$
                        reproject(crs=img$select(0)$projection(), scale=resolution_cldmsk)$
                        rename('clouds')
   
   #(D): Create an inverse of the clouds_buffered band
   clouds_buffered_inv <- clouds_buffered$lt(1)$rename('clouds_inv')
   
   #(E): Return the binary cloud parameter
   return(img$addBands(clouds)$
            addBands(clouds_buffered)$
            addBands(clouds_buffered_inv)$
            addBands(cld_prb))
  
 }
 
#(II): Calculate fraction of pixels within region of interest that are marked as clouds:
 Add_CloudFraction = function(img) {
  
   #(A): Extract band containing the calculated cloud property
    clouds <- img$select('clouds')
  
   #(B): Calculate the fraction of pixels within aoi that are marked as clouds
    CloudFraction = clouds$reduceRegion(
      reducer = ee$Reducer$mean(),
      geometry = aoi_Shapefile,
      scale = resolution_cldmsk
      )$get('clouds')
  
   #(C): Add CloudFraction as an image property to each image
    return(img$set('CloudFraction', CloudFraction))
    }

#(III): Mask all cloudy pixels in all images in the image collection.
 Add_CloudMask = function(img){
   
   #(A): Extract cloud information
   clouds <- img$select('clouds')
   
   #(B): Mask cloudy pixels with value zero
   return(img$updateMask(clouds$lt(1)))
 } 
 
#(IV): If there are empty images within the collection, then the property CloudFraction is missing for those images. We manually 
#have to add an empty CloudFraction property to those images. We therefore set the CloudFraction value to a no data value 
#of -9999 for all images where the CloudFraction value is NULL.  
 Add_NULL_CloudFraction = function(img){
   
   #(A): Set the CloudFraction value to a no data value of -9999 if equal to NULL
    CloudFraction <- ee$List(list(img$get('CloudFraction'), -9999))$reduce(ee$Reducer$firstNonNull())
    
   #(B): Return the updated clouds fraction properties  
    return(img$set("CloudFraction", CloudFraction))
    }

#(V): After merging s2col with the cloudprobability dataset some images might not have any cloud information available. For those
#images the property s2cloudless is missing. We manually have to add an empty s2cloudless property to those images. Thus, For  
#those images where the property s2_cloudless does not exist, we add this property as a constant(0) image with the bandname 
#"probability"
 Add_NULL_s2cloudless = function(img){

   #(A): Create a constant image and reproject to s2_projection
   img_constant_prob <- ee$Image$constant(0)$rename("probability")
   img_constant_prob <- img_constant_prob$reproject(s2_projection)
   
   #(B): If img does not yet have a s2cloudless property, add img_constant_prob as this property to img
   s2cloudless <- ee$List(list(img$get('s2cloudless'), img_constant_prob))$reduce(ee$Reducer$firstNonNull())

   #(C): Return the updated clouds fraction properties  
   return(img$set("s2cloudless", s2cloudless))
 }
 

#######################################################################################################################################