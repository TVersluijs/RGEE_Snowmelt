####################################################################################################################################
#
#                                       Sentinel 2 - Calculate a relative band value
#
####################################################################################################################################

#DESCRIPTION: Create new bands within each image where the band value (NDSI, NDVI or NDMI) of each pixel is subtracted from the 
#average value for that image (NDSI_mean, NDVI_mean or NDMI_mean). This will result in a relative phenology of each pixel compared
#to the average phenology within that image at each point in time. Subtraction only works between images so we first have to 
#transform the mean band value to a constant image (but separate for each point in time).

####################################################################################################################################

#(I): Calculate average NDSI, NDVI and NDMI value for each image in the image collection:
AddMeanBandValue = function(img) {
    
    #(A): Extract bands containing the calculated cloud property
     Bands <- img$select('NDSI', 'NDVI', 'NDMI')
    
    #(B): Calculate the mean band values for each image
     BandsMean <- Bands$reduceRegion(
        reducer = ee$Reducer$mean(),
        geometry = aoi_Shapefile,
        scale = resolution
        )
    
    #(C): Set the band value to a no data value of -9999 for all features where they equal NULL
     NDSI_mean <- ee$List(list(BandsMean$get("NDSI"), -9999))$reduce(ee$Reducer$firstNonNull())
     NDVI_mean <- ee$List(list(BandsMean$get("NDVI"), -9999))$reduce(ee$Reducer$firstNonNull())
     NDMI_mean <- ee$List(list(BandsMean$get("NDMI"), -9999))$reduce(ee$Reducer$firstNonNull())
    
    #(D): Add Mean Band value as an image property to each image
     return(img$
               set('NDSI_mean', NDSI_mean)$
               set('NDVI_mean', NDVI_mean)$
               set('NDMI_mean', NDMI_mean))
     }

#(II): Relative NDSI:
get_Relative_NDSI <- function(img){
    
    #(A): Create a constant image equal to the mean band value
      mean_NDSI <- ee$Number(img$get("NDSI_mean"))
      image_mean_NDSI <- ee$Image(mean_NDSI)
    
    #(B): Add relative NDSI band to the image
     return(img$addBands(img$
                            select("NDSI")$
                            subtract(image_mean_NDSI)$
                            rename('NDSI_Relative')$
                            copyProperties(img, img$propertyNames())))
     }

#(III): Relative NDVI:
get_Relative_NDVI <- function(img){
    
    #(A): Create a constant image equal to the mean band value
     mean_NDVI <- ee$Number(img$get("NDVI_mean"))
     image_mean_NDVI <- ee$Image(mean_NDVI)
    
    #(B): Add relative NDVI band to the image 
     return(img$addBands(img$
                            select("NDVI")$
                            subtract(image_mean_NDVI)$
                            rename('NDVI_Relative')$
                            copyProperties(img, img$propertyNames())))
     } 

#(III): Relative NDMI:
get_Relative_NDMI <- function(img){
    
    #(A): Create a constant image equal to the mean band value
     mean_NDMI <- ee$Number(img$get("NDMI_mean"))
     image_mean_NDMI <- ee$Image(mean_NDMI)
    
    #(B): Add relative NDMI band to the image
     return(img$addBands(img$
                            select("NDMI")$
                            subtract(image_mean_NDMI)$
                            rename('NDMI_Relative')$
                            copyProperties(img, img$propertyNames())))
     }


###################################################################################################################################