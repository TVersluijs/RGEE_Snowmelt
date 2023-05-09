####################################################################################################################################
#
#                              Sentinel 2 - Extract band values at points within the area of interest
#
####################################################################################################################################

#DESCRIPTION: The function Extract_BandValuesAtPoints() can be used to iterate through all images of a Sentinel-2 image collection. 
#For each image, the value of the band of interest (e.g. NDSI) is extracted for each point of interest. The resulting point-specific
#band values are subsequently stored as feature properties in a feature collection together with each point's latitude, longitude
#and the datetime of the image. Within this feature collection, each point thus has a separate feature with its properties. This 
#results in a feature collection of all points for the current image. This feature collection is then appended to a list of feature
#collections from previous images in the iteration.

#The function Extract_BandValuesAtPoints() takes two arguments. The first argument is the current element of the image collection
#(i.e. the current image within the iteration) and the second element takes the output FeatureCollection from the iteration that 
#preceded it. The latter is not possible for the first iteration, that's why an initial empty feature collection to start the 
#iteration with must be defined. 

#Note that inside $map() functions all processing has to be done in the language of the server (javascript Api of google earth 
#engine). Thus, inside $map() functions, client-side functions such as $getInfo() or print(), paste0() cannot be used.

####################################################################################################################################

#(I): Function to extract BandValues at all points of interest within an image and store these values in an expanding Feature Collection:
Extract_BandValuesAtPoins = function(img, FC_initial) {
    
    #(A): Extract Band values within the current image at each point in the Locations feature collection. We add this band
    #value as a property to each feature (i.e. location). The resulting output is a feature collection for the current img
    #where the band values are added as property for each location.
     FC_image <- img$reduceRegions(collection = Locations, 
                                   reducer = ee$Reducer$mean(), #$setOutputs(list('NDSI'))
                                   scale = resolution, #implicitly specified in crs below
                                   crs=crs, #or, crs=S2Projection
                                   crsTransform=NULL)
    
    #(B): reduceRegions does not return any band-values when the point is masked by clouds or does not fall within
    #the boundaries of the image. Therefore we manually have to add an empty NDSI (or NDVI, NDMI) property to those features.
    #We therefore set the band value to a no data value of -9999 for all features where the band value is NULL.
     FC_image <- FC_image$map(function(feature){
         ndsi <- ee$List(list(feature$get('NDSI'), -9999))$reduce(ee$Reducer$firstNonNull())
         ndvi <- ee$List(list(feature$get('NDVI'), -9999))$reduce(ee$Reducer$firstNonNull())
         ndmi <- ee$List(list(feature$get('NDMI'), -9999))$reduce(ee$Reducer$firstNonNull())
         return(feature$
                   set("NDSI", ndsi)$
                   set("NDVI", ndvi)$
                   set("NDMI", ndmi))
         })
    
    #(C): Add date and day of year of current image as a property to each feature within the feature collection FC_image
     date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
     doy <- ee$Date(img$get('system:time_start'))$getRelative('day', 'year')
     FC_image <- FC_image$map(function(feature){return(feature$set("Date", date)$set("doy", doy))})
    
    #(D): Merge the feature collection of the current image (FC_image) onto the feature collection FC_initial.
     return (ee$FeatureCollection(FC_initial)$merge(FC_image))
     #FC_initial is thus updated at each iteration. This corresponds to base R code: FC_intitial <- merge(FC_initial, FC_image)
}

###################################################################################################################################
