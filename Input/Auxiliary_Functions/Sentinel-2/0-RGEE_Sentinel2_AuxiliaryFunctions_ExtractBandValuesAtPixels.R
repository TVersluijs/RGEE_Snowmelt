####################################################################################################################################
#
#                            Sentinel 2 - Extract band values at all pixels within the area of interest
#
####################################################################################################################################

#DESCRIPTION: The function Extract_BandValuesAtPixels() can be used to iterate through all images of a Sentinel-2 image collection. 
#For each image, the value of the band of interest (e.g. NDSI) is extracted for each pixel in the image. The resulting pixel-specific
#band values are subsequently stored as feature properties in a feature collection together with each pixel's latitude, longitude
#and the datetime of the image. Within this feature collection, each pixel thus has a separate feature with its properties. This 
#results in a feature collection of all pixels for the current image. This feature collection is then appended to a list of feature
#collections from previous images in the iteration.

#The function Extract_BandValuesAtPixels() takes two arguments. The first argument is the current element of the image collection
#(i.e. the current image within the iteration) and the second element takes the output FeatureCollection from the iteration that 
#preceded it. The latter is not possible for the first iteration, that's why an initial empty feature collection to start the 
#iteration with must be defined. 

#Note that inside $map() functions all processing has to be done in the language of the server (javascript Api of google earth 
#engine). Thus, inside $map() functions, client-side functions such as $getInfo() or print(), paste0() cannot be used.

####################################################################################################################################

#(I): Function to extract BandValues at all pixels within an image and store these values in an expanding Feature Collection:
Extract_BandValuesAtPixels = function(img, FC_initial) {
    
    #for debugging
    #img <- s2_clouds_filtered$first()$select("NDSI")
    
    #(A): Calculate the Bandvalue (e.g. NDSI) at each pixel within the current image and store each pixel value as a separate feature.
    #The resulting output is a feature collection of all features (pixels) for the current image
     FC_image <- img$sample( #sampling is done automatically for all Bands of the image
        region=aoi_Shapefile, #All pixels within aoi_Shapefile will be stored as a separate feature
        geometries=TRUE,  #if TRUE, add center of sampled pixel as the geometry property of the output feature
        projection=S2Projection, #Set to native projection of S2 image
        #scale=resolution, #Set to native resolution of satellite image (implicitly set using projection above)
        seed=23, #Create reproducible results using the same random seed
        dropNulls=TRUE) #If TRUE, the result is post-filtered to drop features that have a NULL value for NDSI
    
    # #(B): Make sure there is a NDSI value at each feature within the feature collection (now redundant due to dropNulls=T above):
    # #Set the band value to a no data value of -9999 for all features where the band value is NULL.
    #  FC_image <- FC_image$map(function(feature){
    #    ndsi <- ee$List(list(feature$get('NDSI'), -9999))$reduce(ee$Reducer$firstNonNull())
    #    return(feature$set("NDSI", ndsi))})  
    
    #(C): Add date and day of year of current image as a property to each feature within the feature collection FC_image
     date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
     doy <- ee$Date(img$get('system:time_start'))$getRelative('day', 'year')
     FC_image <- FC_image$map(function(feature){return(feature$set("Date", date)$set("doy", doy))})
    
    #(D): Add latitude and longitude of each pixel as a property to each feature
     FC_image <- FC_image$map(function(feature){
        coordinates <- feature$geometry()$coordinates()
        lon <- coordinates$get(0)
        lat <- coordinates$get(1)
        return(feature$set('lon', lon)$set('lat', lat))
        })
    
    #(E): Add a pixel_ID property to each feature
     FC_image <- FC_image$map(function(feature){
        #feature=FC_image$first() #for debugging
        lon_tmp <- feature$get("lon")
        lat_tmp <- feature$get("lat")
        lon_tmp <- ee$Number$format(lon_tmp, '%.5f')
        lat_tmp <- ee$Number$format(lat_tmp, '%.5f')
        string_tmp <- ee$String("_")
        pixel_ID <- lat_tmp$cat(string_tmp)$cat(lon_tmp)
        #feature <- feature$set('pixel_ID', pixel_ID) #for debugging
        return(feature$set('pixel_ID', pixel_ID))
        })  
    
    #(F): Merge the feature collection of the current image (FC_image) onto the feature collection FC_initial.
     return (ee$FeatureCollection(FC_initial)$merge(FC_image))
     #FC_initial is updated at each iteration. This corresponds to base R code: FC_intitial <- merge(FC_initial, FC_image)
    
}

###################################################################################################################################
