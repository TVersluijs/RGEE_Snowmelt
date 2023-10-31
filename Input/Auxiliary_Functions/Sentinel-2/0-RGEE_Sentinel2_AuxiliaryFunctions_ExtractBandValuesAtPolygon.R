####################################################################################################################################
#
#                              Sentinel 2 - Extract band values at area of interest
#
####################################################################################################################################

#DESCRIPTION: The function Extract_BandValuesAtPolygon() can be used to iterate through all images of a Sentinel-2 image collection.
#For each image, the average value of the band of interest (e.g. NDSI) is extracted. The resulting band values are subsequently 
#stored as feature properties in a feature collection together with the datetime of the image. This feature collection is then 
#appended to a list of feature collections from previous images in the iteration. Thus, at each iteration we extract band values 
#of interest for the polygon within the current image, store this as a  feature collection and add this to an expanding list of 
#feature collections from previous iterations.

#The function Extract_BandValuesAtPolygon() takes two arguments. The first argument is the current element of the image collection
#(i.e. the current image within the iteration) and the second element takes the output FeatureCollection from the iteration that
#preceded it. The latter is not possible for the first iteration, that's why an initial empty feature collection to start the
#iteration with must be defined.

#Note that inside $map() functions all processing has to be done in the language of the server (javascript Api of google earth
#engine). Thus, inside $map() functions, client-side functions such as $getInfo() or print(), paste0() cannot be used.

####################################################################################################################################

#(I): Function to extract BandValues for a single polygon within an image and store these values in an expanding Feature Collection:
Extract_BandValuesAtPolygon = function(img, FC_initial) {

    ##For debugging
    #img <- s2_snow_masked$select("NDSI", "NDVI", "NDMI")$filterDate(paste0(year_ID, "-06-05"), paste0(year_ID, "-07-31"))$first()

    #(A): Calculate the mean NDSI, NDVI, and NDMI value for the current image within the polygon aoi_Shapefile.
    #We add these mean values as properties to a feature. The resulting output is a feature collection for the current
    #image where these means are added as properties for aoi_Shapefile. ReduceRegions does not include masked pixels 
    #(i.e. pixels defined as cloud or water) when calculating these means.
      FC_image <- img$select("NDSI", "NDVI", "NDMI")$
          reduceRegions(collection = aoi_Shapefile,
                        reducer = ee$Reducer$mean(), #$setOutputs(list('NDSI'))
                        scale = resolution,
                        crs=crs,
                        crsTransform=NULL)

    #(B): Add datetime of current image as a property to each feature within the feature collection FC_image
      date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
      FC_image <- FC_image$map(function(feature){return(feature$set("Date", date))})

    #(C): Make sure there is a value at each feature within the feature collection:

      #reduceRegion does not return any band values when all pixels in the sub area are masked by clouds, or do not fall within
      #the boundaries of the image. Therefore we manually have to add an empty band property to those features.
      #We therefore set the band value to a no data value of -9999 for all features where the band value is NULL.
       FC_image <- FC_image$map(function(feature){
          ndsi <- ee$List(list(feature$get('NDSI'), -9999))$reduce(ee$Reducer$firstNonNull())
          ndvi <- ee$List(list(feature$get('NDVI'), -9999))$reduce(ee$Reducer$firstNonNull())
          ndmi <- ee$List(list(feature$get('NDMI'), -9999))$reduce(ee$Reducer$firstNonNull())
          date <- ee$List(list(feature$get('Date'), -9999))$reduce(ee$Reducer$firstNonNull())
          return(feature$
                   set("NDSI", ndsi)$
                   set("NDVI", ndvi)$
                   set("NDMI", ndmi)$
                   set("Date", date)
                 )})
   
    #(G): Merge the feature collection of the current image (FC_image) onto the feature collection FC_initial.
     return (ee$FeatureCollection(FC_initial)$merge(FC_image))
     #FC_initial is thus updated at each iteration. This corresponds to base R code: FC_intitial <- merge(FC_initial, FC_image)
}

###################################################################################################################################