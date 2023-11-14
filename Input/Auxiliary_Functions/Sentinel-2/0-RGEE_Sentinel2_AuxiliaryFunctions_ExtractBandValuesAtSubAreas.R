####################################################################################################################################
#
#                              Sentinel 2 - Extract band values at SubAreas within the area of interest
#
####################################################################################################################################

#DESCRIPTION: The function Extract_BandValuesAtSubAreas() can be used to iterate through all images of a Sentinel-2 image collection.
#For each image, the value of the band of interest (e.g. NDSI) is extracted for each SubArea within the feature collection aoi_SubAreas.
#The resulting SubArea-specific band values are subsequently stored as feature properties in a feature collection together with
#the datetime of the image. Within this feature collection, each SubArea thus has a separate feature with its properties. This
#results in a feature collection of all SubAreas for the current image. This feature collection is then appended to a list of feature
#collections from previous images in the iteration. Thus, at each iteration we extract band values of interest for all SubAreas within
#the current image, store this as a  feature collection and add this to an expanding list of feature collections from previous iterations.

#In addition, based on the binary band "SNOW" in each image the fraction of snow-covered pixels is calculated by dividing the total
#number of snow-covered pixels by the total number of pixels within the SubArea. This excludes any pixels that are masked due to e.g.
#cloud cover. Thus, note that this SNOW band and masking of unwanted pixels should occur before running this function.

#The function Extract_BandValuesAtSubAreas() takes two arguments. The first argument is the current element of the image collection
#(i.e. the current image within the iteration) and the second element takes the output FeatureCollection from the iteration that
#preceded it. The latter is not possible for the first iteration, that's why an initial empty feature collection to start the
#iteration with must be defined.

#Note that inside $map() functions all processing has to be done in the language of the server (javascript Api of google earth
#engine). Thus, inside $map() functions, client-side functions such as $getInfo() or print(), paste0() cannot be used.

####################################################################################################################################

#(I): Function to extract BandValues for all SubAreas within an image and store these values in an expanding Feature Collection:
Extract_BandValuesAtSubAreas = function(img, FC_initial) {

    ##For debugging
    #img <- s2_snow_masked$select("FSC_Gascoin2020", "FSC_Aalstad2020", NDSI", "NDVI", "NDMI", "SNOW")$filterDate(paste0(year_ID, "-06-05"), paste0(year_ID, "-07-31"))$first()

    #(A): Calculate the mean NDSI, NDVI, NDMI and SNOW value for the current image within each subarea in the aoi_SubAreas feature collection.
    #We add these mean values as properties to each feature (i.e. sub area). The resulting output is a feature collection for the current
    #image where these means are added as properties for each sub area (i.e. different features). Note that taking the mean of the binary
    #SNOW band results in a fraction of snow covered pixels for each image. ReduceRegions does not include masked pixels (i.e. pixels
    #defined as cloud or water) when calculating these means.
      FC_image <- img$reduceRegions(collection = aoi_SubAreas,
                        reducer = ee$Reducer$mean(), #$setOutputs(list('NDSI'))
                        scale = resolution,
                        crs=crs,
                        crsTransform=NULL)

    #(B): Add datetime of current image as a property to each feature within the feature collection FC_image
      date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
      FC_image <- FC_image$map(function(feature){return(feature$set("Date", date))})

    #THE COMMENTED CODE BELOW SHOULD BE EQUIVALENT TO CALCULATING THE MEAN NUMBER OF PIXELS DEFINED AS SNOW COVERED AS DONE ABOVE
    #
    # #(C): Add total number of pixels within each SubArea that are not masked by clouds nor water
    #
    #     #(I): Select band indicating which pixels are snow-covered (note that some of the pixels in these images might be cloud/water-masked):
    #      snow <- img$select('SNOW')
    #
    #     #(II): Count number of pixels within each subarea that are not masked by clouds nor water (i.e. both snow covered pixels and pixels without snow cover):
    #      TotalPixels <- snow$reduceRegions(
    #         reducer = ee$Reducer$count(),
    #         collection = aoi_SubAreas,
    #         scale = resolution,
    #         crs=crs,
    #         crsTransform=NULL
    #         )
    #
    #     #(III): Only select 'Cluster_ID' (SubArea) and 'count' (renamed to Totalpixels) columns:
    #      TotalPixels = TotalPixels$map(function(feature){
    #        return(ee$Feature(NULL, list(TotalPixels = feature$get('count'), Cluster_ID = feature$get('Cluster_ID'))))
    #        })
    #
    #     #(IV): Add TotalPixels as a new property to the feature collection FC_image
    #
    #       #Use an equals filter to specify how the collections match.
    #        Filter_Cluster_ID <- ee$Filter$equals(
    #          leftField= 'Cluster_ID',
    #          rightField= 'Cluster_ID')
    #
    #       #Define the join.
    #        innerJoin <- ee$Join$inner()
    #
    #       #Apply the join.
    #        Cluster_ID_Join <- innerJoin$apply(FC_image, TotalPixels,  Filter_Cluster_ID)
    #
    #       #Add features of second ('secondary') feature collection to those of the first ('primary') feature collection
    #        FC_image <- Cluster_ID_Join$map(function(pair) {
    #           f1 <- ee$Feature(pair$get('primary'))
    #           f2 <- ee$Feature(pair$get('secondary'))
    #           return(f1$set(f2$toDictionary()))
    #           })
    #
    # #(D): Add total number of snow-covered pixels within each subregion that are not masked by clouds nor water
    #
    #  #Note that all reducer functions do not include pixels where the mask is equal to 0. We have specified a cloud/water mask previously
    #  #to each image in the image collection, thus pixels containing clouds/water (with a mask equal to 0) are not used in the reduction
    #  #(i.e. in calculating the mean, i.e. fraction of snow covered pixels). We can now also mask pixels without snow, which thus
    #  #leads to images where all unmasked pixels are snow covered (and not covered by clouds/water). We can then simply count the number
    #  #of unmasked pixels to get the number of pixels that are snow covered.
    #
    #     #(I): Mask pixels without snow (i.e. SNOW=1)
    #      snow <- snow$updateMask(img$select('SNOW')$gt(0.99))
    #
    #     #(II): Count number of snow covered pixels within each subarea
    #      SnowyPixels <- snow$reduceRegions(
    #         reducer = ee$Reducer$count(),
    #         collection = aoi_SubAreas,
    #         scale = resolution,
    #         crs=crs,
    #         crsTransform=NULL
    #         )
    #
    #     #(III): Only select Cluster_ID and count (renamed to SnowyPixels) columns:
    #      SnowyPixels = SnowyPixels$map(function(feature){
    #        return(ee$Feature(NULL, list(SnowyPixels = feature$get('count'), Cluster_ID = feature$get('Cluster_ID'))))
    #        })
    #
    #     #(IV): Add Snowy pixels as a new property to the feature collection FC_image
    #
    #       #Use an equals filter to specify how the collections match.
    #        Filter_Cluster_ID <- ee$Filter$equals(
    #          leftField= 'Cluster_ID',
    #          rightField= 'Cluster_ID')
    #
    #       #Define the join.
    #        innerJoin <- ee$Join$inner()
    #
    #       #Apply the join.
    #        Cluster_ID_Join <- innerJoin$apply(FC_image, SnowyPixels,  Filter_Cluster_ID)
    #
    #       #Add features of second ('secondary') feature collection to those of the first ('primary') feature collection
    #        FC_image <- Cluster_ID_Join$map(function(pair) {
    #           f1 <- ee$Feature(pair$get('primary'))
    #           f2 <- ee$Feature(pair$get('secondary'))
    #           return(f1$set(f2$toDictionary()))
    #           })

    #(E): Make sure there is a value at each feature within the feature collection:

      #reduceRegion does not return any band values when all pixels in the sub area are masked by clouds, or do not fall within
      #the boundaries of the image. Therefore we manually have to add an empty band property to those features.
      #We therefore set the band value to a no data value of -9999 for all features where the band value is NULL.
       FC_image <- FC_image$map(function(feature){
          fsc_gascoin2020 <- ee$List(list(feature$get('FSC_Gascoin2020'), -9999))$reduce(ee$Reducer$firstNonNull())
          fsc_aalstad2020 <- ee$List(list(feature$get('FSC_Aalstad2020'), -9999))$reduce(ee$Reducer$firstNonNull())
          ndsi <- ee$List(list(feature$get('NDSI'), -9999))$reduce(ee$Reducer$firstNonNull())
          ndvi <- ee$List(list(feature$get('NDVI'), -9999))$reduce(ee$Reducer$firstNonNull())
          ndmi <- ee$List(list(feature$get('NDMI'), -9999))$reduce(ee$Reducer$firstNonNull())
          #snow <- ee$List(list(feature$get('SNOW'), -9999))$reduce(ee$Reducer$firstNonNull())
          date <- ee$List(list(feature$get('Date'), -9999))$reduce(ee$Reducer$firstNonNull())
          #totalpixels <- ee$List(list(feature$get('TotalPixels'), -9999))$reduce(ee$Reducer$firstNonNull())
          #snowypixels <- ee$List(list(feature$get('SnowyPixels'), -9999))$reduce(ee$Reducer$firstNonNull())
          return(feature$
                   set("FSC_Gascoin2020", fsc_gascoin2020)$
                   set("FSC_Aalstad2020", fsc_aalstad2020)$
                   set("NDSI", ndsi)$
                   set("NDVI", ndvi)$
                   set("NDMI", ndmi)$
                   #set("SNOW", snow)$
                   set("Date", date)#$
                   #set("TotalPixels", totalpixels)$
                   #set("SnowyPixels", snowypixels)
                 )})

    # #(F): Calculate fraction of pixels that are snow covered based on counts of pixels
    #
    #    #Map over the feature collection for this image:
    #     FC_image <- FC_image$map(function(feature){
    #       #feature <- FC_image$first() #For debugging
    #       TotalPixels <- ee$Number(feature$get('TotalPixels'))
    #       SnowyPixels <- ee$Number(feature$get('SnowyPixels'))
    #       FractionSnowyPixels <- SnowyPixels$divide(TotalPixels)
    #       FractionSnowyPixels <- ee$Number(FractionSnowyPixels) #ee$List(list(FractionSnowyPixels))
    #       return(feature$set("Snowfraction_count", FractionSnowyPixels))})

    #(G): Merge the feature collection of the current image (FC_image) onto the feature collection FC_initial.
     return (ee$FeatureCollection(FC_initial)$merge(FC_image))
     #FC_initial is thus updated at each iteration. This corresponds to base R code: FC_intitial <- merge(FC_initial, FC_image)
}

###################################################################################################################################