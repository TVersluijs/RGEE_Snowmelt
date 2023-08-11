#####################################################################################################################################

#In this script the timing of snow melt is calculated based on Sentinel-2 data for all locations specified in an input file. The
#user can specify a bufferzone (radius) to depict the area in which snow melt will be analysed per location. All locations are analysed 
#consecutively (using a loop). First, a location specific bounding box is drawn per point location (taking into account the specified 
#buffer zone) and Sentinel-2 satellite data is extracted within this bounding box. Second, clouds and permanent water bodies are
#filtered within this bounding box. Third, if this bounding box overlaps with multiple satellite tiles for a certain day, a composite
#image is created (picking the pixel with least cloud cover). Finally, snow melt is analysed within each locations's buffer zone based 
#on one of the following methods (specified by the user by setting the parameter 'method'):

# (I): 'snowfraction': Calculates the fraction of pixels within each buffer zone over time where NDSI > 'NDSI_threshold', fits a GAM
#       through these data and extract the moment when this model passes a user-specified 'Snowfraction_threshold'.
# (II): 'avg_NDSI': calculates the average NDSI value over time within each point's buffer zone, fits a GAM through these data and 
#       calculates when this model passes the specified NDSI threshold defining the moment of snow melt. In addition, time series
#       of the average NDVI and NDMI are extracted within each point's buffer zone.

#The 'snowfraction' method is preferred, because it intuitively makes sense to look a the fraction of snow-covered pixels
#over time. It is harder to justify the avg_NDSI method, because it is rather unclear what this average NDSI value entails.

#This script is similar to the script '2-RGEE_TomVersluijs_Points.R'. However, in the latter script all points are analysed
#simultaneously. This has the restriction that it only works for rather small areas (<100km2) and that the user must load a 
#shapefile to specify the outline of the study area in which all points of interest should be located. This works well when looking 
#at a small area like Zackenberg. However, when points of interest are spaced further apart (like tracking data of migratory birds) 
#the shapefile required to cover all these points is so large that this likely results in computation errors. The current script
#circumvents this issue because no shapefile is required as input. The downside of the current script is that it might take 
#significantly longer to run than script '2-RGEE_TomVersluijs_Points.R'.

#Copyright Tom Versluijs 2023-07-25. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R". 
#Note that a GoogleDrive is required.

#####################################################################################################################################

      #Clear workspace and set python environment
       rm(list=ls())
       library(here)
       if(file.exists(paste0(here::here(), "/Input/rgee_environment_dir.rds"))){
         rgee_environment_dir <- readRDS(paste0(here::here(), "/Input/rgee_environment_dir.rds"))
         reticulate::use_python(rgee_environment_dir, required=T)
         reticulate::py_config()
         }
       
      #Load packages
       #renv::restore() #revert to last version of R-packages used to successfully run this script (optional).
       library(pacman)
       p_load(sf, rgee, ggplot2, mgcv, dplyr)

      #Define ggplot2 plotting theme
       theme_tom <- function(){
         theme_classic() %+replace%
           theme(axis.title = element_text(size=18),
                 axis.text = element_text(size=16),
                 legend.position = "none",
                 strip.background = element_rect(fill = "white", colour = "black", size = rel(2)),
                 complete = TRUE)}    

      #Load auxiliary functions
       source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "Sentinel2_AuxiliaryFunctions")
       sapply(source_files, source, chdir = TRUE) ; rm(source_files)

      #Initialize earth engine
       rgee::ee_Initialize(user = "tom.versluijs@gmail.com", drive = TRUE)

##################################################################################################################################
       
#II: Specify parameters of interest
       
##################################################################################################################################       
      
 #Manually specify parameters of interest

   #(a): Sentinel-2 satellite

     #Sentinel-2 dataset ("S2_HARMONIZED" or "S2_SR_HARMONIZED")
     s2_dataset <- "S2_SR_HARMONIZED"
     #Use S2_HARMONIZED (level 1c product) for data from 2016 - 2018
     #Use S2_SR_HARMONIZED (level 2a product) for data starting from 2019

     #Spatial resolution of satellite image
     #Note that this parameter only affects the calculation of the fraction of snowcover using the function
     #'AddSnowFraction' and the calculation of average bandvalues using the function 'Extract_BandValuesAtPoins'.
     #Detection of cloud-pixels occurs at resolution 'resolution_cldmsk' specified above. Detection of water-pixels
     #occurs at the native resolution of the sentinel-2 dataset (10m).
     resolution=10 # default maximum resolution for Sentinel-2 = 10m

   #(b) Area of interest

     #Specify name of study area (used as prefix in output files)
     area_name="ZAC"

     #Coordinate reference system used for calculations
     #EPSG:4326 is recommended for areas spanning multiple UTM zones, but increased computation time (i.e. spherical coordinate system).
     #EPSG:326XX is results in reduced computation time for areas located within a single UTM zone (i.e. planar coordinate system).
     crs <- "EPSG:4326"

   #(c) Point locations

     #Name of file with Locations of interest
     input_locations <- "TestLocations_Zackenberg.csv"
     #Make sure it has the columns "LON_x", "LAT_y" and "DateTime"
     #DateTime should be in format 'dd/mm/yyyy hh:mm:ss'

     #Buffer radius around each point location (in meters)
     Buffer_radius_m=250

   #(c) Dates

     #Specify the year of interest:
     year_ID <- "2022"

     #Date range of all images considered for analysis
     start_date <- paste0(year_ID, "-03-15") #choose date (well) before the first snowmelt occurs within the study site
     end_date <- paste0(year_ID, "-09-15") #choose date (well) after last date of tracking

   #(d) Snow detection

     #NDSI threshold for snow detection
     NDSI_threshold=0.42

     #Define the snowcover fraction within the aoi for which the date of its occurrence will be calculated
     Snowfraction_threshold=0.5

     #Define the preferred method for the analysis of snowmelt
     method=c("snowfraction", "avg_NDSI") #either "snowfraction", "avg_NDSI", or both using c()
     #snowfraction: Counts the fraction of pixels within the aoi with NDSI > 'NDSI_threshold' over time and extracts the moment when this
     #              fraction is equal to 'Snowfraction_threshold'.
     #avg_NDSI:     Calculates the average NDSI, NDVI and NDMI values within the aoi over time and extracts the moment when the NDSI value is equal to 'NDSI_threshold'.

   #(e): Cloud masking

     #Should clouds be masked from the analysis (default=TRUE).
     #Note that cloud masking will automatically be set to TRUE when create_composite=TRUE and/or when mask_water=TRUE
     #Highly recommended to set to TRUE.
     mask_clouds=TRUE

     #Cloud probability threshold (ranging 0-100) for S2_CLOUD_PROBABILITY dataset
     #Defines above which value of the cloud probability parameter a pixel is considered cloudy. Higher values
     #are thus more conservative (i.e. pixels are less likely to be classified as cloudy)
     cldprb_threshold=70

     #Define radius of cloud erosion and dilation (in meters)
     #First, all cloudy pixels without cloudy pixels within radius cld_erosion will be removed.
     #Second, all remaining pixels will be dilated by radius cld_buffer. Make sure that cld_buffer
     #is larger than cld_erosion. Larger values result in a more course cloud-filtering process.
     cld_erosion=250 #250 for conservative cloud estimation
     cld_buffer=270 #270 for conservative cloud estimation

     #Define resolution of cloud-filtering (in meters)
     resolution_cldmsk=60 #(default = 60)

     #Maximum fraction of cloud cover allowed in each image
     max_cloud_fraction=0.75 #(default = 0.75)

   #(f): Water masking

     #Should permanent waterbodies be masked from the analysis (default=TRUE).
     mask_water=TRUE
     mask_water_type="water_mask_ESA" #either "water_mask_ESA", "water_mask_Manual", or "both"
     #"water_mask_ESA" uses the world wide ESA dataset with a 10m resolution which works well with large areas (default)
     #"water_mask_Manual" uses a manual approach based on NDWI, NDSI and NIR bands and generally works well for small details.
     #"both" employs both methods and sets pixels to water if one or both of these methods indicates so.

     #Parameters for type "water_mask_Manual"

       #NDWI threshold above which is pixel is perceived as water (or snow!)
       NDWI_threshold=0.0

       #NDSI and NIR thresholds above which a pixel is perceived as snow
       NIR_threshold=1100

       #Date range to detect permanent waterbodies based on NDWI
       start_date_NDWI <- paste0(year_ID, "-07-15")
       end_date_NDWI <- paste0(year_ID, "-08-15")

       #Date range to detect permanent waterbodies based on NDSI (used to differentiate between water and snow in the detection of
       #permanent waterbodies).This date range should occur c.a. 2 weeks before the date range from start_date_NDWI to end_date_NDWI.
       start_date_NDSI <- paste0(year_ID, "-07-01")
       end_date_NDSI <- paste0(year_ID, "-08-01")

   #(g): Composite image

     #Specify whether a composite image will be generated for all images (i.e. tiles) per unique doy (default=TRUE). Note that setting this
     #to TRUE might result in an increase in computation time.
     create_composite=TRUE

   #(h): GAM sequential outlier filtering

     #parameters for auxiliary function to filter outliers from GAM using a two-step approach
     outlier_removal=TRUE #should sequential outlier removal be employed when fitting GAMs to the data
     outlier_thresh_1=0.4 #first threshold in outlier removal (relative to y-range of data, default=0.4)
     outlier_thresh_2=0.3 #second threshold in outlier removal (relative to y-range of data, default=0.3)
     
     #Specify the degree of GAM smoothing by setting the 'k' parameter
     gam_k_outlier=10 #Number of knots when filtering outliers (default=10)
     gam_k=50 #Number of knots when making model predictions (default=50). 
     #Larger values result in a more precise GAM-fit, but might result in overfitting.

     
##################################################################################################################################
          
#III: Define some additional parameters (automated)
          
##################################################################################################################################
     
  #Automatically define some additional parameters 
        
   #Create a unique data_ID
    data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_S2")           
  
   #Set mask_clouds to TRUE if a composite image needs to be generated         
    if(create_composite==TRUE){mask_clouds=TRUE}      
   
   #Set mask_clouds to TRUE if water masking is TRUE and mask_water_type==water_mask_Manual
    if(mask_water==TRUE & (mask_water_type=="water_mask_Manual" | mask_water_type=="both")){mask_clouds=TRUE}      
   
   #Create output folder
    if(dir.exists(paste0(here(), "/Output/S2/Points_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/S2/Points_Snowmelt"), recursive = TRUE)}
     
          
##################################################################################################################################
 
#IV: Read dataframe with points of interest
                    
##################################################################################################################################
          
  #Read a dataframe with points for which Sentinel-2 data needs to be extracted
   df_locations <- read.csv(paste0("Input/", input_locations), header=T)[,c("LON_x", "LAT_y")]
   colnames(df_locations) <- c("LON_x", "LAT_y")
          
  #Extract only unique locations and add a LocationID for every unique lat/lon combination
   df_locations <- unique(df_locations[,c("LON_x", "LAT_y")])
   df_locations$LocationID <- paste0("Location_", 1:nrow(df_locations))
          
  #plot all locations including bufferzone on a map
   Locations_ee  <-  st_as_sf(df_locations, coords = c("LON_x", "LAT_y"), crs = crs)
   Locations_ee <- st_transform(Locations_ee, crs=crs)
   Locations_ee <- sf_as_ee(Locations_ee[,c("LocationID", "geometry")])
   coordinates_point <- Locations_ee$getInfo()$features[[1]]$geometry$coordinates
   if(Buffer_radius_m>0){
         bufferBy <- function(Buffer_radius_m) {
           return(function(feature) {
             return(feature$buffer(Buffer_radius_m))   
           })
         }
         Locations_ee <- Locations_ee$map(bufferBy(Buffer_radius_m))
       }
   Map$setCenter(coordinates_point[1], coordinates_point[2], 8)
   Map$addLayer(Locations_ee, list(color="red"), paste0("Locations_", year_ID))
   rm(Locations_ee, coordinates_point)
    
       
################################################################################################################################################################################           
       
#V: Iterate through all locations of interest and extract either (I) the average NDSI, NDVI and NDMI or (II) the fraction of snowcover within each point's bufferzone over time
       
################################################################################################################################################################################          
          
 #Create an empty dataframe for storing the mean Bandvalue results per location
  Locations_BandValues <- data.frame(LocationID=character(),
                                     Date=character(),
                                     doy=numeric(),
                                     NDSI=numeric(),
                                     NDVI=numeric(),
                                     NDMI=numeric(),
                                     LON_x=numeric(),
                                     LAT_y=numeric())
       
 #Create an empty dataframe for storing the Snowfraction results per location 
  aoi_SnowCover <- data.frame(LocationID=character(),
                              Date=character(),
                              doy=numeric(),
                              SnowFraction=numeric(),
                              LON_x=numeric(),
                              LAT_y=numeric())
   
 #Loop through all locations and extract the average band value and/or snowfraction within the location including buffer zone
  for(i in 1:nrow(df_locations)){
    
  #(A): Print progress:
     print(paste0("Starting analysis for location ", i, " out of ", nrow(df_locations), " locations"))
   
  #(B): Load current point of interest and add required buffer zone:
    
     #Store current location in separate dataframe and store its locationID
      Locations <- df_locations[i,]
      Location_ID <- Locations[,'LocationID']
    
     #Transform current location to an ee_object
      Locations  <-  st_as_sf(Locations, coords = c("LON_x", "LAT_y"), crs = crs)
      Locations <- st_transform(Locations, crs=crs)
      Locations <- sf_as_ee(Locations[,c("LocationID", "geometry")])
    
     #Store the central point as an ee_Geometry object
      coordinates_point <- Locations$getInfo()$features[[1]]$geometry$coordinates
      point <- ee$Geometry$Point(coordinates_point[1], coordinates_point[2])
    
      #Locations is now a feature collection with only one feature (point). We are now going to add a buffer around this single feature
      #which will transform it into a polygon with radius "Buffer_radius_m" and as a center the original point location. This feature 
      #has as property the LocationID and as geometry the outline of the polygon.
    
     #Define the buffer function
      bufferBy <- function(Buffer_radius_m){
        return(function(feature){return(feature$buffer(Buffer_radius_m))})
        }
    
     #Define bounding box function
      bounding_box_func = function(feature){
        bounding_box = feature$bounds()
        return(bounding_box)
        }
    
     #Add buffer around the point in Locations, and define aoi (i.e. bounding box around buffer zone for Locations)
      if(Buffer_radius_m>0){
      
        #Add buffer radius around point location (creating a polygon)
         Locations <- Locations$map(bufferBy(Buffer_radius_m))
      
        #Define area of interest (aoi) by drawing a bounding box around the polygon 'Locations'
         aoi <- Locations$map(bounding_box_func)
         }
      if(!Buffer_radius_m>0){
      
        #Add buffer radius around point location (creating a polygon)
         Locations <- Locations$map(bufferBy(resolution))
      
        #Define area of interest (aoi) by drawing a bounding box around the polygon 'Locations'
         aoi <- Locations$map(bounding_box_func)
         }
    
     # #Plot feature on a map (i.e. polygon of bufferzone surrounding Locations)
     #  Map$setCenter(coordinates_point[[1]], coordinates_point[[2]], 10)
     #  Map$addLayer(point, list(color="red"))+
     #  Map$addLayer(Locations, list(color="green"))+
     #  Map$addLayer(aoi$style(color="blue", fillColor='00000000'))
    
     #Transform aoi to a geometry object  
      aoi <- aoi$geometry()
    
     #Create aoi_shapefile (corresponding to buffer zone of Location i) as this is required input for several auxiliary functions.  
      aoi_Shapefile <- ee$FeatureCollection(Locations)
    
 #(C): Extract sentinel-2 data for the point of interest
    
    #Extract Sentinel-2 Surface Reflectance data for a given area and date range:
     s2_col <- ee$ImageCollection(paste0('COPERNICUS/', s2_dataset))
     s2_col <- s2_col$
       filterBounds(aoi)$
       filterDate(start_date, end_date)
    
    #Add Sentinel-2 cloud probability as a new property to s2_col
     if(mask_clouds==TRUE){
     
       #Load cloud information from S2_CLOUD_PROBABILITY
       s2_cloudless_col <-ee$ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
       s2_cloudless_col <- s2_cloudless_col$
         filterBounds(aoi)$
         filterDate(start_date, end_date) 
     
       #Store sentinel-2 projection (of band "probability")
       s2_projection <- s2_cloudless_col$first()$select("probability")$projection() 
     
       #Add Sentinel-2 cloud probability as a new property to s2_col
       s2_col <- Add_CloudProb_Property(s2_col=s2_col, s2_cloudless_col=s2_cloudless_col)
      
      #Make sure that each image contains the s2cloudless property (create an empty property if it does not)
       s2_col <- s2_col$map(Add_NULL_s2cloudless)
     
      # #Check that s2cloudless was added as a property (for debugging)
      #  s2_col$first()$propertyNames()$getInfo()  
     
      }
   
    # #Clip all images to the shapefile area (aoi_Shapefile). This is not necessary because the auxilliary functions we use
    # #later in this script make calculations specifically for aoi_Shapefile (buffer zone) and not for aoi (bounding box).
    #  s2_col <- s2_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})    
 
    #Add a NDSI, NDVI, NDMI and NDWI band to the clipped image collection
     s2_col <- s2_col$
       map(getNDSI)$
       map(getNDVI)$
       map(getNDMI)$
       map(getNDWI)      
    
    # #Plot the first image in s2_col (for Debugging):
    #  image <- s2_col$filterDate(paste0(year_ID, "-06-17"), end_date)$first()
    # 
    # #Plot an RGB, NDSI and NDVI, NDMI and NDWI image for the first extracted satellite image
    # Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
    # Map$addLayer(image, list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
    #   Map$addLayer(image, list(bands=c("NDSI"), min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
    #   Map$addLayer(image, list(bands=c("NDVI"),min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
    #   Map$addLayer(image, list(bands=c("NDMI"),min=-1, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDMI')+
    #   Map$addLayer(image, list(bands=c("NDWI"),min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')
    
    #Create a timeseries gif of RGB images for the area of interest (For debugging)
    
      # #Check number of images in collection
      # s2_col$size()$getInfo()

      # #Create a timelapse video of RGB band
      # videoArgs <- list(dimensions=200, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
      # tryCatch({browseURL(s2_col$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
    
 #(D): Filter and mask clouds within the image collection
    
    #There are many images that are almost fully covered in dense clouds. We thus first filter and mask clouds in the image collection.
    #We delete all images with a cloud cover fraction >= max_cloud_fraction. In the remaining images we mask all pixels that are
    #covered by clouds. 
     if(mask_clouds==TRUE){ 
    
       #print message   
       print("Cloud masking = TRUE")
      
       #Add the cloud fraction within the buffer zone surrounding Location to each separate image by mapping the cloud functions over the image collection
       s2_col <- s2_col$
         #Determine which pixels are clouds using the cloud_probability property
         map(compute_Clouds_cldprb)$
         #Add the fraction of cloud-covered pixels within the shapefile area as image property
         map(Add_CloudFraction)$
         #Add NULL to those images in which cloudfraction could not be calculated
         map(Add_NULL_CloudFraction)$
         #Add date and time characteristics to each image
         map(add_Date)$
         #Add 'seconds since 1st of january' as a new band to each image to make sure that each pixel in the composite image can be traced back to the original image.
         map(add_Time)
    
        # #Check if CloudFraction property has been added to each image (for debugging)
        # s2_col$first()$propertyNames()$getInfo()
    
        # #Visual check of cloud mask
        #  image <- s2_col$filterDate(paste0(year_ID, "-05-31"), end_date)$first()$
        #    clipToCollection(aoi_Shapefile)$
        #    select("clouds_unadjusted", "clouds", "clouds_inv", "clouds_prob", "B4", "B3", "B2")      
        #  Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
        #  Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
        #    #Map$addLayer(image, list(bands='clouds_prob', min=0, max=100, opacity=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'clouds_prob')+
        #    Map$addLayer(image, list(bands='clouds_unadjusted', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'clouds unadjusted')+
        #    Map$addLayer(image, list(bands='clouds', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'clouds buffered')
        #  #Map$addLayer(image, list(bands='clouds_inv', min=0, max=1, opacity=1), 'clouds_inv')
        #  
        #We have now calculated the fraction of cloud covered pixels within our buffer zone for each image in the image collection.
        #However, Sentinel-2 satellite data is recorded as different tiles, these tiles are generally much bigger than a single study site.
        #Each Sentinel-2 tile also contains a native measure of cloud cover. This is the property 'CLOUDY_PIXEL_PERCENTAGE'. This property
        #thus refers to the percentage of cloud cover within the whole tile. We have now calculated a more appropriate measure for our
        #area of interest only called 'CloudFraction'. Below we compare both measures, as we expect they should be correlated. Differences
        #between both measures can for example arise when the area of interest was covered by clouds while most of the larger region (tile)
        #was not or vice versa.
    
        # #Compare calculated cloud fraction to that as calculated for the complete tile
        # #Extract cloud fraction of all images in image collection for the area of interest
        # cloud_aoi <- unlist(s2_col$aggregate_array('CloudFraction')$getInfo())
        # 
        # #Extract cloud fraction of all images in image collection for the complete tile to which our aoi belongs
        # cloud_tile <- unlist(s2_col$aggregate_array('CLOUDY_PIXEL_PERCENTAGE')$getInfo())
        # 
        # #Replace -9999 values by NA
        # cloud_aoi[cloud_aoi < -9000] <- NA
        # cloud_tile[cloud_tile < -9000] <- NA
        # 
        # #Assess relationship between both measures:
        # ggplot()+
        #   geom_point(aes(x=cloud_aoi, y=cloud_tile))+
        #   stat_smooth(method='lm', formula=y~x)+
        #   theme_tom()
        # #Thus, this indicates that our local cloud measure does not necessarily correspond to that as measured on the whole tile as
        # #explained above.
    
       #Exclude all images from the image collection that have a CloudFraction value > max_cloud_fraction, and mask all cloud pixels within the remaining images
        s2_clouds_filtered <- s2_col$
          #Only remove images that are fully cloud-covered
          filter(ee$Filter$lt('CloudFraction', max_cloud_fraction))$
          #Apply cloudmask for individual pixels
          map(Add_CloudMask)
    
       # #Create timelapse video of the cloud filtered/masked RGB images (for debugging)
       # videoArgs <- list(dimensions=200, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
       # tryCatch({browseURL(s2_clouds_filtered$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
    
       }
     if(mask_clouds==FALSE){
      
       #print message   
        print("Cloud masking = FALSE")
      
       #Add Date and Time characteristics to each image
       s2_clouds_filtered <- s2_col$
         map(add_Date)$
         #Add 'seconds since 1st of january' as a new band to each image to make sure that each pixel in the composite image can be traced back to the original image.
         map(add_Time)
      
       }
  
 #(E): Mask permanent waterbodies (ponds, lakes, rivers, sea) within the image collection
    
    #Mask permanent waterbodies if mask_water==TRUE
     if(mask_water==TRUE){
      
       #(A): print message   
       print("Water masking = TRUE")
       
       #(B): Method 1: Extract permanent waterbodies from the ESA WorldCover map
        if(mask_water_type=="water_mask_ESA" | mask_water_type=="both"){
           
           #print message   
           print("Water masking type = ESA")
           
           #Load auxilliary function
           compute_Water_ESA=compute_Water_ESA #sourced
           
           #Extract pixels corresponding to permanent waterbodies 
           water_mask1 <- compute_Water_ESA(ESA_index=80)
           
         }
       
       #(C): Method 2: Manual detection of permanent waterbodies based on NDWI, NDSI and NIR bands
        if(mask_water_type=="water_mask_Manual" | mask_water_type=="both"){
       
           #print message   
           print("Water masking type = Manual")
         
           #Load auxilliary function
           compute_Water_Manual=compute_Water_Manual #sourced
         
           #Extract pixels corresponding to permanent waterbodies
           water_mask2 <- compute_Water_Manual(img_col=s2_col, start_date_NDWI=start_date_NDWI, end_date_NDWI=end_date_NDWI, NDWI_threshold=NDWI_threshold,
                                               start_date_NDSI=start_date_NDSI, end_date_NDSI=end_date_NDSI, NDSI_threshold=NDSI_threshold, NIR_threshold=NIR_threshold)

           
        }
           
       #(D): Apply the water masking function to each image in the collection:
        
         #Define the final water_mask:
         if(mask_water_type=="water_mask_ESA"){water_mask <- water_mask1}
         if(mask_water_type=="water_mask_Manual"){water_mask <- water_mask2}
         if(mask_water_type=="both"){water_mask <- water_mask1$Or(water_mask2)}     
       
         # #Plot permanent waterbodies(for debugging)
         # Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
         # Map$addLayer(s2_col$filterDate(paste0(year_ID, "-08-04"), end_date)$first(), list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
         # Map$addLayer(water_mask,list(min=0, max = 1, palette = c('ffffff', '000000')), 'Water_mask')
         
         #Load auxilliary function to mask water-pixels
         Add_WaterMask=Add_WaterMask #sourced
       
         #Apply the final watermask:
         s2_clouds_filtered <- s2_clouds_filtered$map(Add_WaterMask)
      
       # #Create a timeseries GIF of RGB images of the water and cloud filtered image collection (for debugging)
       # videoArgs <- list(dimensions=200, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
       # tryCatch({browseURL(s2_clouds_filtered$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
      
      }
     if(mask_water==FALSE){
       
       #print message   
       print("Water masking = FALSE")
       
     }
     
 #(F): Create a composite image for each day by mosaicking all images from the that day
    
    #There might be days for which multiple satellite photos are available. We thus need to make a decision on what to do with overlapping 
    #values per pixel. We use a qualityMosaic that can select the overlapping pixel with the least cloudcover by tapping into the band 
    #"cloud_inv". If two pixels have equal cloudcover, then the pixelvalue from this first image is chosen.
    
    #Note that cloud/water masking is taken into account when making the composite image. The reducer function "qualitymosaic" will  for every pixel only use
    #the unmasked pixels to calcualte the output value (i.e. if for a certain pixel there are 4 overlapping images available for the same date and one
    #of them is covered by clouds (and thus masked), then only the remaining cloud-free (and unmasked) pixels are used to calculate the band values
    #at that pixel). See, https://gis.stackexchange.com/questions/337666/applying-google-earth-engine-cloud-mask-on-each-individual-image-before-mosaicin,
    #or: https://gis.stackexchange.com/questions/371356/combining-single-cloud-masked-images-using-google-earth-engine, or
    #https://open-mrv.readthedocs.io/en/latest/image_composite_Web.html.
    
   #F.1: Make a composite image for a single date (doy) - USED FOR DEBUGGING ONLY
    {
    # #Pick a date and check whether there are multiple tiles spanning the aoi
    # doy=145
    # s2_img <- s2_clouds_filtered$filterMetadata('year', 'equals', as.numeric(year_ID))$filterMetadata('doy', 'equals', doy)
    # 
    # #Number of images in s2_clouds_filtered
    # s2_img$size()$getInfo()
    # 
    # #Inspect whether there is any overlap between the tiles for this day:
    # Map$setCenter(coordinates_point[1], coordinates_point[2], 8)
    # Map$addLayers(s2_img, list(bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7)), 'S2_Original')+
    # #Map$addLayers(s2_img, list(bands=c("seconds"), min=15053255, max=15053280), opacity=0.75, 'image_ID')
    # Map$addLayer(aoi)
    # 
    # #Create a composite image for this day
    # s2_composite <- ee$Image(s2_img$qualityMosaic('clouds_inv'))#$reproject(crs=crs, crsTransform=NULL, scale=resolution)
    # s2_composite <- s2_composite$clip(aoi)
    # #s2_composite$select('NDSI')$projection()$getInfo()
    # #Note that there is no need to reproject as this is done automatically when plotting the data, or when extracting other output.
    # 
    # #Inspect the RGB band of the Mosaicked image
    # Map$setCenter(coordinates_point[1], coordinates_point[2], 8)
    # Map$addLayers(s2_img, list(bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7)), 'S2_Tiles')+
    # Map$addLayer(s2_composite, list(bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7)), 'S2_Median')
    # 
    # #Inspect the "clouds" band of the Mosaicked image
    # Map$setCenter(coordinates_point[1], coordinates_point[2], 8)
    # Map$addLayers(s2_img$select("clouds_inv"), list(palette=c("red", "green"), min=0, max=1), opacity=1, 'clouds_Tiles')+
    # Map$addLayer(s2_composite, list(bands=c("clouds_inv"), palette=c("red", "green"), min=0, max=1), opacity=1, 'Clouds_median')
    # #note that this image is only shown properly when fully zoomed in
    # 
    # # #Inspect to which original image each pixel in the Mosaicked image belongs (using the "seconds" band)
    # # Map$addLayer(s2_composite, list(bands=c("seconds"), min=15053255, max=15053280), opacity=1, 'image_ID')
    }
      
   #F.2: Create a median composite image for each day of year (doy) per year
   #Create a list of unique doy per year, loop through this list and add all images belonging to each doy to this list, group
   #images and create a composite image by reducing each pixel value using a qualitymosaic based on cloudcover.
   if(create_composite==TRUE){   
     
    #print message   
     print("Composite image = TRUE")
     
    #Define a list of unique years and doy from the image collection.
     years <- ee$List(s2_clouds_filtered$aggregate_array('year'))$distinct()$sort()
     days <- ee$List(s2_clouds_filtered$aggregate_array('doy'))$distinct()$sort()
    
    #Map over a the list of each distinct doy, and year to build a list of image mosaics
     s2_colList <- years$map(ee_utils_pyfunc(function(year){
      
      return(days$map(ee_utils_pyfunc(function(doy){
        
        #Select all images for current doy and year combination
        selection <- s2_clouds_filtered$
          #Filter image collection by year
          filterMetadata('year', 'equals', year)$
          #Filter image collection by doy
          filterMetadata('doy', 'equals', doy)$
          #Select NDSI band only
          select(c('B4', 'B3', 'B2', 'NDSI', 'NDVI', 'NDMI', 'NDWI', 'clouds', 'clouds_inv'))

        #Create a composite image for the current doy and year combination
        selection_mosaic <- ee$Image(selection$
                                       #Create Quality mosaic
                                       qualityMosaic('clouds_inv')$ #median()$ #mosaic() #reduce(ee$Reducer$qualityMosaic('clouds_inv'))
                                       #Reproject mosaicked image
                                       #reproject(crs=crs, crsTransform=NULL, scale=resolution)$
                                       #Copy timestamp properties of first image to mosaic image
                                       copyProperties(selection$first(), c('system:id', 'system:time_start')))
        
        #Return composite median image (which is stored in list)  
        return(selection_mosaic$
                 #Set doy as an image property
                 set('doy', doy)$
                 #Set year as an image property
                 set('year', year)
               )})))}))$flatten()
    
    #Convert the image List to an ImageCollection.
     s2_col_composite <- ee$ImageCollection$fromImages(s2_colList)
    
    #Inspect the new composite image collection
     #s2_clouds_filtered$size()$getInfo() #before compositing
     #s2_col_composite$size()$getInfo() #after compositing
     #s2_col_composite$first()$bandNames()$getInfo()
     #s2_col_composite$first()$propertyNames()$getInfo()
     #s2_col_composite$aggregate_stats("system:time_start")$getInfo()
     #s2_col_composite$first()$select('NDSI')$projection()$getInfo()
    
     #Note that the default projection of a composite image is "EPSG:4326" with a scale of 1 degree. This is NOT a  
     #problem as we specify  the required output projection at all reduceRegions and export steps below. See:
     #https://developers.google.com/earth-engine/guides/projections#the-default-projection
     #Alternatively, we can specify the reproject function several lines above (which is now commented)
    
   }
   if(create_composite==FALSE){
     
     #print message   
     print("Composite image = FALSE")
     
     #set s2_col_composite to s2_clouds_filtered (i.e. to unmosaicked image)
     s2_col_composite <- s2_clouds_filtered
     
   }
       
    # #Create a timelapse video of RGB band
    # videoArgs <- list(dimensions=400, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
    # tryCatch({browseURL(s2_col_composite$getVideoThumbURL(videoArgs)) }, error = function(cond){return("Too many pixels. Reduce dimensions.")})
    # #Note that missing pixels are actually not recorded by the satellite and are NOT caused by coding errors    

  #(G): Extract the fraction of snow covered pixels within the buffer zone of 'Location' for all images in the image collection
    
    if("snowfraction" %in% method){ 
     
    #(G.1): Add a binary snow cover band to each image and calculate the fraction of snow covered pixels within the buffer zone (aoi_shapefile)
     
     #Map Snow computation functions over the mosaicked image collection
     s2_col_composite <- s2_col_composite$
       #Determine which pixels are snow-covered (NDSI > NDSI threshold)
       map(computeSnow)$
       #add the fraction of snow covered pixels within aoi_Shapefile as an image property (excluding cloud masked pixels from the calculations)
       map(AddSnowFraction)$
       #Add a NULL value to images for which the snow fraction could not be calculated
       map(AddNULLSNOW)
     
       # #Mask and display the binary layer (for debugging).
       #  Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
       #  Map$addLayer(s2_col_composite$first(),list(bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
       #  Map$addLayer(s2_col_composite$select('SNOW')$first())

    #(G.2): Extract fraction of snowcover within the buffer zone for each image in the image collection   
     
     #Create an empty feature collection:
      FC_initial <- ee$FeatureCollection(ee$Feature(NULL))
     
     #Define an Iteration function to extract the fraction of snowcover within the buffer zone of Location for each image
      extract_snowcover <- function(img, FC_initial){

       #Extract snowfraction and doy from current image
       SnowFraction <- img$get('SnowFraction')
       date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
       
       #Store snowFraction, doy and date as properties in a featurecollection with one feature
       Feature_tmp <- ee$FeatureCollection(ee$Feature(NULL, c(SnowFraction=SnowFraction, Date=date)))
       
       #Merge the feature collection of the current image (Feature_tmp) onto the feature collection FC_initial.
       return (ee$FeatureCollection(FC_initial)$merge(Feature_tmp))
       #FC_initial is thus updated at each iteration. This corresponds to base R code: FC_intitial <- merge(FC_initial, FC_image)
     
       }
     
     #iterate this function over all images in the collection (output a feature collection)
      FC_merged <- ee$FeatureCollection(s2_col_composite$iterate(extract_snowcover, FC_initial))

     #Export the feature collection as a .csv table 
     #We export the data instead of using aggregate_array() as the latter might fail due to computation timeouts.
     
       #Setup task
        task_vector <- ee_table_to_drive(
          collection= FC_merged,
          description = paste0(data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_", Location_ID, "_Data_SnowFraction"),
          folder="RGEE_tmp",
          fileFormat="CSV",
          selectors=c('SnowFraction', 'Date')
          )
     
       #Run and monitor task
        print(paste0("calculating the fraction of snowcover for location ", i))
        task_vector$start() 
        ee_monitoring(task_vector, max_attempts=1000000) #250s at 100m
     
       #Import results
        exported_stats <- ee_drive_to_local(task=task_vector, dsn=paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_", Location_ID, "_Data_SnowFraction"))
        aoi_SnowCover_new <- read.csv(exported_stats)
        unlink(exported_stats)

       #Add day of year
        aoi_SnowCover_new$Date <- as.POSIXct(aoi_SnowCover_new$Date, format="%Y-%m-%d %H:%M:%S")
        aoi_SnowCover_new$doy <- as.numeric(strftime(aoi_SnowCover_new$Date, format = "%j"))
        
       #Remove NAs in the SnowFraction variable  
        aoi_SnowCover_new$SnowFraction[aoi_SnowCover_new$SnowFraction < -9000] <- NA #replace -9999 by NA
        aoi_SnowCover_new <- aoi_SnowCover_new[!is.na(aoi_SnowCover_new$SnowFraction), ]

       #Add LocationID
        aoi_SnowCover_new$LocationID <- Location_ID
       
       #Sort dataframe by LocationID and doy
        index <- with(aoi_SnowCover_new, order(LocationID, doy))
        aoi_SnowCover_new <- aoi_SnowCover_new[index,]
        aoi_SnowCover_new <- aoi_SnowCover_new[,c("LocationID", "Date", "doy", "SnowFraction")]
        
       #Add coordinates to each location  
        aoi_SnowCover_new <- left_join(aoi_SnowCover_new, df_locations[,c("LON_x", "LAT_y", "LocationID")], by=c("LocationID"))
        
       #Add dataframe for current location to dataframe from previous iterations:
        aoi_SnowCover <- rbind(aoi_SnowCover, aoi_SnowCover_new)
        
       # #Save dataframe for current location
       #  write.csv(aoi_SnowCover_new, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_", Location_ID, "_Data_SnowFraction.csv"), row.names = FALSE)
     
    }
 
  #(H): Extract Sentinel-2 mean band values (NDSI, NDVI, NDMI) within the buffer zone of 'Location' for all images in the image collection 
    
    if("avg_NDSI" %in% method){
     
    #Create an iteration function that we will use to iterate through all images of the image collection. For each image, the  
    #average value of certain bands is extracted within the polygon of Location (i.e. equivalent to aoi_Shapefile). The resulting
    #band values (and the datetime of the image) are added as a property to each feature (i.e. location). This results
    #in an updated feature collection specific for the current image. This feature collection is then appended to a list
    #of feature collections from previous image iterations. After iterating through n-images, the result is thus a feature
    #collection list of length n * the number of features (in this case 1 feature) within the feature collection.
    
    #Thus, at each iteration we extract the mean band values of interest within the buffer zone of Location (i.e. aoi_Shapefile) for 
    #the current image, store this as a feature collection and add this to an expanding list of feature collections from previous 
    #iterations for previous images (dates).
    
    #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
     FC_initial <- ee$FeatureCollection(ee$List(list()))
    
    #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
    #(in this case the current iteration image) and the second element takes the output value from the iteration that preceded it. The latter
    #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with 
    #has been defined. This function calculates mean NDSI, NDVI and NDMI within aoi_shapefile (buffer zone of Location i)
     Extract_BandValuesAtPoins = Extract_BandValuesAtPoins #sourced
    
    #Iterate over the ImageCollection and select appropriate band values
     FC_merged <- ee$FeatureCollection(s2_col_composite$select("NDSI", "NDVI", "NDMI")$iterate(Extract_BandValuesAtPoins, FC_initial))
    
 #(H): Transform the feature collection with average Band values per day per location to a dataframe by exporting it
    
    #(I): Google Drive (preferred option, slightly faster than using Asset folder)
    
      #Setup task
       task_vector <- ee_table_to_drive(
        collection= FC_merged,
        description = paste0(data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_", Location_ID, "_Data_MeanBandvalues"),
        folder="RGEE_tmp",
        fileFormat="CSV",
        selectors=c('LocationID', 'Date', 'doy', 'NDSI', 'NDVI', 'NDMI')
        )
    
      #Run and monitor task
       print(paste0("Calculating the mean NDSI, NDVI and NDMI for location ", i))
       task_vector$start() 
       ee_monitoring(task_vector, max_attempts=1000000) #250s at 100m
    
      #Import results
       exported_stats <- ee_drive_to_local(task=task_vector, dsn=paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_", Location_ID, "_Data_MeanBandvalues"))
       Locations_BandValues_new <- read.csv(exported_stats)
       unlink(exported_stats)
    
    # #(II): Asset folder (alternative option to Google Drive, but slightly slower)
    # 
    #   #Create an asset folder on local machine
    #    assetid <- paste0(ee_get_assethome(), '/RGEE_tmp')
    # 
    #   #Setup task
    #    task_vector <- ee_table_to_asset(
    #     collection = FC_merged,
    #     assetId = assetid,
    #     overwrite = TRUE
    #     )
    # 
    #   #Run and monitor task
    #    task_vector$start()
    #    ee_monitoring(task_vector, max_attempts=1000000)
    # 
    #   #Import results
    #    ee_fc <- ee$FeatureCollection(assetid)
    #    Locations_BandValues_new <- data.frame(LocationID=unlist(ee_fc$aggregate_array('LocationID')$getInfo()),
    #                                           Date=unlist(ee_fc$aggregate_array('Date')$getInfo()),
    #                                           doy=unlist(ee_fc$aggregate_array('doy')$getInfo()),
    #                                           NDSI=round(as.numeric(unlist(ee_fc$aggregate_array('NDSI')$getInfo())), 5), 
    #                                           NDVI=round(as.numeric(unlist(ee_fc$aggregate_array('NDVI')$getInfo())), 5), 
    #                                           NDMI=round(as.numeric(unlist(ee_fc$aggregate_array('NDMI')$getInfo())), 5)
    #                                           )
    
      #Change -9999 to NA
       Locations_BandValues_new$NDSI[Locations_BandValues_new$NDSI < -9000] <- NA
       Locations_BandValues_new$NDVI[Locations_BandValues_new$NDVI < -9000] <- NA
       Locations_BandValues_new$NDMI[Locations_BandValues_new$NDMI < -9000] <- NA
    
      #Sort dataframe by LocationID and doy
       index <- with(Locations_BandValues_new, order(LocationID, doy))
       Locations_BandValues_new <- Locations_BandValues_new[index,]
    
      #Add coordinates to each location  
       Locations_BandValues_new <- left_join(Locations_BandValues_new, df_locations[,c("LON_x", "LAT_y", "LocationID")], by=c("LocationID"))
    
      #Add dataframe for current location to dataframe from previous iterations:
       Locations_BandValues <- rbind(Locations_BandValues, Locations_BandValues_new)
      
      # #Save dataframe for current location
      #  write.csv(Locations_BandValues_new, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_", Location_ID, "_Data_MeanBandvalues.csv"), row.names = FALSE)
  
    }
       
  }
  
 #Save combined dataframe for all locations
  write.csv(Locations_BandValues, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_MeanBandvalues_polygon.csv"), row.names = FALSE)
  write.csv(aoi_SnowCover, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_SnowFraction_polygon.csv"), row.names = FALSE)    
  #Locations_BandValues <- read.csv(paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_MeanBandvalues_polygon.csv"), header = T)
  #aoi_SnowCover <- read.csv(paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_SnowFraction_polygon.csv"), header=T)
  
  
##################################################################################################################################
        
#VI: Fit a Generalized Additive Model (GAM) through the Sentinel-2 band value data for each location
        
##################################################################################################################################       
  
  #We fit a Generalized Additive Model (GAM) through the data for each point location. We do this using a sequential
  #outlier-filtering process. We first fit a GAM through the data, calculating model predictions and residuals. We then 
  #exclude all rows from the dataframe where the residual >= (0.25 * the range of the data). We then re-fit a GAM to 
  #this reduced dataset, make predictions and calculate residuals. We then exclude all rows from the reduced dataframe 
  #where the residual >= (0.1 * the range of the data). This gives us a final dataframe through which we fit a third
  #GAM and store its model predictions, given the dataset in which outliers were thus removed through two subsequent 
  #steps. This whole process is executed using the function f_gam_SeqRemOutliers. Note that a sequential step is
  #required because initially some datapoints might falsely be assigned a large residual because of one extreme
  #outlier. After removal of this extreme outlier and refitting of the GAM, it can be better assessed which data
  #points truly have a large residual and can thus be assigned as 'true' outliers.
  f_gam_SeqRemOutliers <- f_gam_SeqRemOutliers #sourced
  
  #Specify sequential outlier thresholds:
  outlier_thresh_1=outlier_thresh_1
  outlier_thresh_2=outlier_thresh_2
  outlier_removal=outlier_removal
  
###################################################################################################################################  
  
 #(SNOWFRACTION) - Fit a Generalized Additive Model (GAM) through the extracted SnowFraction data  
  
  if("snowfraction" %in% method){
  
   #(A) Specify for which snowfraction the corresponding date needs to be extracted
    Snowfraction_threshold = Snowfraction_threshold
  
   #(B) Create an empty dataframe
     Locations_SnowFraction <- data.frame(SnowFraction=numeric(),
                                          Date=factor(),
                                          doy=numeric(),
                                          LocationID=factor(),
                                          outliers=logical())
  
     Locations_SnowFraction_predictions <- data.frame(LocationID=character(), 
                                                      doy=numeric(), 
                                                      SnowFraction_gam_predict=numeric(), 
                                                      stringsAsFactors=FALSE)   
  
   #(C) Loop through all Locations and fit a separate gam with sequential outlier removal to the location specific SnowFraction data
    for(i in unique(aoi_SnowCover$LocationID)){
    
      #For debugging  
      #i=unique(aoi_SnowCover$LocationID)[1]
    
      #Select Location-specific subset of data:
       df_Location_SnowFraction <- aoi_SnowCover[aoi_SnowCover$LocationID==i & !is.na(aoi_SnowCover$SnowFraction), 
                                                 c("SnowFraction", "Date", "doy", "LocationID")]
    
      #Fit a gam through the Location-specific mean_NDSI ~ doy data and employ sequential outlier removal
       df_Location_SnowFraction <- f_gam_SeqRemOutliers(data=df_Location_SnowFraction, y="SnowFraction", x="doy", outlier_removal=outlier_removal, 
                                                        outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                        default_k=gam_k_outlier)

      #Sort df_Location_SnowFraction by doy:
       df_Location_SnowFraction <- df_Location_SnowFraction[order(df_Location_SnowFraction$doy),]
    
      #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
       Locations_SnowFraction <- rbind(Locations_SnowFraction, df_Location_SnowFraction)
    
      #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves
    
        #Refit GAM through data
         index <- which(df_Location_SnowFraction$outliers==FALSE)
         mod_gam <- with(df_Location_SnowFraction[index,], mgcv::gam(SnowFraction ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))
    
        #Use gam to make predictions on a more detailed (1-day) day of year interval
         aoi_SnowFraction_predicted <- data.frame(LocationID=i, doy=seq(min(df_Location_SnowFraction$doy), max(df_Location_SnowFraction$doy), 1))
         aoi_SnowFraction_predicted$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=aoi_SnowFraction_predicted, type="response")
         aoi_SnowFraction_predicted <- aoi_SnowFraction_predicted[order(aoi_SnowFraction_predicted$doy),]
         aoi_SnowFraction_predicted$year <- year_ID
   
      #Add predictions to Locations_SnowFraction_predictions dataframe:  
       Locations_SnowFraction_predictions <- rbind(Locations_SnowFraction_predictions, aoi_SnowFraction_predicted)
    
    }
  
      #Change column LocationID to a factor:
       Locations_SnowFraction$LocationID <- as.factor(as.character(Locations_SnowFraction$LocationID))
       Locations_SnowFraction_predictions$LocationID <- as.factor(as.character(Locations_SnowFraction_predictions$LocationID))
  
   #(D) Plot the raw SnowFraction datapoints and gam predictions for each Location:
  
    # #Plot SnowFraction and model predictions for all locations in a single plot
    #   p_Locations_SnowFraction = ggplot()+ 
    #    geom_point(data=Locations_SnowFraction, aes(x=doy, y=SnowFraction, fill=LocationID, col=LocationID))+
    #    geom_line(data=Locations_SnowFraction_predictions, aes(x=doy, y=SnowFraction_gam_predict, col=LocationID)) + 
    #    xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) + 
    #    ylab(paste0("Snowcover fraction within study area in ", year_ID)) +
    #    theme_tom()
    # 
    #   ggsave(plot=p_Locations_SnowFraction, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_SnowFraction_polygon.pdf"), width=10, height = 8)
  
    #Plot SnowFraction and model predictions in a separate plot per location
      p_Locations_SnowFraction_grid = ggplot()+ 
       geom_point(data=Locations_SnowFraction[Locations_SnowFraction$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
       geom_point(data=Locations_SnowFraction[Locations_SnowFraction$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
       geom_line(data=Locations_SnowFraction_predictions, aes(x=doy, y=SnowFraction_gam_predict), col="#1620de", lwd=1.25)+
       geom_vline(xintercept = 150, colour="grey", lty=2)+
       geom_hline(yintercept = Snowfraction_threshold, colour="grey", lty=2)+
       facet_wrap(~LocationID, ncol=ceiling(length(unique(Locations_SnowFraction$LocationID))^0.5))+
       xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
       ylab(paste0("Snowcover fraction within study area in ", year_ID)) +
       theme_tom()
  
      # ggsave(plot=p_Locations_SnowFraction_grid, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_SnowFraction_grid_polygon.pdf"), width=12, height = 10)
  
   #(E) Calculate at which day of year the SnowFraction is equal to Snowfraction_threshold for each location using predictions from mod_gam

     #Create an empty dataframe     
      df_SnowFraction_Locations <- data.frame(LocationID=character(), 
                                              doy=numeric(), 
                                              stringsAsFactors=FALSE)   
  
     #Loop through all locations in the Locations_SnowFraction_predictions dataframe and extract doy for the required Snowfraction for each location:
      for(i in unique(Locations_SnowFraction_predictions$LocationID)){
    
       #For debugging  
       #i=unique(Locations_SnowFraction_predictions$LocationID)[1] #for debugging
    
       #Select location-specific subset of data:
        df_SnowFraction_gam <- Locations_SnowFraction_predictions[Locations_SnowFraction_predictions$LocationID==i & !is.na(Locations_SnowFraction_predictions$SnowFraction),]
        df_SnowFraction_tmp <- Locations_SnowFraction[Locations_SnowFraction$LocationID==i & !is.na(Locations_SnowFraction$SnowFraction),]
    
       #Detect cutoff points where curve goes from above SnowFraction threshold to below SnowFraction threshold
        df_SnowFraction_gam$cutoff <- ifelse(df_SnowFraction_gam$SnowFraction_gam_predict >= Snowfraction_threshold, 1, 0)
        df_SnowFraction_gam$dif <- c(0, diff(df_SnowFraction_gam$cutoff))
        #the column 'cutoff' indicates whether the gam prediction is above (1) or below (0) the SnowFraction threshold
        #the column 'dif' indicates when there is a change from 1 to 0 (-1) or 0 to 1 (1) in the column cutoff
        #Thus, those rows where the column 'dif' is equal to -1 indicate moments where the SnowFraction value changes from above
        #the threshold to below the threshold. It might be possible that this happens multiple times within a season due to
        #measurement errors or cloud effects. We therefore need to determine which 'cutoff' most likely corresponds to the 
        #actual moment of snowmelt 
    
       #If SnowFraction-threshold was at least crossed once:
        if(any(df_SnowFraction_gam$dif<0)){
      
          #Select all moments (cutoffs) where dif==-1
           cutoffs <- data.frame(index=which(df_SnowFraction_gam$dif<0))
      
          #For the period 30 days after each cutoff point, sum the number of days that have a SnowFraction value larger than Snowfraction_threshold. If a 
          #If a cutoff represents actual snowmelt, then we do not expect any days after this moment with SnowFraction > Snowfraction_threshold. Thus, the
          #closer this sum is to 0, the more likely this cutoff corresponds to the actual moment of snowmelt.
           cutoffs$min <- cutoffs$index -30
           cutoffs$min[cutoffs$min<1] <- 1
           cutoffs$max <- cutoffs$index + 29
           cutoffs$max[cutoffs$max>nrow(df_SnowFraction_gam)] <- nrow(df_SnowFraction_gam)
           cutoffs$sum_cutoff_plus_30 <- apply(cutoffs, 1, function(x){sum(df_SnowFraction_gam$cutoff[x['index']:(x['max'])])})
           cutoff_best <- cutoffs[cutoffs$sum_cutoff_plus_30==min(cutoffs$sum_cutoff_plus_30),'index'][1]
      
          #Approximate day of snowmelt in period from (cutoff_best-1 : cutoff_best) 
           newdata_subset <- df_SnowFraction_gam[max(0, cutoff_best-2) : min(cutoff_best+1, nrow(df_SnowFraction_gam)),]
           doy_approx <- stats::approx(x = newdata_subset$SnowFraction_gam_predict, y = newdata_subset$doy, xout = Snowfraction_threshold)$y[1]
      
           #For debugging:
           # p_tmp <- ggplot()+
           #   geom_point(data=df_SnowFraction_tmp[df_SnowFraction_tmp$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
           #   geom_point(data=df_SnowFraction_tmp[df_SnowFraction_tmp$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
           #   geom_line(data=df_SnowFraction_gam, aes(x=doy, y=SnowFraction_gam_predict), col = "red") +
           #   geom_point(aes(x=doy_approx, y=Snowfraction_threshold), col="blue", size=3)+
           #   geom_hline(yintercept=Snowfraction_threshold, lty=2, col="grey")+
           #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
           #   ylab("SnowFraction-value at pixel") +
           #   theme_classic()+
           #   ggtitle(i)
      
        }
    
       #If threshold was not crossed:
        if(!any(df_SnowFraction_gam$dif<0)){
      
          #No date of snowmelt could be defined
           doy_approx <- NA
      
           #For debugging:
           # p_tmp <- ggplot()+
           #   geom_point(data=df_SnowFraction_tmp[df_SnowFraction_tmp$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
           #   geom_point(data=df_SnowFraction_tmp[df_SnowFraction_tmp$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
           #   geom_line(data=df_SnowFraction_gam, aes(x=doy, y=SnowFraction_gam_predict), col = "red") +
           #   geom_hline(yintercept=Snowfraction_threshold, lty=2, col="grey")+
           #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
           #   ylab("SnowFraction-value at pixel") +
           #   theme_classic()+
           #   ggtitle(i)
      
        }
    
    #Add day of snowmelt and LocationID to dataframe df_SnowFraction_Locations    
     df_SnowFraction_Locations[which(unique(Locations_SnowFraction_predictions$LocationID)==i),'LocationID'] <- as.character(i)
     df_SnowFraction_Locations[which(unique(Locations_SnowFraction_predictions$LocationID)==i),'doy'] <- doy_approx
  
     }
  
   #(F): Plot SnowFraction, model predictions and date of snowmelt in a separate plot per location
     p_Locations_SnowFraction_Snowmelt_grid <- p_Locations_SnowFraction_grid +
       geom_point(data=df_SnowFraction_Locations[!is.na(df_SnowFraction_Locations$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)
  
     ggsave(plot=p_Locations_SnowFraction_Snowmelt_grid, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_SnowFraction_Snowmelt_grid_polygon.pdf"), width=12, height = 10)
  
   #(G): Save the dataframe with snowfraction dates per location   
     
     #Add coordinates to each location  
      df_SnowFraction_Locations <- left_join(df_SnowFraction_Locations, df_locations, by=c("LocationID"))
  
     #Save date when Snowfraction within aoi is equal to 'Snowfraction_threshold' per location as a csv file
      write.csv(df_SnowFraction_Locations, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Snowmelt_SnowFraction_", Snowfraction_threshold, "_polygon.csv"), row.names = FALSE)
  
  }
  
###############################################################################################################################################################
  
 #(MEAN NDSI) - Fit a Generalized Additive Model (GAM) through the NDSI data for each Location

  if("avg_NDSI" %in% method){
       
   #(A) Specify at which NDSI value, the point is considered snow-free
    NDSI_threshold = NDSI_threshold
        
   #(B) Create an empty dataframe
    Locations_NDSI <- data.frame(NDSI=numeric(),
                                 Date=factor(),
                                 doy=numeric(),
                                 LocationID=factor(),
                                 outliers=logical())
        
    Locations_NDSI_predictions <- data.frame(LocationID=character(), 
                                             doy=numeric(), 
                                             NDSI_gam_predict=numeric(), 
                                             stringsAsFactors=FALSE)   
        
   #(C) Loop through all Locations and fit a separate gam with sequential outlier removal to the location specific mean NDSI data
    for(i in unique(Locations_BandValues$LocationID)){
          
      #For debugging  
      #i=unique(Locations_BandValues$LocationID)[1]
          
      #Select Location-specific subset of data:
       df_Location_NDSI <- Locations_BandValues[Locations_BandValues$LocationID==i & !is.na(Locations_BandValues$NDSI), 
                                           c("NDSI", "Date", "doy", "LocationID")]
          
      #Fit a gam through the Location-specific mean_NDSI ~ doy data and employ sequential outlier removal
       df_Location_NDSI <- f_gam_SeqRemOutliers(data=df_Location_NDSI, y="NDSI", x="doy", outlier_removal=outlier_removal, 
                                           outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                           default_k=gam_k_outlier)
          
      #Sort df_Location_NDSI by doy:
       df_Location_NDSI <- df_Location_NDSI[order(df_Location_NDSI$doy),]
          
      #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
       Locations_NDSI <- rbind(Locations_NDSI, df_Location_NDSI)
          
      #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves
          
       #Refit GAM through data
        index <- which(df_Location_NDSI$outliers==FALSE)
        mod_gam <- with(df_Location_NDSI[index,], mgcv::gam(NDSI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))
          
       #Use gam to make predictions on a more detailed (1-day) day of year interval
        aoi_NDSI_predicted <- data.frame(LocationID=i, doy=seq(min(df_Location_NDSI$doy), max(df_Location_NDSI$doy), 1))
        aoi_NDSI_predicted$NDSI_gam_predict <- stats::predict(mod_gam, newdata=aoi_NDSI_predicted, type="response")
        aoi_NDSI_predicted <- aoi_NDSI_predicted[order(aoi_NDSI_predicted$doy),]
          
       #Add predictions to Locations_NDSI_predictions dataframe:  
        Locations_NDSI_predictions <- rbind(Locations_NDSI_predictions, aoi_NDSI_predicted)
          
        }
        
     #Change column LocationID to a factor:
      Locations_NDSI$LocationID <- as.factor(as.character(Locations_NDSI$LocationID))
      Locations_NDSI_predictions$LocationID <- as.factor(as.character(Locations_NDSI_predictions$LocationID))
        
   #(D) Plot the raw NDSI datapoints and gam predictions for each Location:
        
     # #Plot NDSI and model predictions for all locations in a single plot
     #  p_Locations_NDSI = ggplot()+ 
     #    geom_point(data=Locations_NDSI, aes(x=doy, y=NDSI, fill=LocationID, col=LocationID))+
     #    geom_line(data=Locations_NDSI_predictions, aes(x=doy, y=NDSI_gam_predict, col=LocationID)) + 
     #    xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) + 
     #    ylab("NDSI") +
     #    theme_tom()
     #    
     #  ggsave(plot=p_Locations_NDSI, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDSI_polygon.pdf"), width=10, height = 8)
         
     #Plot NDSI and model predictions in a separate plot per location
      p_Locations_NDSI_grid = ggplot()+ 
        geom_point(data=Locations_NDSI[Locations_NDSI$outliers==FALSE,], aes(x=doy, y=NDSI))+
        geom_point(data=Locations_NDSI[Locations_NDSI$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
        geom_line(data=Locations_NDSI_predictions, aes(x=doy, y=NDSI_gam_predict), col="#1620de", lwd=1.25)+
        geom_vline(xintercept = 150, colour="grey", lty=2)+
        geom_hline(yintercept = NDSI_threshold, colour="grey", lty=2)+
        facet_wrap(~LocationID, ncol=ceiling(length(unique(Locations_NDSI$LocationID))^0.5))+
        xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
        ylab("NDSI") +
        theme_tom()
        
      # ggsave(plot=p_Locations_NDSI_grid, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDSI_grid_polygon.pdf"), width=12, height = 10)
        
   #(E) Calculate at which day of year the average NDSI is equal to NDSI_threshold for each location using predictions from mod_gam
        
     #IMPORTANT: Note that these calculations of the average date of snowmelt within the buffer zone around each point of interest are based  
     #on the average NDSI value within this zone. Thus, only a single GAM is fitted through these mean data. A better alternative would be to 
     #fit a separate GAM (with outlier removal) through all pixels in the buffer zone around each point, calculate the date of snowmelt 
     #for each pixel and average those values to get the average date of snowmelt in the bufferzone around each point. The latter approach 
     #can be found in script 5-RGEE_TomVersluijs_PixelDateOfSnowmelt. However, this approach is not feasible for areas larger than 50km2
     #due to computation limits.
        
     #Create an empty dataframe     
      df_Snowmelt_Locations <- data.frame(LocationID=character(), 
                                          doy=numeric(), 
                                          stringsAsFactors=FALSE)   
        
     #Loop through all locations in the Locations_NDSI_predictions dataframe and extract doy of snowmelt for each location:
      for(i in unique(Locations_NDSI_predictions$LocationID)){
          
          #For debugging  
          #i=unique(Locations_NDSI_predictions$LocationID)[1] #for debugging
          
          #Select location-specific subset of data:
          df_NDSI_gam <- Locations_NDSI_predictions[Locations_NDSI_predictions$LocationID==i & !is.na(Locations_NDSI_predictions$NDSI),]
          df_NDSI_tmp <- Locations_NDSI[Locations_NDSI$LocationID==i & !is.na(Locations_NDSI$NDSI),]
          
          #Detect cutoff points where curve goes from above NDSI threshold to below NDSI threshold
          df_NDSI_gam$cutoff <- ifelse(df_NDSI_gam$NDSI_gam_predict >= NDSI_threshold, 1, 0)
          df_NDSI_gam$dif <- c(0, diff(df_NDSI_gam$cutoff))
          #the column 'cutoff' indicates whether the gam prediction is above (1) or below (0) the ndsi threshold
          #the column 'dif' indicates when there is a change from 1 to 0 (-1) or 0 to 1 (1) in the column cutoff
          #Thus, those rows where the column 'dif' is equal to -1 indicate moments where the NDSI value changes from above
          #the threshold to below the threshold. It might be possible that this happens multiple times within a season due to
          #measurement errors or cloud effects. We therefore need to determine which 'cutoff' most likely corresponds to the 
          #actual moment of snowmelt 
          
          #If NDSI-threshold was at least crossed once:
          if(any(df_NDSI_gam$dif<0)){
            
            #Select all moments (cutoffs) where dif==-1
            cutoffs <- data.frame(index=which(df_NDSI_gam$dif<0))
            
            #For the period 30 days after each cutoff point, sum the number of days that have a NDSI value larger than NDSI_threshold. If a 
            #cutoff represents actual snowmelt, then we do not expect any days after this moment with NDSI > NDSI_threshold. Thus, the closer
            #this sum is to 0, the more likely this cutoff corresponds to the actual moment of snowmelt.
            cutoffs$min <- cutoffs$index -30
            cutoffs$min[cutoffs$min<1] <- 1
            cutoffs$max <- cutoffs$index + 29
            cutoffs$max[cutoffs$max>nrow(df_NDSI_gam)] <- nrow(df_NDSI_gam)
            cutoffs$sum_cutoff_plus_30 <- apply(cutoffs, 1, function(x){sum(df_NDSI_gam$cutoff[x['index']:(x['max'])])})
            cutoff_best <- cutoffs[cutoffs$sum_cutoff_plus_30==min(cutoffs$sum_cutoff_plus_30),'index'][1]
            
            #Approximate day of snowmelt in period from (cutoff_best-1 : cutoff_best) 
            newdata_subset <- df_NDSI_gam[max(0, cutoff_best-2) : min(cutoff_best+1, nrow(df_NDSI_gam)),]
            doy_snowmelt <- stats::approx(x = newdata_subset$NDSI_gam_predict, y = newdata_subset$doy, xout = NDSI_threshold)$y[1]
            
            #For debugging:
            # p_tmp <- ggplot()+
            #   geom_point(data=df_NDSI_tmp[df_NDSI_tmp$outliers==FALSE,], aes(x=doy, y=NDSI))+
            #   geom_point(data=df_NDSI_tmp[df_NDSI_tmp$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
            #   geom_line(data=df_NDSI_gam, aes(x=doy, y=NDSI_gam_predict), col = "red") +
            #   geom_point(aes(x=doy_snowmelt, y=NDSI_threshold), col="blue", size=3)+
            #   geom_hline(yintercept=NDSI_threshold, lty=2, col="grey")+
            #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            #   ylab("NDSI-value at pixel") +
            #   theme_classic()+
            #   ggtitle(i)
            
          }
          
          #If threshold was not crossed:
          if(!any(df_NDSI_gam$dif<0)){
            
            #No date of snowmelt could be defined
            doy_snowmelt <- NA
            
            #For debugging:
            # p_tmp <- ggplot()+
            #   geom_point(data=df_NDSI_tmp[df_NDSI_tmp$outliers==FALSE,], aes(x=doy, y=NDSI))+
            #   geom_point(data=df_NDSI_tmp[df_NDSI_tmp$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
            #   geom_line(data=df_NDSI_gam, aes(x=doy, y=NDSI_gam_predict), col = "red") +
            #   geom_hline(yintercept=NDSI_threshold, lty=2, col="grey")+
            #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            #   ylab("NDSI-value at pixel") +
            #   theme_classic()+
            #   ggtitle(i)
            
          }
          
          #Add day of snowmelt and LocationID to dataframe df_Snowmelt_Locations    
          df_Snowmelt_Locations[which(unique(Locations_NDSI_predictions$LocationID)==i),'LocationID'] <- as.character(i)
          df_Snowmelt_Locations[which(unique(Locations_NDSI_predictions$LocationID)==i),'doy'] <- doy_snowmelt
        }
        
    #(F): Plot NDSI, model predictions and date of snowmelt in a separate plot per location
      p_Locations_NDSI_Snowmelt_grid <- p_Locations_NDSI_grid +
          geom_point(data=df_Snowmelt_Locations[!is.na(df_Snowmelt_Locations$doy),], aes(x=doy, y=NDSI_threshold), col="red", size=3)
        
      ggsave(plot=p_Locations_NDSI_Snowmelt_grid, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDSI_Snowmelt_grid_polygon.pdf"), width=12, height = 10)
    
    #(G): Save the dataframe with average snowmelt dates per location     
          
       #Add coordinates to each location  
        df_Snowmelt_Locations <- left_join(df_Snowmelt_Locations, df_locations, by=c("LocationID"))
        
       #Save date of average NDSI within aoi equal to NDSI_threshold per location as a csv file
        write.csv(df_Snowmelt_Locations, paste0("Output/S2/Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Snowmelt_NDSI_", NDSI_threshold, "_polygon.csv"), row.names = FALSE)
  
  }    

##########################################################################################################################################################################

 #The End   

##########################################################################################################################################################################
        