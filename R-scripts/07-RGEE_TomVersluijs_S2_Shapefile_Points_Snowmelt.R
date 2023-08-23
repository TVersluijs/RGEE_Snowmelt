##################################################################################################################################

#In this script the timing of snow melt is calculated based on Sentinel-2 data for all locations specified in an input file. The
#user can specify a bufferzone (radius) to depict the area in which snow melt will be analysed per location. All locations are
#required to be located within a single shapefile and are then analysed simultaneously. First, clouds and permanent water bodies 
#are filtered within the shapefile. Second, if the shapefile overlaps with multiple satellite tiles for a certain day, a composite
#image is created (picking the pixel with least cloudcover). Finally, snow melt is analysed for each location's buffer zone based 
#on one of the following methods (specified by the user by setting the parameter 'method'):

# (I): 'snowfraction': Calculates the fraction of pixels within each buffer zone over time where NDSI > 'NDSI_threshold', fits a
#       GAM through these data and extract the moment when this model passes a user-specified 'Snowfraction_threshold'.
# (II): 'avg_NDSI': calculates the average NDSI value over time within each point's buffer zone, fits a GAM through these data 
#       and calculates when this model passes the specified NDSI threshold defining the moment of snow melt. In addition, time 
#       series of the average NDVI and NDMI are extracted within each point's buffer zone.

#The 'snowfraction' method is preferred, because it intuitively makes sense to look a the fraction of snow-covered pixels
#over time. It is harder to justify the avg_NDSI method, because it is rather unclear what this average NDSI value entails.

#This script is similar to the script '04-RGEE_TomVersluijs_S2_Points_Snowmelt.R'. However, in the latter script all points  
#are analysed consecutively using a loop, which makes that script significantly slower to run. It does however not rely on a 
#shapefile and thus works for points spaced much further apart (i.e. tracking data of birds migratory birds). The current
#script only works for small areas of c.a. 50-100 km2 (larger areas might result in computation errors).

#Copyright Tom Versluijs 2023-07-31. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R". 
#Note that a GoogleDrive is required.

##################################################################################################################################

#I: Setup workspace

##################################################################################################################################

      #(0): Clear workspace and set python environment
       rm(list=ls())
       utils::install.packages("here")
       library(here)
       if(file.exists(paste0(here::here(), "/Input/rgee_environment_dir.rds"))){
         rgee_environment_dir <- readRDS(paste0(here::here(), "/Input/rgee_environment_dir.rds"))
         reticulate::use_python(rgee_environment_dir, required=T)
         reticulate::py_config()
         }
        
      #(1): Load packages
       #renv::restore() #revert to last version of R-packages used to successfully run this script (optional).
       utils::install.packages("pacman")
       library(pacman)
       p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr)

      #(2): Define ggplot2 plotting theme
        theme_tom <- function(){
           theme_classic() %+replace%
             theme(axis.title = element_text(size=18),
                   axis.text = element_text(size=16),
                   legend.position = "none",
                   strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2)),
                   complete = TRUE)}    
        
      #(3): Load auxiliary functions
        source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "Sentinel2_AuxiliaryFunctions")
        sapply(source_files, source, chdir = TRUE) ; rm(source_files)
        
      #(4): Initialize earth engine
        rgee::ee_Initialize(user = "tom.versluijs@gmail.com", drive = TRUE)
        
##################################################################################################################################
       
#II: Specify parameters of interest
       
##################################################################################################################################       
       
#(4): Manually specify parameters of interest

  #(a): Sentinel-2 satellite

    #Sentinel-2 dataset ("S2_HARMONIZED" or "S2_SR_HARMONIZED")
    s2_dataset <- "S2_SR_HARMONIZED"
    #Use S2_HARMONIZED (level 1c product) for data from 2016 - 2018
    #Use S2_SR_HARMONIZED (level 2a product) for data starting from 2019

    #Resolution of sampling in meters
    resolution=20 #default maximum resolution for Sentinel-2 = 10m

  #(b) Area of interest

    #Specify name of study area (used as prefix in output files)
    area_name="ZAC"

    #Name of Shapefile located in Input folder (specific area of interest)
	#Follow the guide "Manual_CreateShapefilePolygons.docx" when creating this shapefile.
    shapefile <- "ZAC_Outline_EPSG4326.shp"

    #Coordinate reference system used for calculations
    #EPSG:4326 is recommended for areas spanning multiple UTM zones, but increased computation time (i.e. spherical coordinate system).
    #EPSG:326XX is results in reduced computation time for areas located within a single UTM zone (i.e. planar coordinate system).
    crs <- "EPSG:4326"

  #(c) Point locations

    #Name of file with point locations located in Input folder
    input_locations <- "TestLocations_Zackenberg.csv"
    #Make sure it has the columns "LON_x", "LAT_y" and "DateTime"
    #DateTime should be in format 'dd/mm/yyyy hh:mm:ss'

    #Buffer radius around each point location (in meters)
    Buffer_radius_m=250

  #(c) Dates

    #Specify the year of interest:
    year_ID <- "2022"

    #Date range of all images considered for analysis
    start_date <- paste0(year_ID, "-03-15") #choose date (well) before the first snowmelt occurs within the study site (default -03-15)
    end_date <- paste0(year_ID, "-09-15") #choose date (well) after last date of tracking (default -09-15)

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
    #Note that cloud masking will automatically be set to TRUE when mask_water=TRUE
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
    resolution_cldmsk=150 #(default = 150)

    #Maximum fraction of cloud cover allowed in each image
    max_cloud_fraction=1.0

  #(f): Water masking

    #Should permanent waterbodies be masked from the analysis (default=TRUE).
    mask_water=TRUE
    mask_water_type="both" #either "water_mask_ESA", "water_mask_Manual", or "both"
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
    
  #(5): Automatically define some additional parameters 
    
    #Create a unique data_ID
    data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_S2")           
    
    #Set mask_clouds to TRUE if a composite image needs to be generated         
    if(create_composite==TRUE){mask_clouds=TRUE}      
    
    #Set mask_clouds to TRUE if water masking is TRUE and mask_water_type==water_mask_Manual
    if(mask_water==TRUE & (mask_water_type=="water_mask_Manual" | mask_water_type=="both")){mask_clouds=TRUE}      
    
    #Create output folder
    if(dir.exists(paste0(here(), "/Output/S2/Shapefile_Points_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/S2/Shapefile_Points_Snowmelt"), recursive = TRUE)}
    
    
##################################################################################################################################

#III: Read and display the unfiltered data

##################################################################################################################################

    #(6) Read aoi_Shapefile and convert to a feature collection
      root_fldr <- here()
      aoi_Shapefile <- st_read(paste0(root_fldr, "/Input/Shapefiles/", shapefile), quiet=T)
      aoi_Shapefile <- st_transform(aoi_Shapefile, crs="EPSG:4326")
      aoi_Shapefile <- sf_as_ee(aoi_Shapefile)
      aoi_Shapefile <- ee$FeatureCollection(aoi_Shapefile)
      #This feature collection can be used to clip the satellite imagery to the study area
      
       #Calculate bounding box for aoi_Shapefile
        aoi <- aoi_Shapefile$geometry()$bounds()
        
       #Calculate central point
        coordinates_point <- aoi_Shapefile$geometry()$centroid()
        
       #Plot shapefile and bounding box (for debugging)
        Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 6)
        Map$addLayer(aoi_Shapefile)+
          Map$addLayer(ee$FeatureCollection(aoi)$style(color='red', fillColor='00000000'))+
          Map$addLayer(coordinates_point, list(color="blue"))
        
       #Calculate the size of the study area in km2:
        img <- ee$Image$pixelArea()$divide(1000000)
        area2 <- img$reduceRegion(
          reducer= ee$Reducer$sum(),
          geometry= aoi_Shapefile,
          crs=crs,
          scale= resolution,
          maxPixels= 1E13)
        paste0('Size of study area of calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')
      
    #(7) Plot an RGB, NDSI and NDVI image for the first extracted satellite image

      #Extract Sentinel-2 Surface Reflectance data for a given area and date range:
      s2_col <- ee$ImageCollection(paste0('COPERNICUS/', s2_dataset))
      s2_col <- s2_col$
        filterBounds(aoi)$ #For areas that do not span multiple tiles, change aoi to coordinates_point to greatly speed-up code.
        filterDate(start_date, end_date)
  
      #Add Sentinel-2 cloud probability as a new property to s2_col
      if(mask_clouds==TRUE){
        
        #Load cloud information from S2_CLOUD_PROBABILITY
        s2_cloudless_col <-ee$ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
        s2_cloudless_col <- s2_cloudless_col$
          filterBounds(aoi)$ #For areas that do not span multiple tiles, change aoi to coordinates_point to greatly speed-up code.
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
      
      #Clip all images to the area depicted by aoi_Shapefile
      s2_col <- s2_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})
      
    #(8): Plot an RGB, NDSI and NDVI image for a single extracted satellite image (for debugging)
      
      #Select a single image for initial plot:
      image <- s2_col$filterDate(paste0(year_ID, "-07-15"), end_date)$first()
      
      #Normalized difference snow index:
      ndsi_s2 <-  image$expression('(B3 - B11) / (B3 + B11)',
                                   list('B11'=image$select('B11'),
                                        'B3'=image$select('B3')))
      
      #Normalized difference vegetation index:
      ndvi_s2 <- image$expression('(B8 - B4) / (B8 + B4)',
                                  list('B8'=image$select('B8'),
                                       'B4'=image$select('B4')))
  
      #Add Normalized difference water index to this image:
      ndwi_s2 <- image$expression('(B3 - B8) / (B3 + B8)',
                                  list('B3'=image$select('B3'),
                                       'B8'=image$select('B8')))
      
      #Plot all layers
      Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
      Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
      Map$addLayer(ndsi_s2,list(min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
      Map$addLayer(ndvi_s2,list(min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
      Map$addLayer(ndwi_s2,list(min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')
      
    #(9): Define a featurecollection of points for which Sentinel-2 data needs to be extracted
      ZAC_Paths <- read.csv(paste0("Input/", input_locations), header=T)[,c("LON_x", "LAT_y", "DateTime")]
      colnames(ZAC_Paths) <- c("LON_x", "LAT_y", "DateTime")
      
      #Add a unique LocationID to every unique lat/lon combination
      Locations_sf <- unique(ZAC_Paths[,c("LON_x", "LAT_y")])
      Locations_sf$LocationID <- paste0("Location_", 1:nrow(Locations_sf))
      ZAC_Paths <- left_join(ZAC_Paths, Locations_sf, by=c("LON_x", "LAT_y"))
      
      Locations  <-  st_as_sf(Locations_sf, coords = c("LON_x", "LAT_y"), crs = crs)
      Locations <- st_transform(Locations, crs=crs)
      Locations <- sf_as_ee(Locations[,c("LocationID", "geometry")])
      
      #Add buffer around all point locations
      if(Buffer_radius_m>0){
        bufferBy <- function(Buffer_radius_m) {
        return(function(feature) {
          return(feature$buffer(Buffer_radius_m))   
          })
        }
        Locations <- Locations$map(bufferBy(Buffer_radius_m))
        }
      
      #Locations is now a feature collection. A feature in our analysis comprises the feature type (point in case no buffer
      #was added, or a circular polygon with radius Buffer_radius_m with a center the original point location), Location_ID  
      #and the feature's geometry as lat and lon coordinates. A feature collection corresponds to a collection of such features 
      #(i.e. a collection of different Locations, with their respective feature types, LocationIDs, DateTimes and coordinates).
      
      #Plot feature collection on a map (i.e. locations as points)
      image <- s2_col$filterDate(paste0(year_ID, "-06-15"), end_date)$first()
      Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
      Map$addLayer(aoi_Shapefile)+
      Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
      Map$addLayer(Locations, list(color="red"), paste0("Locations_", year_ID))    
      
    #(10) Add NDSI, NDVI, NDMI and NDWI to the clipped image collection
      
     #Map over image collection
      s2_col <- s2_col$
        map(getNDSI)$
        map(getNDVI)$
        map(getNDMI)$
        map(getNDWI)
    
    #(11) Create a timeseries gif of RGB images for the shapefile area
      
      #Check number of images in collection (for debugging)
      #s2_col$size()$getInfo()
      
      #Create a timelapse video of RGB band (for debugging)
      # videoArgs <- list(dimensions=350, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
      # tryCatch({browseURL(s2_col$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
        

##################################################################################################################################
        
#IV: Filter and mask clouds within the image collection
        
##################################################################################################################################
        
   #There are many images that are almost fully covered in dense clouds. We thus first filter and mask clouds in the image collection.
   #We delete all images with a cloud cover fraction >= max_cloud_fraction. In the remaining images we mask all pixels that are 
   #covered by clouds.
    if(mask_clouds==TRUE){   
        
      #print message   
      print("Cloud masking = TRUE")
      
      #(11) Add cloud information within the aoi_Shapefile to the image collection: 

        #Add the cloud fraction within the shapefile area to each separate image by mapping the cloud functions over the image collection
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
        
        #Check if CloudFraction property has been added to each image
        #s2_col$first()$propertyNames()$getInfo()
        
        # #Visual check of cloud mask (for debugging)
        # image <- s2_col$filterDate(paste0(year_ID, "-05-31"), paste0(year_ID, "-09-14"))$first()$
        #   clipToCollection(aoi_Shapefile)$
        #   select("clouds_unadjusted", "clouds", "clouds_inv", "clouds_prob", "B4", "B3", "B2")      
        # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
        # Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
        #   #Map$addLayer(image, list(bands='clouds_prob', min=0, max=100, opacity=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'clouds_prob')+
        #   Map$addLayer(image, list(bands='clouds_unadjusted', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'clouds unadjusted')+
        #   Map$addLayer(image, list(bands='clouds', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'clouds buffered')
        # #Map$addLayer(image, list(bands='clouds_inv', min=0, max=1, opacity=1), 'clouds_inv')
        
        #We have now calculated the fraction of cloud covered pixels within our shapefile area for each image in the image collection.
        #However, Sentinel-2 satellite data is recorded as different tiles, these tiles are generally much bigger than a single study site.
        #Each Sentinel-2 tile also contains a native measure of cloud cover. This is the property 'CLOUDY_PIXEL_PERCENTAGE'. This property
        #thus refers to the percentage of cloud cover within the whole tile. We have now calculated a more appropriate measure for our
        #shapefile area only called 'CloudFraction'. Below we compare both measures, as we expect they should be correlated. Differences
        #between both measures can for example arise when the shapefile area was covered by clouds while most of the larger region (tile)
        #was not or vice versa.
        
        # #Extract cloud fraction of all images in image collection for the aoi_Shapefile
        # cloud_aoi <- unlist(s2_col$aggregate_array('CloudFraction')$getInfo())
        # 
        # #Extract cloud fraction of all images in image collection for the complete tile to which our aoi_Shapefile belongs
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
        #   #Thus, this indicates that our local cloud measure does not necessarily correspond to that as measured on the whole tile as
        #   #explained above.
      
       #(12): Exclude all images from the image collection that have a CloudFraction value > max_cloud_fraction, and mask all cloud pixels within the remaining images
    
        #In the analysis for the whole study area we exclude images with more than 20% snowcover within our shapefile 
        #However, now we're interested in individual pixels. So even when 90% of the shapefile area is snowcovered,
        #the pixels we're interested in might be located in the 10% cloud free area. We sthus only filter images with
        #a 100% cloudcover (cloudFraction=1.0). We apply a cloud mask to mask out cloudy pixels within each image.

        #Mask cloudy pixels within each image: 
         s2_clouds_filtered <- s2_col$
            #Only remove images that are fully cloud-covered
            filter(ee$Filter$lt('CloudFraction', max_cloud_fraction))$
            #Apply cloudmask for individual pixels
            map(Add_CloudMask)
         
         # #Create timelapse video of the cloud filtered/masked RGB images (for debugging)
         # videoArgs <- list(dimensions=400, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
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
             
##################################################################################################################################

#V: Mask permanent waterbodies (ponds, lakes, rivers, sea) within the image collection

##################################################################################################################################

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
         # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
         # Map$addLayer(s2_col$filterDate(paste0(year_ID, "-08-04"), paste0(year_ID, "-09-18"))$first(), list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
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

##################################################################################################################################

#VI: Create a composite image for each day by mosaicking all images from the that day

##################################################################################################################################

  #There might be days for which multiple satellite photos are available. We thus need to make a decision on what to do with overlapping
  #values per pixel. We use a qualityMosaic that can select the overlapping pixel with the least cloudcover by tapping into the band
  #"cloud_inv". If two pixels have equal cloudcover, then the pixelvalue from this first image is chosen.

  #Note that cloud/water masking is taken into account when making the composite image. The reducer function "qualitymosaic" will  for every pixel only use
  #the unmasked pixels to calculate the output value (i.e. if for a certain pixel there are 4 overlapping images available for the same date and one
  #of them is covered by clouds (and thus masked), then only the remaining cloud-free (and unmasked) pixels are used to calculate the band values
  #at that pixel). See, https://gis.stackexchange.com/questions/337666/applying-google-earth-engine-cloud-mask-on-each-individual-image-before-mosaicin,
  #or: https://gis.stackexchange.com/questions/371356/combining-single-cloud-masked-images-using-google-earth-engine, or
  #https://open-mrv.readthedocs.io/en/latest/image_composite_Web.html.

   #F.1: Make a composite image for a single date (doy) and calculate fraction of snowcover - USED FOR DEBUGGING ONLY
    {
    # #Pick a date and check whether there are multiple tiles spanning the aoi
    # doy=165
    # s2_img <- s2_clouds_filtered$filterMetadata('year', 'equals', as.numeric(year_ID))$filterMetadata('doy', 'equals', doy)
    # s2_img <- s2_img$select(c('B4', 'B3', 'B2', 'NDSI', 'clouds_inv'))
    # 
    # # #Number of images in s2_clouds_filtered
    # # s2_img$size()$getInfo()
    # # s2_img$first()$bandNames()$getInfo()
    # # s2_img$first()$propertyNames()$getInfo()
    # 
    # # #Inspect whether there is any overlap between the tiles for this day:
    # # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
    # # Map$addLayers(s2_img, list(bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7)), 'S2_Original')+
    # # #Map$addLayers(s2_img, list(bands=c("seconds"), min=15053255, max=15053280), opacity=0.75, 'image_ID')
    # # Map$addLayer(aoi)
    # 
    # #Create a composite image for this day
    # s2_composite <- ee$Image(s2_img$
    #                            #Create a quality mosaic based on cloud_inv
    #                            qualityMosaic('clouds_inv')$
    #                            #Copy timestamp properties of first image to mosaic image
    #                            copyProperties(s2_img$first(), c('system:id', 'system:time_start')))#$reproject(crs=crs, crsTransform=NULL, scale=resolution)
    # s2_composite <- s2_composite$clip(aoi)
    # #s2_composite$select('NDSI')$projection()$getInfo()
    # #Note that there is no need to reproject as this is done automatically when plotting the data, or when extracting other output.
    # 
    # #s2_composite$bandNames()$getInfo()
    # #s2_composite$propertyNames()$getInfo()
    # #s2_composite$select('NDSI')$projection()$getInfo()
    # 
    # # #Visual checks
    # #  #Inspect the RGB band of the Mosaicked image
    # #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
    # #  Map$addLayer(s2_composite, list(bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7)), 'S2_Median')
    # #
    # #  #Inspect the "clouds" band of the Mosaicked image
    # #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
    # #  Map$addLayer(s2_composite, list(bands=c("clouds_inv"), palette=c("red", "green"), min=0, max=1), opacity=1, 'Clouds_median')
    # #
    # #  #Inspect to which original image each pixel in the Mosaicked image belongs (using the "seconds" band)
    # #  Map$addLayer(s2_composite, list(bands=c("seconds"), min=15053255, max=15053280), opacity=1, 'image_ID')
    # 
    # #Calculate fraction of snowcover for this composite image:
    # 
    #   #Select the NDSI band
    #   ndsi <- s2_composite$select('NDSI')
    # 
    #   #Create a binary layer using logical operations.
    #   snow <- ndsi$gt(NDSI_threshold)$rename('SNOW')
    # 
    #   #Return the binary snow parameter
    #   s2_composite <- s2_composite$addBands(snow)
    # 
    #   # #Mask and display the binary layer (for debugging).
    #   #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
    #   #  Map$addLayer(s2_composite,list(bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
    #   #  Map$addLayer(s2_composite$select('SNOW'))
    # 
    #   #Calculate the fraction of pixels within aoi that are marked as snow
    #   SnowFraction = snow$reduceRegion(
    #     reducer = ee$Reducer$mean(),
    #     geometry = aoi_Shapefile,
    #     scale = resolution
    #     )$get('SNOW')
    # 
    #   # #Print snowfraction
    #   # SnowFraction$getInfo()
    # 
    #   #Add SnowFraction as an image property to each image
    #   s2_composite <- s2_composite$set('SnowFraction', SnowFraction)
    # 
    #   #Set the SnowFraction value to a no data value of -9999 if equal to NULL
    #   SnowFraction <- ee$List(list(s2_composite$get('SnowFraction'), -9999))$reduce(ee$Reducer$firstNonNull())
    # 
    #   #Return the updated SnowFraction property
    #   s2_composite <- s2_composite$set("SnowFraction", SnowFraction)
    # 
    #   #s2_composite$propertyNames()$getInfo()
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
          #Select required bands only
          select(c('B4', 'B3', 'B2', 'NDSI', 'NDVI', 'NDMI', 'NDWI', 'clouds', 'clouds_inv'))
          #select(c('NDSI', 'clouds_inv'))

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
      
##################################################################################################################################

#VI: Extract the fraction of snowcovered pixels within each point's buffer zone for all images in the image collection        

##################################################################################################################################

 #Method 'SnowFraction'
  if("snowfraction" %in% method){

   #Add a binary snow cover band to each image and calculate the fraction of snow covered pixels within each point's buffer zone

    #Map Snow computation functions over the mosaicked image collection
    s2_col_composite <- s2_col_composite$
      #Determine which pixels are snow-covered (NDSI > NDSI threshold)
      map(computeSnow)
   
   #Create an iteration function that we will use to iterate through all images of the image collection. For each image, the 
   #average value of the "SNOW" band is extracted at each feature (i.e. location) within the feature collection Locations.
   #The resulting fractino of snowcover and the datetime of the image are added as a property to each feature (i.e. location).
   #This results in an updated feature collection specific for the current image. This feature collection is then appended to a list
   #of feature collections from previous image iterations. After iterating through n-images, the result is thus a feature
   #collection list of length n * the number of features within the feature collection.
   
   #Thus, at each iteration we extract the average value of the "SNOW" band (i.e. the fraction of snowcover as it comprises of
   #zeros and ones only) for all locations for the current image, store this as a  feature collection and add this to an expanding 
   #list of feature collections from previous iterations.
   
   #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
   FC_initial_snowfraction <- ee$FeatureCollection(ee$List(list()))
    
   #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
   #(in this case the current iteration image) and the second element takes the output value from the iteration that preceded it. The latter
   #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with 
   #has been defined.
   Extract_SnowFractionAtPoints = function(img, FC_initial_snowfraction) {
     
           #(A): Extract Snowfraction within the current image at each point in the Locations feature collection. We add this
           #fraction as a property to each feature (i.e. location). The resulting output is a feature collection for the current
           #img where the snowfraction is added as property for each location.
           FC_image <- img$reduceRegions(collection = Locations, 
                                         reducer = ee$Reducer$mean()$setOutputs(list("SnowFraction")),
                                         scale = resolution, #implicitly specified in crs below
                                         crs=crs, #or, crs=S2Projection
                                         crsTransform=NULL)
           
           #(B): reduceRegions does not return any SnowFraction value when the point is masked by clouds or does not fall within
           #the boundaries of the image. Therefore we manually have to add an empty SnowFraction property to those features. We
           #therefore set the SnowFraction value to a no data value of -9999 for all features where the SnowFraction value is NULL.
           FC_image <- FC_image$map(function(feature){
             snowfraction <- ee$List(list(feature$get('SnowFraction'), -9999))$reduce(ee$Reducer$firstNonNull())
             return(feature$set("SnowFraction", snowfraction))
             })

           #(C): Add date and day of year of current image as a property to each feature within the feature collection FC_image
           date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
           doy <- ee$Date(img$get('system:time_start'))$getRelative('day', 'year')
           FC_image <- FC_image$map(function(feature){return(feature$set("Date", date)$set("doy", doy))})
           
           #(D): Merge the feature collection of the current image (FC_image) onto the feature collection FC_initial_snowfraction
           return (ee$FeatureCollection(FC_initial_snowfraction)$merge(FC_image))
           #FC_initial_snowfraction is thus updated at each iteration. This corresponds to base R code: 
           #FC_initial_snowfraction <- merge(FC_initial_snowfraction, FC_image)
         } 
  
   #Iterate over the ImageCollection and select appropriate band values
   FC_merged_SnowFraction <- ee$FeatureCollection(s2_col_composite$select("SNOW")$iterate(Extract_SnowFractionAtPoints, FC_initial_snowfraction))

   #Export the feature collection as a .csv table
   #We export the data instead of using aggregate_array() as the latter might fail due to computation timeouts.

     #Setup task
      task_vector <- ee_table_to_drive(
        collection= FC_merged_SnowFraction,
        description = paste0(data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_SnowFraction"),
        folder="RGEE_tmp",
        fileFormat="CSV",
        selectors=c('SnowFraction', 'Date', 'LocationID')
        )

     #Run and monitor task
      print(paste0("Uploading the snowfraction data to Google Drive"))
      task_vector$start()
      ee_monitoring(task_vector, max_attempts=1000000) #250s at 100m

     #Import results
      exported_stats <- ee_drive_to_local(task=task_vector, dsn=paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_SnowFraction"))
      Locations_Data_SnowFraction <- read.csv(exported_stats)
      unlink(exported_stats)

     #Add day of year
      Locations_Data_SnowFraction$Date <- as.POSIXct(Locations_Data_SnowFraction$Date, format="%Y-%m-%d %H:%M:%S")
      Locations_Data_SnowFraction$Date_doy <- as.numeric(strftime(Locations_Data_SnowFraction$Date, format = "%j"))

     #Remove NAs in the SnowFraction variable
      Locations_Data_SnowFraction$SnowFraction[Locations_Data_SnowFraction$SnowFraction < -9000] <- NA #replace -9999 by NA
      #Locations_Data_SnowFraction <- Locations_Data_SnowFraction[!is.na(Locations_Data_SnowFraction$SnowFraction), ]

     #Sort dataframe by LocationID and Date_doy
      index <- with(Locations_Data_SnowFraction, order(LocationID, Date_doy))
      Locations_Data_SnowFraction <- Locations_Data_SnowFraction[index,]
      Locations_Data_SnowFraction <- Locations_Data_SnowFraction[,c("LocationID", "Date", "Date_doy", "SnowFraction")]

     #Add coordinates to each location
      Locations_Data_SnowFraction <- left_join(Locations_Data_SnowFraction, Locations_sf[,c("LON_x", "LAT_y", "LocationID")], by=c("LocationID"))

     #Save dataframe
      write.csv(Locations_Data_SnowFraction, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_SnowFraction.csv"), row.names = FALSE)

    }
   
##################################################################################################################################
      
#VII: Extract Sentinel-2 band values for points of interest (with buffer) for all images in the image collection (output a feature collection)  
      
##################################################################################################################################       
     
 #Method 'avg_NDSI'       
  if("avg_NDSI" %in% method){       
         
      #Extract Sentinel-2 band values for points of interest (with buffer) for all images in the image collection (output a feature collection)  
         
       #Create an iteration function that we will use to iterate through all images of the image collection. For each image, the value
       #of certain bands is extracted at each feature (i.e. location) within the feature collection Locations. The resulting
       #band values (and the datetime of the image) are added as a property to each feature (i.e. location). This results
       #in an updated feature collection specific for the current image. This feature collection is then appended to a list
       #of feature collections from previous image iterations. After iterating through n-images, the result is thus a feature
       #collection list of length n * the number of features within the feature collection.
           
       #Thus, at each iteration we extract band values of interest of all locations for the current image, store this as a  
       #feature collection and add this to an expanding list of feature collections from previous iterations.
           
       #Store default Sentinel-2 image projection 
        S2Projection <- s2_col_composite$first()$select("NDSI")$projection()
        #S2Projection$getInfo() 
       
       #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
        FC_initial_meanbandvalues <- ee$FeatureCollection(ee$List(list()))
           
       #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
       #(in this case the current iteration image) and the second element takes the output value from the iteration that preceded it. The latter
       #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with 
       #has been defined.
        Extract_BandValuesAtPoins = Extract_BandValuesAtPoins #sourced
           
       #Iterate over the ImageCollection and select appropriate band values
        FC_merged_meanbandvalues <- ee$FeatureCollection(s2_col_composite$select("NDSI", "NDVI", "NDMI")$iterate(Extract_BandValuesAtPoins, FC_initial_meanbandvalues))
         
      #Transform feature collection with Band values to a dataframe     
           
        #Export the feature collection as a .csv table
        #We export the data instead of using aggregate_array() as the latter might fail due to computation timeouts.
        
        #Setup task
        task_vector_2 <- ee_table_to_drive(
          collection= FC_merged_meanbandvalues,
          description = paste0(data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_MeanBandvalues"),
          folder="RGEE_tmp",
          fileFormat="CSV",
          selectors=c('NDSI', 'NDVI', 'NDMI', 'Date', 'LocationID')
        )
        
        #Run and monitor task
        print(paste0("Uploading mean BandValue data to Google Drive"))
        task_vector_2$start()
        ee_monitoring(task_vector_2, max_attempts=1000000) #250s at 100m
        
        #Import results
        exported_stats <- ee_drive_to_local(task=task_vector_2, dsn=paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_MeanBandvalues"))
        Locations_Data_BandValues <- read.csv(exported_stats)
        unlink(exported_stats)
        
        #Add day of year
        Locations_Data_BandValues$Date <- as.POSIXct(Locations_Data_BandValues$Date, format="%Y-%m-%d %H:%M:%S")
        Locations_Data_BandValues$Date_doy <- as.numeric(strftime(Locations_Data_BandValues$Date, format = "%j"))
        
        #Remove NAs in the NDSI, NDVI and NDMI variables
        Locations_Data_BandValues$NDSI[Locations_Data_BandValues$NDSI < -9000] <- NA
        Locations_Data_BandValues$NDVI[Locations_Data_BandValues$NDVI < -9000] <- NA
        Locations_Data_BandValues$NDMI[Locations_Data_BandValues$NDMI < -9000] <- NA
        
        #Sort dataframe by LocationID and Date_doy
        index <- with(Locations_Data_BandValues, order(LocationID, Date_doy))
        Locations_Data_BandValues <- Locations_Data_BandValues[index,]
        Locations_Data_BandValues <- Locations_Data_BandValues[,c("LocationID", "Date", "Date_doy", "NDSI", "NDVI", "NDMI")]
        
        #Add coordinates to each location
        Locations_Data_BandValues <- left_join(Locations_Data_BandValues, Locations_sf[,c("LON_x", "LAT_y", "LocationID")], by=c("LocationID"))
        
        #Save dataframe
        write.csv(Locations_Data_BandValues, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Data_MeanBandvalues.csv"), row.names = FALSE)
        
    }
      
##################################################################################################################################
           
#VIII: Fit a Generalized Additive Model (GAM) through the Sentinel-2 band value data for each location
           
##################################################################################################################################       
 
 #Method "SnowFraction"        
  if("snowfraction" %in% method){
   
    #(21): SNOW - Fit a Generalized Additive Model (GAM) through the SnowFraction data for each Location

        #We fit a Generalized Additive Model (GAM) through the data for each point location. We do this using a sequential
        #outlier-filtering process. We first fit a GAM through the data, calculating model predictions and residuals. We then
        #exclude all rows from the dataframe where the residual >= (0.4 * the range of the data). We then re-fit a GAM to
        #this reduced dataset, make predictions and calculate residuals. We then exclude all rows from the reduced dataframe
        #where the residual >= (0.2 * the range of the data). This gives us a final dataframe through which we fit a third
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

        #(a) Specify for which SnowFraction value the date of occurence will be calculated
            Snowfraction_threshold = Snowfraction_threshold

        #(b) Create an empty dataframe
            Locations_SnowFraction <- data.frame(SnowFraction=numeric(),
                                                 Date=factor(),
                                                 Date_doy=numeric(),
                                                 LocationID=factor(),
                                                 outliers=logical())

            Locations_SnowFraction_predictions <- data.frame(LocationID=character(),
                                                             Date_doy=numeric(),
                                                             SnowFraction_gam_predict=numeric(),
                                                             stringsAsFactors=FALSE)

        #(c) Loop through all Locations and fit a separate gam with sequential outlier removal to the location specific SnowFraction data
            for(i in unique(Locations_Data_SnowFraction$LocationID)){

              #For debugging
               #i=unique(Locations_Data_SnowFraction$LocationID)[1]

              #Select Location-specific subset of data:
               df_Location <- Locations_Data_SnowFraction[Locations_Data_SnowFraction$LocationID==i & !is.na(Locations_Data_SnowFraction$SnowFraction),
                                                   c("SnowFraction", "Date", "Date_doy", "LocationID")]

              #Fit a gam through the Location-specific SnowFraction ~ Date_doy data and employ sequential outlier removal
               df_Location <- f_gam_SeqRemOutliers(data=df_Location, y="SnowFraction", x="Date_doy", outlier_removal=outlier_removal,
                                                   outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                   default_k=gam_k_outlier)

              #Sort df_Location by Date_doy:
               df_Location <- df_Location[order(df_Location$Date_doy),]

              #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
               Locations_SnowFraction <- rbind(Locations_SnowFraction, df_Location)

              #Create more detailed predictions (not only at the Date_doy present in the datframe) at a 1-day interval to plot more smooth curves

               #Refit GAM through data
                index <- which(df_Location$outliers==FALSE)
                mod_gam <- with(df_Location[index,], mgcv::gam(SnowFraction ~ s(Date_doy, k=min(gam_k, length(index)-1)), method="REML"))

               #Use gam to make predictions on a more detailed (1-day) day of year interval
                aoi_SnowFraction_predicted <- data.frame(LocationID=i, Date_doy=seq(min(df_Location$Date_doy), max(df_Location$Date_doy), 1))
                aoi_SnowFraction_predicted$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=aoi_SnowFraction_predicted, type="response")
                aoi_SnowFraction_predicted <- aoi_SnowFraction_predicted[order(aoi_SnowFraction_predicted$Date_doy),]

               #Add predictions to Locations_SnowFraction_predictions dataframe:
                Locations_SnowFraction_predictions <- rbind(Locations_SnowFraction_predictions, aoi_SnowFraction_predicted)

              }

            #Change column LocationID to a factor:
             Locations_SnowFraction$LocationID <- as.factor(as.character(Locations_SnowFraction$LocationID))
             Locations_SnowFraction_predictions$LocationID <- as.factor(as.character(Locations_SnowFraction_predictions$LocationID))

        #(d) Plot the raw SnowFraction datapoints and gam predictions for each Location:

            #Plot SnowFraction and model predictions for all locations in a single plot
             p_Locations_SnowFraction = ggplot()+
               geom_point(data=Locations_SnowFraction, aes(x=Date_doy, y=SnowFraction, fill=LocationID, col=LocationID))+
               geom_line(data=Locations_SnowFraction_predictions, aes(x=Date_doy, y=SnowFraction_gam_predict, col=LocationID)) +
               xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
               ylab("SnowFraction") +
               theme_tom()

               #ggsave(plot=p_Locations_SnowFraction, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_SnowFraction.pdf"), width=10, height = 8)

            #Plot SnowFraction and model predictions in a separate plot per location
             p_Locations_SnowFraction_grid = ggplot()+
               geom_point(data=Locations_SnowFraction[Locations_SnowFraction$outliers==FALSE,], aes(x=Date_doy, y=SnowFraction))+
               geom_point(data=Locations_SnowFraction[Locations_SnowFraction$outliers==TRUE,], aes(x=Date_doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
               geom_line(data=Locations_SnowFraction_predictions, aes(x=Date_doy, y=SnowFraction_gam_predict), col="#1620de", linewidth=1.25)+
               geom_vline(xintercept = 150, colour="grey", lty=2)+
               geom_hline(yintercept = Snowfraction_threshold, colour="grey", lty=2)+
               facet_wrap(~LocationID, ncol=ceiling(length(unique(Locations_SnowFraction$LocationID))^0.5))+
               xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
               ylab("SnowFraction") +
               theme_tom()

               ggsave(plot=p_Locations_SnowFraction_grid, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_SnowFraction_grid.pdf"), width=12, height = 10)

        #(e) Calculate at which day of year the SnowFraction is equal to Snowfraction_threshold for each location using predictions from mod_gam

            #Create an empty dataframe
             df_Snowmelt_SnowFraction_Locations <- data.frame(LocationID=character(),
                                                              Date_doy=numeric(),
                                                              stringsAsFactors=FALSE)

            #Loop through all locations in the Locations_SnowFraction_predictions dataframe and extract doy of snowmelt for each location:
             for(i in unique(Locations_SnowFraction_predictions$LocationID)){

               #For debugging
                #i=unique(Locations_SnowFraction_predictions$LocationID)[1] #for debugging

               #Select location-specific subset of data:
                df_Location_gam <- Locations_SnowFraction_predictions[Locations_SnowFraction_predictions$LocationID==i & !is.na(Locations_SnowFraction_predictions$SnowFraction),]
                df_Location_tmp <- Locations_SnowFraction[Locations_SnowFraction$LocationID==i & !is.na(Locations_SnowFraction$SnowFraction),]

                #Detect cutoff points where curve goes from above SnowFraction threshold to below SnowFraction threshold
                df_Location_gam$cutoff <- ifelse(df_Location_gam$SnowFraction_gam_predict >= Snowfraction_threshold, 1, 0)
                df_Location_gam$dif <- c(0, diff(df_Location_gam$cutoff))
                #the column 'cutoff' indicates whether the gam prediction is above (1) or below (0) the SnowFraction threshold
                #the column 'dif' indicates when there is a change from 1 to 0 (-1) or 0 to 1 (1) in the column cutoff
                #Thus, those rows where the column 'dif' is equal to -1 indicate moments where the SnowFraction value changes from above
                #the threshold to below the threshold. It might be possible that this happens multiple times within a season due to
                #measurement errors or cloud effects. We therefore need to determine which 'cutoff' most likely corresponds to the
                #actual moment of snowmelt

                #If SnowFraction-threshold was at least crossed once:
                if(any(df_Location_gam$dif<0)){

                  #Select all moments (cutoffs) where dif==-1
                  cutoffs <- data.frame(index=which(df_Location_gam$dif<0))

                  #For the period 30 days after each cutoff point, sum the number of days that have a SnowFraction value larger than Snowfraction_threshold. If a
                  #cutoff represents actual snowmelt, then we do not expect any days after this moment with SnowFraction > Snowfraction_threshold. Thus, the closer
                  #this sum is to 0, the more likely this cutoff corresponds to the actual moment of snowmelt.
                  cutoffs$min <- cutoffs$index -30
                  cutoffs$min[cutoffs$min<1] <- 1
                  cutoffs$max <- cutoffs$index + 29
                  cutoffs$max[cutoffs$max>nrow(df_Location_gam)] <- nrow(df_Location_gam)
                  cutoffs$sum_cutoff_plus_30 <- apply(cutoffs, 1, function(x){sum(df_Location_gam$cutoff[x['index']:(x['max'])])})
                  cutoff_best <- cutoffs[cutoffs$sum_cutoff_plus_30==min(cutoffs$sum_cutoff_plus_30),'index'][1]

                  #Approximate day of snowmelt in period from (cutoff_best-1 : cutoff_best)
                  newdata_subset <- df_Location_gam[max(0, cutoff_best-2) : min(cutoff_best+1, nrow(df_Location_gam)),]
                  doy_snowmelt <- stats::approx(x = newdata_subset$SnowFraction_gam_predict, y = newdata_subset$Date_doy, xout = Snowfraction_threshold)$y[1]

                  # For debugging:
                  # p_tmp <- ggplot()+
                  #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==FALSE,], aes(x=Date_doy, y=SnowFraction))+
                  #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==TRUE,], aes(x=Date_doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
                  #   geom_line(data=df_Location_gam, aes(x=Date_doy, y=SnowFraction_gam_predict), col = "red") +
                  #   geom_point(aes(x=doy_snowmelt, y=Snowfraction_threshold), col="blue", size=3)+
                  #   geom_hline(yintercept=Snowfraction_threshold, lty=2, col="grey")+
                  #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                  #   ylab("SnowFraction-value at pixel") +
                  #   theme_classic()+
                  #   ggtitle(i)

                }

                #If threshold was not crossed:
                if(!any(df_Location_gam$dif<0)){

                  #No date of snowmelt could be defined
                  doy_snowmelt <- NA

                  #For debugging:
                  # p_tmp <- ggplot()+
                  #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==FALSE,], aes(x=Date_doy, y=SnowFraction))+
                  #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==TRUE,], aes(x=Date_doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
                  #   geom_line(data=df_Location_gam, aes(x=Date_doy, y=SnowFraction_gam_predict), col = "red") +
                  #   geom_hline(yintercept=Snowfraction_threshold, lty=2, col="grey")+
                  #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                  #   ylab("SnowFraction-value at pixel") +
                  #   theme_classic()+
                  #   ggtitle(i)

                }

               #Add day of snowmelt and LocationID to dataframe df_Snowmelt_SnowFraction_Locations
                df_Snowmelt_SnowFraction_Locations[which(unique(Locations_SnowFraction_predictions$LocationID)==i),'LocationID'] <- as.character(i)
                df_Snowmelt_SnowFraction_Locations[which(unique(Locations_SnowFraction_predictions$LocationID)==i),'Date_doy'] <- doy_snowmelt
               }

            #Plot SnowFraction, model predictions and date of snowmelt in a separate plot per location
             p_Locations_SnowFraction_Snowmelt_grid <- p_Locations_SnowFraction_grid +
               geom_point(data=df_Snowmelt_SnowFraction_Locations[!is.na(df_Snowmelt_SnowFraction_Locations$Date_doy),], aes(x=Date_doy, y=Snowfraction_threshold), col="red", size=3)

               ggsave(plot=p_Locations_SnowFraction_Snowmelt_grid, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_SnowFraction_Snowmelt_grid.pdf"), width=12, height = 10)

            #Add coordinates to each location
             df_Snowmelt_SnowFraction_Locations <- left_join(df_Snowmelt_SnowFraction_Locations, Locations_sf, by=c("LocationID"))

            #Save date of snowmelt per location as a csv file
             write.csv(df_Snowmelt_SnowFraction_Locations, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Snowmelt_SnowFraction_", Snowfraction_threshold, ".csv"), row.names = FALSE)

    } 

##################################################################################################################################
        
#IV: Fit a Generalized Additive Model (GAM) through the Sentinel-2 band value data for each location
        
##################################################################################################################################       
 
 #Method "avg_NDSI"            
 if("avg_NDSI" %in% method){       
               
     #(22): SNOW - Fit a Generalized Additive Model (GAM) through the NDSI data for each Location

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

         #(a) Specify at which NDSI value, the point is considered snow-free
             NDSI_threshold = NDSI_threshold

         #(b) Create an empty dataframe
             Locations_NDSI <- data.frame(NDSI=numeric(),
                                          Date=factor(),
                                          Date_doy=numeric(),
                                          LocationID=factor(),
                                          outliers=logical())

             Locations_NDSI_predictions <- data.frame(LocationID=character(),
                                                      Date_doy=numeric(),
                                                      NDSI_gam_predict=numeric(),
                                                      stringsAsFactors=FALSE)

         #(c) Loop through all Locations and fit a separate gam with sequential outlier removal to the location specific mean NDSI data
             for(i in unique(Locations_Data_BandValues$LocationID)){

               #For debugging
                #i=unique(Locations_Data_BandValues$LocationID)[1]

               #Select Location-specific subset of data:
                df_Location <- Locations_Data_BandValues[Locations_Data_BandValues$LocationID==i & !is.na(Locations_Data_BandValues$NDSI),
                                                    c("NDSI", "Date", "Date_doy", "LocationID")]

               #Fit a gam through the Location-specific mean_NDSI ~ Date_doy data and employ sequential outlier removal
                df_Location <- f_gam_SeqRemOutliers(data=df_Location, y="NDSI", x="Date_doy", outlier_removal=outlier_removal,
                                                    outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                    default_k=gam_k_outlier)

               #Sort df_Location by Date_doy:
                df_Location <- df_Location[order(df_Location$Date_doy),]

               #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
                Locations_NDSI <- rbind(Locations_NDSI, df_Location)

               #Create more detailed predictions (not only at the Date_doy present in the datframe) at a 1-day interval to plot more smooth curves

                #Refit GAM through data
                 index <- which(df_Location$outliers==FALSE)
                 mod_gam <- with(df_Location[index,], mgcv::gam(NDSI ~ s(Date_doy, k=min(gam_k, length(index)-1)), method="REML"))

                #Use gam to make predictions on a more detailed (1-day) day of year interval
                 aoi_NDSI_predicted <- data.frame(LocationID=i, Date_doy=seq(min(df_Location$Date_doy), max(df_Location$Date_doy), 1))
                 aoi_NDSI_predicted$NDSI_gam_predict <- stats::predict(mod_gam, newdata=aoi_NDSI_predicted, type="response")
                 aoi_NDSI_predicted <- aoi_NDSI_predicted[order(aoi_NDSI_predicted$Date_doy),]

                #Add predictions to Locations_NDSI_predictions dataframe:
                 Locations_NDSI_predictions <- rbind(Locations_NDSI_predictions, aoi_NDSI_predicted)

               }

             #Change column LocationID to a factor:
              Locations_NDSI$LocationID <- as.factor(as.character(Locations_NDSI$LocationID))
              Locations_NDSI_predictions$LocationID <- as.factor(as.character(Locations_NDSI_predictions$LocationID))

         #(d) Plot the raw NDSI datapoints and gam predictions for each Location:

             #Plot NDSI and model predictions for all locations in a single plot
              p_Locations_NDSI = ggplot()+
                geom_point(data=Locations_NDSI, aes(x=Date_doy, y=NDSI, fill=LocationID, col=LocationID))+
                geom_line(data=Locations_NDSI_predictions, aes(x=Date_doy, y=NDSI_gam_predict, col=LocationID)) +
                xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
                ylab("NDSI") +
                theme_tom()

                #ggsave(plot=p_Locations_NDSI, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDSI.pdf"), width=10, height = 8)

             #Plot NDSI and model predictions in a separate plot per location
              p_Locations_NDSI_grid = ggplot()+
                geom_point(data=Locations_NDSI[Locations_NDSI$outliers==FALSE,], aes(x=Date_doy, y=NDSI))+
                geom_point(data=Locations_NDSI[Locations_NDSI$outliers==TRUE,], aes(x=Date_doy, y=NDSI), col="black", pch=16, alpha=0.2)+
                geom_line(data=Locations_NDSI_predictions, aes(x=Date_doy, y=NDSI_gam_predict), col="#1620de", lwd=1.25)+
                geom_vline(xintercept = 150, colour="grey", lty=2)+
                geom_hline(yintercept = NDSI_threshold, colour="grey", lty=2)+
                facet_wrap(~LocationID, ncol=ceiling(length(unique(Locations_NDSI$LocationID))^0.5))+
                xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
                ylab("NDSI") +
                theme_tom()

                ggsave(plot=p_Locations_NDSI_grid, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDSI_grid.pdf"), width=12, height = 10)

         #(e) Calculate at which day of year the average NDSI is equal to NDSI_threshold for each location using predictions from mod_gam

             #Create an empty dataframe
              df_Snowmelt_NDSI_Locations <- data.frame(LocationID=character(),
                                              Date_doy=numeric(),
                                              stringsAsFactors=FALSE)

             #Loop through all locations in the Locations_NDSI_predictions dataframe and extract doy of snowmelt for each location:
              for(i in unique(Locations_NDSI_predictions$LocationID)){

                #For debugging
                 #i=unique(Locations_NDSI_predictions$LocationID)[1] #for debugging

                #Select location-specific subset of data:
                 df_Location_gam <- Locations_NDSI_predictions[Locations_NDSI_predictions$LocationID==i & !is.na(Locations_NDSI_predictions$NDSI),]
                 df_Location_tmp <- Locations_NDSI[Locations_NDSI$LocationID==i & !is.na(Locations_NDSI$NDSI),]

                 #Detect cutoff points where curve goes from above NDSI threshold to below NDSI threshold
                 df_Location_gam$cutoff <- ifelse(df_Location_gam$NDSI_gam_predict >= NDSI_threshold, 1, 0)
                 df_Location_gam$dif <- c(0, diff(df_Location_gam$cutoff))
                 #the column 'cutoff' indicates whether the gam prediction is above (1) or below (0) the ndsi threshold
                 #the column 'dif' indicates when there is a change from 1 to 0 (-1) or 0 to 1 (1) in the column cutoff
                 #Thus, those rows where the column 'dif' is equal to -1 indicate moments where the NDSI value changes from above
                 #the threshold to below the threshold. It might be possible that this happens multiple times within a season due to
                 #measurement errors or cloud effects. We therefore need to determine which 'cutoff' most likely corresponds to the
                 #actual moment of snowmelt

                 #If NDSI-threshold was at least crossed once:
                 if(any(df_Location_gam$dif<0)){

                   #Select all moments (cutoffs) where dif==-1
                   cutoffs <- data.frame(index=which(df_Location_gam$dif<0))

                   #For the period 30 days after each cutoff point, sum the number of days that have a NDSI value larger than NDSI_threshold. If a
                   #cutoff represents actual snowmelt, then we do not expect any days after this moment with NDSI > NDSI_threshold. Thus, the closer
                   #this sum is to 0, the more likely this cutoff corresponds to the actual moment of snowmelt.
                   cutoffs$min <- cutoffs$index -30
                   cutoffs$min[cutoffs$min<1] <- 1
                   cutoffs$max <- cutoffs$index + 29
                   cutoffs$max[cutoffs$max>nrow(df_Location_gam)] <- nrow(df_Location_gam)
                   cutoffs$sum_cutoff_plus_30 <- apply(cutoffs, 1, function(x){sum(df_Location_gam$cutoff[x['index']:(x['max'])])})
                   cutoff_best <- cutoffs[cutoffs$sum_cutoff_plus_30==min(cutoffs$sum_cutoff_plus_30),'index'][1]

                   #Approximate day of snowmelt in period from (cutoff_best-1 : cutoff_best)
                   newdata_subset <- df_Location_gam[max(0, cutoff_best-2) : min(cutoff_best+1, nrow(df_Location_gam)),]
                   doy_snowmelt <- stats::approx(x = newdata_subset$NDSI_gam_predict, y = newdata_subset$Date_doy, xout = NDSI_threshold)$y[1]

                   #For debugging:
                   # p_tmp <- ggplot()+
                   #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==FALSE,], aes(x=Date_doy, y=NDSI))+
                   #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==TRUE,], aes(x=Date_doy, y=NDSI), col="black", pch=16, alpha=0.2)+
                   #   geom_line(data=df_Location_gam, aes(x=Date_doy, y=NDSI_gam_predict), col = "red") +
                   #   geom_point(aes(x=doy_snowmelt, y=NDSI_threshold), col="blue", size=3)+
                   #   geom_hline(yintercept=NDSI_threshold, lty=2, col="grey")+
                   #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                   #   ylab("NDSI-value at pixel") +
                   #   theme_classic()+
                   #   ggtitle(i)

                 }

                 #If threshold was not crossed:
                 if(!any(df_Location_gam$dif<0)){

                   #No date of snowmelt could be defined
                   doy_snowmelt <- NA

                   #For debugging:
                   # p_tmp <- ggplot()+
                   #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==FALSE,], aes(x=Date_doy, y=NDSI))+
                   #   geom_point(data=df_Location_tmp[df_Location_tmp$outliers==TRUE,], aes(x=Date_doy, y=NDSI), col="black", pch=16, alpha=0.2)+
                   #   geom_line(data=df_Location_gam, aes(x=Date_doy, y=NDSI_gam_predict), col = "red") +
                   #   geom_hline(yintercept=NDSI_threshold, lty=2, col="grey")+
                   #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                   #   ylab("NDSI-value at pixel") +
                   #   theme_classic()+
                   #   ggtitle(i)

                 }

                #Add day of snowmelt and LocationID to dataframe df_Snowmelt_NDSI_Locations
                 df_Snowmelt_NDSI_Locations[which(unique(Locations_NDSI_predictions$LocationID)==i),'LocationID'] <- as.character(i)
                 df_Snowmelt_NDSI_Locations[which(unique(Locations_NDSI_predictions$LocationID)==i),'Date_doy'] <- doy_snowmelt
                }

             #Plot NDSI, model predictions and date of snowmelt in a separate plot per location
              p_Locations_NDSI_Snowmelt_grid <- p_Locations_NDSI_grid +
                geom_point(data=df_Snowmelt_NDSI_Locations[!is.na(df_Snowmelt_NDSI_Locations$Date_doy),], aes(x=Date_doy, y=NDSI_threshold), col="red", size=3)

                ggsave(plot=p_Locations_NDSI_Snowmelt_grid, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDSI_Snowmelt_grid.pdf"), width=12, height = 10)

             #Add coordinates to each location
              df_Snowmelt_NDSI_Locations <- left_join(df_Snowmelt_NDSI_Locations, Locations_sf, by=c("LocationID"))

             #Save date of snowmelt per location as a csv file
              write.csv(df_Snowmelt_NDSI_Locations, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Snowmelt_NDSI_", NDSI_threshold, ".csv"), row.names = FALSE)

              
      #(23): VEGETATION - Fit a Generalized Additive Model (GAM) through the NDVI data for each location
              
         #(a) Create an empty dataframe
             Locations_NDVI <- data.frame(NDVI=numeric(),
                                          Date=factor(),
                                          Date_doy=numeric(),
                                          LocationID=factor(),
                                          outliers=logical())
             
             Locations_NDVI_predictions <- data.frame(LocationID=character(), 
                                                      Date_doy=numeric(), 
                                                      NDVI_gam_predict=numeric(), 
                                                      stringsAsFactors=FALSE)   
             
         #(b) Loop through all Locations and fit a separate gam with sequential outlier removal to the Location specific NDVI data
             for(i in unique(Locations_Data_BandValues$LocationID)){
               
               #For debugging  
               #i=unique(Locations_Data_BandValues$LocationID)[1]
               
               #Select Location-specific subset of data:
               df_Location <- Locations_Data_BandValues[Locations_Data_BandValues$LocationID==i & !is.na(Locations_Data_BandValues$NDVI), 
                                                   c("NDVI", "Date", "Date_doy", "LocationID")]
               
               #Fit a gam through the location-specific NDVI ~ Date_doy data and emply sequential outlier removal
               df_Location <- f_gam_SeqRemOutliers(data=df_Location, y="NDVI", x="Date_doy", outlier_removal=outlier_removal, 
                                               outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                               default_k=gam_k_outlier)
               
               #Sort df_Location by Date_doy
               df_Location <- df_Location[order(df_Location$Date_doy), ]
               
               #Bind the location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
               Locations_NDVI <- rbind(Locations_NDVI, df_Location)
               
               #Create more detailed predictions (not only at the Date_doy present in the datframe) at a 1-day interval to plot more smooth curves
               
               #Refit GAM through data
               index <- which(df_Location$outliers==FALSE)
               mod_gam <- with(df_Location[index,], mgcv::gam(NDVI ~ s(Date_doy, k=min(gam_k, length(index)-1)), method="REML"))
               
               #Use gam to make predictions on a more detailed (1-day) day of year interval
               aoi_NDVI_predicted <- data.frame(LocationID=i, Date_doy=seq(min(df_Location$Date_doy), max(df_Location$Date_doy), 1))
               aoi_NDVI_predicted$NDVI_gam_predict <- stats::predict(mod_gam, newdata=aoi_NDVI_predicted, type="response")
               aoi_NDVI_predicted <- aoi_NDVI_predicted[order(aoi_NDVI_predicted$Date_doy), ]
               
               #Add predictions to Locations_NDVI_predictions dataframe:  
               Locations_NDVI_predictions <- rbind(Locations_NDVI_predictions, aoi_NDVI_predicted)
               
             }
             
             #Change column LocationID to a factor:
             Locations_NDVI$LocationID <- as.factor(as.character(Locations_NDVI$LocationID))
             Locations_NDVI_predictions$LocationID <- as.factor(as.character(Locations_NDVI_predictions$LocationID))
             
         #(c) Plot the raw NDVI datapoints and gam predictions for each location:
              
             #Plot NDVI and model predictions for all locations in a single plot
             p_Locations_NDVI = ggplot()+ 
                geom_point(data=Locations_NDVI, aes(x=Date_doy, y=NDVI, fill=LocationID, col=LocationID))+
                geom_line(data=Locations_NDVI_predictions, aes(x=Date_doy, y=NDVI_gam_predict, col=LocationID)) + 
                xlab("Day of year") + 
                ylab("NDVI") +
                theme_tom()
             
                #ggsave(plot=p_Locations_NDVI, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDVI.pdf"), width=10, height = 8)
              
             #Plot NDVI and model predictions in a separate plot per location
             p_Locations_NDVI_grid = ggplot()+ 
                geom_point(data=Locations_NDVI[Locations_NDVI$outliers==FALSE,], aes(x=Date_doy, y=NDVI), col="black")+
                geom_point(data=Locations_NDVI[Locations_NDVI$outliers==TRUE,], aes(x=Date_doy, y=NDVI), col="black", pch=16, alpha=0.2)+
                geom_line(data=Locations_NDVI_predictions, aes(x=Date_doy, y=NDVI_gam_predict), col="#08a31a", lwd=1.25) + 
                geom_vline(xintercept = 150, colour="grey", lty=2)+
                facet_wrap(~LocationID, ncol=ceiling(length(unique(Locations_NDVI$LocationID))^0.5))+
                xlab("Day of year") + 
                ylab("NDVI") +
                theme_tom()
             
                ggsave(plot=p_Locations_NDVI_grid, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDVI_grid.pdf"), width=12, height = 10)
             
      #(24): MOISTURE - Fit a Generalized Additive Model (GAM) through the NDMI data for each location
              
         #(a) Create an empty dataframe
             Locations_NDMI <- data.frame(NDMI=numeric(),
                                          Date=factor(),
                                          Date_doy=numeric(),
                                          LocationID=factor(),
                                          outliers=logical())
             
             Locations_NDMI_predictions <- data.frame(LocationID=character(), 
                                                      Date_doy=numeric(), 
                                                      NDMI_gam_predict=numeric(), 
                                                      stringsAsFactors=FALSE)   
             
         #(b) Loop through all locations and fit a separate gam with sequential outlier removal to the location specific NDMI data
             for(i in unique(Locations_Data_BandValues$LocationID)){
               
               #For debugging  
               #i=unique(Locations_Data_BandValues$LocationID)[1]
               
               #Select location-specific subset of data:
               df_Location <- Locations_Data_BandValues[Locations_Data_BandValues$LocationID==i & !is.na(Locations_Data_BandValues$NDMI), 
                                                   c("NDMI", "Date", "Date_doy", "LocationID")]
               
               #Fit a gam through the locations-specific NDMI ~ Date_doy data and emply sequential outlier removal
               df_Location <- f_gam_SeqRemOutliers(data=df_Location, y="NDMI", x="Date_doy", outlier_removal=outlier_removal, 
                                                   outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                   default_k=gam_k_outlier)
               
               #Sort df_Location by Date_doy
               df_Location <- df_Location[order(df_Location$Date_doy), ]
               
               #Bind the locations-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
               Locations_NDMI <- rbind(Locations_NDMI, df_Location)
               
               #Create more detailed predictions (not only at the Date_doy present in the datframe) at a 1-day interval to plot more smooth curves
               
               #Refit GAM through data
               index <- which(df_Location$outliers==FALSE)
               mod_gam <- with(df_Location[index,], mgcv::gam(NDMI ~ s(Date_doy, k=min(gam_k, length(index)-1)), method="REML"))
               
               #Use gam to make predictions on a more detailed (1-day) day of year interval
               aoi_NDMI_predicted <- data.frame(LocationID=i, Date_doy=seq(min(df_Location$Date_doy), max(df_Location$Date_doy), 1))
               aoi_NDMI_predicted$NDMI_gam_predict <- stats::predict(mod_gam, newdata=aoi_NDMI_predicted, type="response")
               aoi_NDMI_predicted <- aoi_NDMI_predicted[order(aoi_NDMI_predicted$Date_doy), ]
               
               #Add predictions to Locations_NDMI_predictions dataframe:  
               Locations_NDMI_predictions <- rbind(Locations_NDMI_predictions, aoi_NDMI_predicted)
               
             }
             
              #Change column LocationID to a factor:
               Locations_NDMI$LocationID <- as.factor(as.character(Locations_NDMI$LocationID))
               Locations_NDMI_predictions$LocationID <- as.factor(as.character(Locations_NDMI_predictions$LocationID)) 
             
         #(c) Plot the raw NDMI datapoints and gam predictions for each location:
              
             #Plot NDMI and model predictions for all locations in a single plot
              p_Locations_NDMI = ggplot()+ 
                geom_point(data=Locations_NDMI, aes(x=Date_doy, y=NDMI, fill=LocationID, col=LocationID))+
                geom_line(data=Locations_NDMI_predictions, aes(x=Date_doy, y=NDMI_gam_predict, col=LocationID)) + 
                xlab("Day of year") + 
                ylab("NDMI") +
                theme_tom()
              
                #ggsave(plot=p_Locations_NDMI, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDMI.pdf"), width=12, height = 10)
              
             #Plot NDMI and model predictions in a separate plot per location
              p_Locations_NDMI_grid = ggplot()+ 
                geom_point(data=Locations_NDMI[Locations_NDMI$outliers==FALSE,], aes(x=Date_doy, y=NDMI), col="black")+
                geom_point(data=Locations_NDMI[Locations_NDMI$outliers==TRUE,], aes(x=Date_doy, y=NDMI), col="black", pch=16, alpha=0.2)+
                geom_line(data=Locations_NDMI_predictions, aes(x=Date_doy, y=NDMI_gam_predict), col="#16acde", lwd=1.25)+
                geom_vline(xintercept = 150, colour="grey", lty=2)+
                facet_wrap(~LocationID, ncol=ceiling(length(unique(Locations_NDMI$LocationID))^0.5))+
                xlab("Day of year") + 
                ylab("NDMI") +
                theme_tom()
              
                ggsave(plot=p_Locations_NDMI_grid, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_NDMI_grid.pdf"), width=12, height = 10)
              
      #(25): Plot NDSI, NDVI and NDMI together in a plot per location
              p_Locations_BANDS_grid = ggplot()+ 
                geom_point(data=Locations_NDSI[Locations_NDSI$outliers==FALSE,], aes(x=Date_doy, y=NDSI))+
                geom_point(data=Locations_NDSI[Locations_NDSI$outliers==TRUE,], aes(x=Date_doy, y=NDSI), col="black", pch=16, alpha=0.2)+
                geom_point(data=Locations_NDVI[Locations_NDVI$outliers==FALSE,], aes(x=Date_doy, y=NDVI), col="black")+
                geom_point(data=Locations_NDVI[Locations_NDVI$outliers==TRUE,], aes(x=Date_doy, y=NDVI), col="black", pch=16, alpha=0.2)+
                geom_point(data=Locations_NDMI[Locations_NDMI$outliers==FALSE,], aes(x=Date_doy, y=NDMI), col="black")+
                geom_point(data=Locations_NDMI[Locations_NDMI$outliers==TRUE,], aes(x=Date_doy, y=NDMI), col="black", pch=16, alpha=0.2)+
                geom_line(data=Locations_NDSI_predictions, aes(x=Date_doy, y=NDSI_gam_predict), col="#1620de", lwd=1.25)+
                geom_line(data=Locations_NDVI_predictions, aes(x=Date_doy, y=NDVI_gam_predict), col="#08a31a", lwd=1.25)+ 
                geom_line(data=Locations_NDMI_predictions, aes(x=Date_doy, y=NDMI_gam_predict), col="#16acde", lwd=1.25)+
                geom_vline(xintercept = 150, colour="grey", lty=2)+
                facet_wrap(~LocationID, ncol=ceiling(length(unique(Locations_NDMI$LocationID))^0.5))+
                xlab("Day of year") + 
                ylab("Normalized Difference Band Index") +
                theme_tom()
              
                ggsave(plot=p_Locations_BANDS_grid, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_Locations_Plot_AllBands_grid.pdf"), width=14, height =12)
 }  
      
####################################################################################################################################################
              
#X: OPTIONAL: Assess how the extracted day of snowmelt, NDVI and NDMI values change along the foraging trajectories (Locations) of chicks in ZAC            
              
####################################################################################################################################################              
              
  # #Add Date of snowmelt, NDSI, NDMI and NDVI measures per location to the corresponding date-location combination in ZAC_Paths
  #   ZAC_Paths <- read.csv(paste0("Input/", input_locations), header=T)[,c("ChickID", "LON_x", "LAT_y", "DateTime")]
  #   colnames(ZAC_Paths) <- c("ChickID", "LON_x", "LAT_y", "DateTime")
  #   #Add a unique LocationID to every unique lat/lon combination
  #   Locations <- unique(ZAC_Paths[,c("LON_x", "LAT_y")])
  #   Locations$LocationID <- paste0("Location_", 1:nrow(Locations))
  #   ZAC_Paths <- left_join(ZAC_Paths, Locations, by=c("LON_x", "LAT_y"))
  #   ZAC_Paths$Date_doy <- as.numeric(strftime(ZAC_Paths$DateTime , format = "%j"))
  #   ZAC_Paths <- left_join(ZAC_Paths, Locations_NDSI_predictions, by=c("LocationID", "Date_doy"))
  #   ZAC_Paths <- left_join(ZAC_Paths, Locations_NDVI_predictions, by=c("LocationID", "Date_doy"))
  #   ZAC_Paths <- left_join(ZAC_Paths, Locations_NDMI_predictions, by=c("LocationID", "Date_doy"))
  #   ZAC_Paths <- left_join(ZAC_Paths, df_Snowmelt_NDSI_Locations[,c("LocationID", "Date_doy")], by=c("LocationID"))
  #   colnames(ZAC_Paths) <- c("ChickID", "LON_x", "LAT_y", "DateTime", "LocationID", "Date_doy", "NDSI_gam_predict", "NDVI_gam_predict",
  #                            "NDMI_gam_predict", "Snowmelt_doy")
  #   write.csv(ZAC_Paths, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_ChickMovement_Data_NDSI_NDVI_NDMI.csv"), row.names = FALSE)
  # 
  # #Plot each individuals path through a landscape of Snowmelt doy values
  #   p1 <- ggplot()+
  #     geom_line(data=ZAC_Paths, aes(x=Date_doy, y=Snowmelt_doy ), col="#1620de", lwd=1.25)+
  #     geom_abline()+
  #     facet_wrap(~ChickID, ncol=ceiling(length(unique(ZAC_Paths$ChickID))^0.5))+
  #     geom_smooth(data=ZAC_Paths, aes(x=Date_doy, y=Snowmelt_doy), method="lm", formula=y~x, col="red", size=1.25, se=T)+
  #     geom_point(data=ZAC_Paths, aes(x=Date_doy, y=Snowmelt_doy ), col="black")+
  #     xlab("Day of year") +
  #     ylab("Day of snowmelt") +
  #     theme_tom()
  # 
  #     ggsave(plot=p1, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_MovementPaths_Plot_Snowmelt_doy_grid.pdf"), width=10, height = 8)
  # 
  # #Plot each individuals path through a landscape of NDVI values
  #   p2 <- ggplot()+
  #     geom_point(data=ZAC_Paths, aes(x=Date_doy, y=NDVI_gam_predict ), col="black")+
  #     geom_line(data=ZAC_Paths, aes(x=Date_doy, y=NDVI_gam_predict ), col="#08a31a", lwd=1.25)+
  #     facet_wrap(~ChickID, ncol=ceiling(length(unique(ZAC_Paths$ChickID))^0.5))+
  #     geom_smooth(data=ZAC_Paths, aes(x=Date_doy, y=NDVI_gam_predict), method="lm", formula=y~x, col="red", size=1.25, se=T)+
  #     xlab("Day of year") +
  #     ylab("NDVI") +
  #     theme_tom()
  # 
  #     ggsave(plot=p2, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_MovementPaths_Plot_NDVI_grid.pdf"), width=10, height = 8)
  # 
  # #Plot each individuals path through a landscape of NDMI values
  #   p3 <- ggplot()+
  #     geom_point(data=ZAC_Paths, aes(x=Date_doy, y=NDMI_gam_predict ), col="black")+
  #     geom_line(data=ZAC_Paths, aes(x=Date_doy, y=NDMI_gam_predict ), col="#16acde", lwd=1.25)+
  #     facet_wrap(~ChickID, ncol=ceiling(length(unique(ZAC_Paths$ChickID))^0.5))+
  #     geom_smooth(data=ZAC_Paths, aes(x=Date_doy, y=NDMI_gam_predict), method="lm", formula=y~x, col="red", size=1.25, se=T)+
  #     xlab("Day of year") +
  #     ylab("NDMI") +
  #     theme_tom()
  # 
  #     ggsave(plot=p3, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_MovementPaths_Plot_NDMI_grid.pdf"), width=10, height = 8)
  # 
  # #Calculate the average encountered band value for each individual:
  #  #There do not seem to be any clear trajectories in individuals moving up to higher NDVI or NDMI values. We therefore calculate
  #  #the average encountered day of snowmelt, NDVI and NDMI values along the foraging trajectory
  #  ZAC_Chicks_MeanBandvalue <- ZAC_Paths %>% group_by(ChickID) %>% summarise(mean_Snowmelt_doy=mean(Snowmelt_doy),
  #                                                                            mean_NDVI=mean(NDVI_gam_predict),
  #                                                                            mean_NDMI=mean(NDMI_gam_predict),
  #                                                                            mean_doy=mean(Date_doy))
  # 
  # #Plot average encountered band-values along growth parameter of each chick
  # 
  #   #Load estimated growth parameters per chick
  #    ZAC_ChickGrowth <- read.csv(paste0("Input/ZAC19_ChickGrowth.csv"), header=T)
  #    ZAC_ChickGrowth <- ZAC_ChickGrowth[ZAC_ChickGrowth$max_age>=5,]
  # 
  #   #Combine growth parameters with mean Band value during foraging
  #    ZAC_Chicks_MeanBandvalue <- left_join(ZAC_Chicks_MeanBandvalue, ZAC_ChickGrowth, by="ChickID")
  # 
  #   #Plot average encountered day of snowmelt, NDVI and NDMI versus the upper growth asymptote per chick
  # 
  #     #Mean Snowmelt_doy vs upper asymptote
  #      p4 <- ggplot()+
  #        geom_smooth(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_Snowmelt_doy, y=d, size=max_age), method="lm", formula=y~x, col="red", size=1.25)+
  #        geom_point(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_Snowmelt_doy, y=d, size=max_age), col="black")+
  #        xlab("Mean day of snowmelt across foraging locations") +
  #        ylab("Upper asymptote bodymass (gram)") +
  #        theme_tom()
  # 
  #        ggsave(plot=p4, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_GrowthParameter_Plot_Mean_Snowmelt_doy.pdf"), width=10, height = 8)
  # 
  #     #Mean NDVI vs upper asymptote
  #      p5 <- ggplot()+
  #        geom_point(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_NDVI, y=d, size=max_age), col="black")+
  #        geom_smooth(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_NDVI, y=d, size=max_age), method="lm", formula=y~x, col="red", size=1.25)+
  #        xlab("Mean NDVI across foraging locations") +
  #        ylab("Upper asymptote bodymass (gram)") +
  #        theme_tom()
  # 
  #        ggsave(plot=p5, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_GrowthParameter_Plot_MeanNDVI.pdf"), width=10, height = 8)
  # 
  #     #Mean NDMI vs upper asymptote
  #      p6 <- ggplot()+
  #        geom_point(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_NDMI, y=d, size=max_age), col="black")+
  #        geom_smooth(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_NDMI, y=d, size=max_age), method="lm", formula=y~x, col="red", size=1.25)+
  #        xlab("Mean NDMI across foraging locations") +
  #        ylab("Upper asymptote bodymass (gram)") +
  #        theme_tom()
  # 
  #        ggsave(plot=p6, paste0("Output/S2/Shapefile_Points_Snowmelt/", data_ID, "_Buffer", Buffer_radius_m, "_Resolution", resolution, "_GrowthParameter_Plot_MeanNDMI.pdf"), width=10, height = 8)
  # 
  #  #Construct linear mixed models to test if growth is related to average encountered band values
  #     lmm0 <- lme(d ~ 1, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     lmm1 <- lme(d ~ mean_Snowmelt_doy, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     lmm2 <- lme(d ~ mean_Snowmelt_doy + I(mean_Snowmelt_doy^2), random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     MuMIn::AICc(lmm0, lmm1, lmm2) #lmm0 is best supported
  #     rm(lmm0, lmm1, lmm2)
  # 
  #     lmm0 <- lme(d ~ 1, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     lmm1 <- lme(d ~ mean_NDVI, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     lmm2 <- lme(d ~ mean_NDVI + I(mean_NDVI^2), random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     MuMIn::AICc(lmm0, lmm1, lmm2) #lmm1 is best supported
  #     rm(lmm0, lmm1, lmm2)
  # 
  #     lmm0 <- lme(d ~ 1, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     lmm1 <- lme(d ~ mean_NDMI, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     lmm2 <- lme(d ~ mean_NDMI + I(mean_NDMI^2), random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #     MuMIn::AICc(lmm0, lmm1, lmm2) #lmm0 is best supported
  #     rm(lmm0, lmm1, lmm2)
  # 
  #     #Calculate the difference between the mean visiting date vs mean date of snowmelt for all chicks
  #     mean(ZAC_Chicks_MeanBandvalue$mean_doy - ZAC_Chicks_MeanBandvalue$mean_Snowmelt_doy)
  #     #Chicks on average forage at a certain location 58 days after snowmelt!
  
###############################################################################################################################################
###############################################################################################################################################        
                   