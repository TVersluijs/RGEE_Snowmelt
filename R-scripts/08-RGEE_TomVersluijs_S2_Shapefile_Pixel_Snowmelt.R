##################################################################################################################################

#Use Sentinel-2 data to create pixel-level maps (10m resolution) of the timing of snowmelt in for an area smaller than c.a. 50km2. 
#Snowmelt is calculated per pixel by fitting a GAM through the average NDSI data and extracting the moment this GAM crosses a user 
#specified NDSI threshold. This script requires a single shapefile of the study area as input.

#Copyright Tom Versluijs 2023-11-01. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

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
       p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr, foreach, parallel, doSNOW, gridExtra)

      #(2): Define ggplot2 plotting theme
       theme_tom <- function(){
         theme_classic() %+replace%
         theme(axis.title = element_text(size=18),
             axis.text = element_text(size=16),
             legend.position = "none",
             strip.background = element_rect(fill = "white", colour = "black", size = rel(2)),
             complete = TRUE)}
       
      #(3): Load auxiliary functions
       source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "Sentinel2_AuxiliaryFunctions")
       sapply(source_files, source, chdir = TRUE) ; rm(source_files)
       
      #(4): Initialize earth engine and google drive
       #ee_clean_pyenv()
       #ee_install()
       #ee_check()
       ##rgee:::ee_create_credentials_drive(email="tom.versluijs@gmail.com", quiet=F)
       ##rgee:::ee_create_credentials_earthengine(email="tom.versluijs@gmail.com", display=T)
       ##drive_auth(email = "tom.versluijs@gmail.com")
       rgee::ee_Initialize(user = "tom.versluijs@gmail.com", drive = TRUE)
       #if this results in a credentials error, then delete the folder "C:\Users\tomve\.config\earthengine" and rerun the whole installation.
       
       
##################################################################################################################################
       
#II: Specify parameters of interest
       
##################################################################################################################################       
       
 #(4): Manually specify parameters of interest

   #(a): Sentinel-2 satellite

     #Sentinel-2 dataset ("S2_HARMONIZED" or "S2_SR_HARMONIZED")
     #Use S2_HARMONIZED (level 1c product) for data from 2016 - 2018
     #Use S2_SR_HARMONIZED (level 2a product) for data starting from 2019
     s2_dataset <- "S2_SR_HARMONIZED"

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
     crs <- "EPSG:32627"

   #(c) Dates

     #Specify the year of interest:
     year_ID <- "2022"

     #Date range of all images considered for analysis
     start_date <- paste0(year_ID, "-03-15") #choose date (well) before the first snowmelt occurs within the study site
     end_date <- paste0(year_ID, "-09-15") #choose date (well) after last date of tracking

   #(d) Snow detection

     #NDSI threshold for snow detection
     NDSI_threshold=0.42

   #(e): Cloud masking

     #Should clouds be masked from the analysis (default=TRUE).
     #Note that cloud masking will automatically be set to TRUE when mask_water="water_mask_Manual"
     #Recommended to set to TRUE.
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

   #(g): GAM sequential outlier filtering

     #parameters for auxiliary function to filter outliers from GAM using a two-step approach
     outlier_removal=TRUE #should sequential outlier removal be employed when fitting GAMs to the data
     outlier_thresh_1=0.4 #first threshold in outlier removal (relative to y-range of data, default=0.4)
     outlier_thresh_2=0.3 #second threshold in outlier removal (relative to y-range of data, default=0.3)
       
     #Specify the degree of GAM smoothing by setting the 'k' parameter
     gam_k_outlier=10 #Number of knots when filtering outliers (default=10)
     gam_k=50 #Number of knots when making model predictions (default=50). 
     #Larger values result in a more precise GAM-fit, but at a cost of computation time

          
##################################################################################################################################
         
#III: Create a unique asset folder
         
##################################################################################################################################
         
    #(5) Create a unique data_ID and asset folder for storing the generated datafiles
     
      #Create a unique data_ID
       data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_S2")
        
      #Create a timestamp variable
       timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
     
      #Store NDSI_threshold as a character (used for naming of outputs)
       NDSI_threshold_char <- gsub("\\.", "_", as.character(NDSI_threshold))  
       
      #Create a unique Asset folder (delete this folder if already present) 
       path_asset <- paste0(ee_get_assethome(), "/", data_ID)
       #tryCatch(ee_manage_assetlist(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
       tryCatch(ee_manage_delete(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
       ee_manage_create(path_asset=path_asset, asset_type="Folder")
       ee_manage_assetlist()
    
      #Set mask_clouds to TRUE if water masking is TRUE and mask_water_type==water_mask_Manual
       if(mask_water==TRUE & (mask_water_type=="water_mask_Manual" | mask_water_type=="both")){mask_clouds=TRUE}       
         
      #Create output folder
       if(dir.exists(paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt"), recursive = TRUE)}
       
##################################################################################################################################
         
#IV: Read and display the unfiltered data
         
##################################################################################################################################

    #(6) Read study area shapefile and convert to a feature collection.
       
       #Read aoi_Shapefile shapefile in root folder
       root_fldr <- here()
       aoi_Shapefile <- st_read(paste0(root_fldr, "/Input/Shapefiles/", shapefile), quiet=T)
       aoi_Shapefile <- st_transform(aoi_Shapefile, crs="EPSG:4326")
       aoi_Shapefile <- sf_as_ee(aoi_Shapefile)
       aoi_Shapefile <- ee$FeatureCollection(aoi_Shapefile)
       #This feature collection will be used to clip the satellite imagery to the study area
       
       #Calculate bounding box for aoi_Shapefile
       aoi <- aoi_Shapefile$geometry()$bounds()
       
       #Calculate central point
       coordinates_point <- aoi_Shapefile$geometry()$centroid()
       
       #Plot shapefile and bounding box (for debugging)
       Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
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
       paste0('Size of study area calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')
       
    #(7) Extract Sentinel-2 satellite data including cloud probability band:
     
      #Specify starting and ending date  
       start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
       end_date_doy <- as.numeric(strftime(end_date, format = "%j"))
     
      #Extract Sentinel-2 Surface Reflectance satellite data
       s2_col <- ee$ImageCollection(paste0('COPERNICUS/', s2_dataset))
       s2_col <- s2_col$
         filterBounds(aoi)$ #For areas that do not span multiple tiles, change aoi to coordinates_point to greatly speed-up code.
         filterDate(start_date, end_date)
    
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

      #Clip all images to the area depicted by the shapefile (aoi_Shapefile)
       s2_col <- s2_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})
      
      #Plot a single clipped image: 
       
        #Select a single image for initial plot:
         image <- s2_col$filterDate(paste0(year_ID, "-07-15"), end_date)$first()
      
        #Add Normalized difference snow index to this image:
         ndsi_s2 <-  image$expression('(B3 - B11) / (B3 + B11)',
                                      list('B11'=image$select('B11'),
                                           'B3'=image$select('B3')))
      
        #Add Normalized difference vegetation index to this image:
         ndvi_s2 <- image$expression('(B8 - B4) / (B8 + B4)',
                                     list('B8'=image$select('B8'),
                                          'B4'=image$select('B4')))
         
        #Add Normalized difference water index to this image:
         ndwi_s2 <- image$expression('(B3 - B8) / (B3 + B8)',
                                     list('B3'=image$select('B3'),
                                          'B8'=image$select('B8')))
  
        #Plot all image Bands
         Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
         Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
         Map$addLayer(ndsi_s2,list(min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
         Map$addLayer(ndvi_s2,list(min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
         Map$addLayer(ndwi_s2,list(min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')
     
    #(8) Add NDSI, NDVI, NDMI and NDWI to the clipped image collection
      
     #Map Normalized Difference functions over image collection
      s2_col <- s2_col$
        map(getNDSI)$
        map(getNDVI)$
        map(getNDMI)$
        map(getNDWI)
    
    #(9) Create a timeseries gif of RGB images for the shapefile area
      
      #Check number of images in collection
      s2_col$size()$getInfo()
      
      #Create a timelapse video of RGB band
      videoArgs <- list(dimensions=200, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
      tryCatch({browseURL(s2_col$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
      
        
##################################################################################################################################
        
#IV: Filter and mask clouds within the image collection
        
##################################################################################################################################
        
   #There are many images that are almost fully covered in dense clouds. We thus first filter and mask clouds in the image collection.
   #We delete all images with a cloud cover fraction >= max_cloud_fraction. In the remaining images we mask all pixels that are 
   #covered by clouds.
    if(mask_clouds==TRUE){
      
      #print message   
      print("Cloud masking = TRUE")  
        
      #(10) Add cloud information within the aoi to the image collection: 

        #Add the cloud fraction within the shapefile area to each separate image by mapping the cloud functions over the image collection
        s2_col <- s2_col$
          #Determine which pixels are clouds using the cloud_probability property
          map(compute_Clouds_cldprb)$
          #Add the fraction of cloud-covered pixels within the shapefile area as image property
          map(Add_CloudFraction)$
          #Add NULL to those images in which cloudfraction could not be calculated
          map(Add_NULL_CloudFraction)$
          #Add date characteristics to each image
          map(add_Date)
        
        # #Check if CloudFraction property has been added to each image
        # s2_col$first()$propertyNames()$getInfo()
        # 
        # #Visual check of cloud mask (for debugging)
        # image <- s2_col$filterDate(paste0(year_ID, "-05-31"), end_date)$first()$
        #   clipToCollection(aoi_Shapefile)$
        #   select("clouds_unadjusted", "clouds", "clouds_inv", "clouds_prob", "B4", "B3", "B2")
        # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
        # Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
        #   #Map$addLayer(image, list(bands='clouds_prob', min=0, max=100, opacity=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'clouds_prob')+
        #   Map$addLayer(image, list(bands='clouds_unadjusted', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'clouds unadjusted')+
        #   Map$addLayer(image, list(bands='clouds', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'clouds buffered')
        #   #Map$addLayer(image, list(bands='clouds_inv', min=0, max=1, opacity=1), 'clouds_inv')
  
        #We have now calculated the fraction of cloud covered pixels within our shapefile area for each image in the image collection.
        #However, Sentinel-2 satellite data is recorded as different tiles, these tiles are generally much bigger than a single study site.
        #Each Sentinel-2 tile also contains a native measure of cloud cover. This is the property 'CLOUDY_PIXEL_PERCENTAGE'. This property
        #thus refers to the percentage of cloud cover within the whole tile. We have now calculated a more appropriate measure for our
        #shapefile area only, called 'CloudFraction'. Below we compare both measures, as we expect they should be correlated. Differences
        #between both measures can for example arise when the study area was covered by clouds while most of the larger region (tile)
        #was not or vice versa.
        
        # #Extract cloud fraction of all images in image collection for the shapefile area
        # cloud_aoi <- unlist(s2_col$aggregate_array('CloudFraction')$getInfo())
        # 
        # #Extract cloud fraction of all images in image collection for the complete tile to which our shapefile area belongs
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
        #   theme_classic()
        #   #Thus, this indicates that our local cloud measure does not necessarily correspond to that as measured on the whole tile as
        #   #explained above.
          
      #(11): Exclude all images from the image collection that have a CloudFraction value >= max_cloud_fraction, and mask all cloud pixels within the remaining images
        
        #In the analysis for the whole study area we exclude images with more than 20% snowcover within our shapefile
        #area However, now we're interested in individual pixels. So even when 90% of the shapefile area is snowcovered, 
        #the pixels we're interested in might be located in the 10% cloud free area. We thus only filter images with
        #a 100% cloudcover (cloudFraction=1.0). For all remaining images we apply a cloud mask to mask out cloudy 
        #pixels within each image.
        
        #Apply a manually constructed cloudmask to each image in the collection:
         s2_clouds_filtered <- s2_col$
           #Filter all images with 100% cloudcover (CloudFraction=1.0)
           filter(ee$Filter$lt('CloudFraction', max_cloud_fraction))$
           #Apply cloudmask for individual pixels
           map(Add_CloudMask)
         
        # #Create timelapse video of the cloud filtered/masked RGB images (for debugging)
        #  videoArgs <- list(dimensions=200, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
        #  tryCatch({browseURL(s2_clouds_filtered$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
        }
    if(mask_clouds==FALSE){
    
      #print message   
      print("Cloud masking = FALSE")
      
      #Add Date characteristics to each image
      s2_clouds_filtered <- s2_col$
        map(add_Date)
      
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

##################################################################################################################################
        
#VI: Calculate the date of snowmelt for every 10mx10m pixel within the study area by fitting a GAM through the NDSI data
        
##################################################################################################################################            
         
  #(20): Calculate the date of snowmelt (NDSI <= NDSI_threshold) for every pixel within the study area 
     
     #Note that this is a very computationally demanding process due to the large number of pixels (>500,000) within the
     #study area in combination with a temporal range spanning +100 days. I have used parallel processing methods to speed
     #up the computational time. This has resulted in a reduction from c.a. 80 hours to c.a. 10 hours of processing time.
     
     #(A): Transform each image to a feature Collection of NDSI values for all pixels
     
       #Create an iteration function that we will use to iterate through all images of the image collection. For each image, the
       #value of the NDSI band is extracted for each pixel of the image. The resulting NDSI values (+lat/lon, datetime of the image)
       #of each pixel are stored as feature properties in a feature collection (i.e. where each pixel has a separate feature with
       #properties). This results in a feature collection of all pixels for the current image. This feature collection is then
       #appended to a list of feature collections from previous image iterations.

       #Store default Sentinel-2 image projection 
        S2Projection <- s2_clouds_filtered$first()$select("NDSI")$projection()
        #S2Projection$getInfo() 
     
       #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
        FC_initial <- ee$FeatureCollection(ee$List(list()))

       #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
       #(in this case the current iteration image) and the second element takes the output value from the iteration that preceeded it. The latter
       #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with
       #has been defined. Note that inside $map() functions all processing has to be done in the language of the server (javascript
       #Api of google earth engine). Thus, inside $map() functions, client-side functions such as $getInfo() or print(), paste0()
       #cannot be used.
        Extract_BandValuesAtPixels = Extract_BandValuesAtPixels #note that region=aoi_Shapefile

       #Iterate over the ImageCollection (output is a large feature collection)
        FC_merged <- ee$FeatureCollection(s2_clouds_filtered$select("NDSI")$iterate(Extract_BandValuesAtPixels, FC_initial))
        #FC_merged$first()$getInfo() #for debugging
        
       #Note that only pixels within the shapefile are included because we have already clipped all our images by this shapefile.  
       
     #(B): Transform feature collection to a dataframe:

       #!IMPORTANT! THIS STEP TAKES C.A. 8 - 16 HOURS
       #!IMPORTANT! MAKE SURE THERE IS SUFFICIENT SPACE ON YOUR GOOGLE DRIVE TO EXPORT THE TABLE (3.5 GB)
        
       #We use ee_table_to_drive() to prevent memory limits
        a=Sys.time()
        task_vector1 <- ee_table_to_drive(
          collection = FC_merged,
          description = paste0(timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_NDSI_Shapefile"),
          fileFormat = "CSV",
          selectors = c('NDSI', 'Date', 'lat', 'lon')
          )

        task_vector1$start()
        ee_monitoring(task_vector1, task_time=300, max_attempts = 1000000)
        #ee$data$getTaskList()
        #ee$data$cancelTask()
        
        exported_stats <- ee_drive_to_local(task = task_vector1, dsn=paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_NDSI_Shapefile"))
        df_pixel_ndsi <- read.csv(exported_stats)
        b=Sys.time()
        print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))

       # #Load dataframe (takes ca 2 minutes):
       #  df_pixel_ndsi <- read.csv(paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_NDSI_Shapefile.csv"))

       #Add day of year
        df_pixel_ndsi$doy <- as.numeric(format(as.POSIXct(df_pixel_ndsi$Date, format = "%Y-%m-%d %H:%M:%S"), "%j"))
        
       #Make sure each latitude/longitude combination gets its own pixel_ID (takes ca 1 minute):
        df_pixel_ndsi$pixel_ID <- paste0(format(df_pixel_ndsi$lat, nsmall = 5), "_", format(df_pixel_ndsi$lon, nsmall = 5))

     #(C): Calculate the date of snowmelt for each pixel in the dataframe
        
       #Loop through all pixel_ID's, select dataframe for that pixel containing NDSI values
       #as measured in all images in the image collection (52 in this case), fit gam through
       #the NDSI ~ doy data, determine date when NDSI<NDSI_threshold and store this date of snowmelt
       #together with the pixel_ID in a new dataframe.

       #Specify NDSI threshold (values larger than this threshold are considered snow covered):
        NDSI_threshold=NDSI_threshold
        
       #When running this code as a sequential process on a single core it takes c.a. 72 hours to run. To speed-up
       #this process we run the code in parallel on multiple local computer cores (4). We will add a progress-bar
       #to be able to track the progress. In addition, the progress is slowed down by having to work with this
       #gigantic dataset (df_pixel_ndsi) as one chunk. The code becomes MUCH faster if we split this up into
       #several sub-datasets (e.g. based on subsets of 10000 unique pixel_IDs).

       #Specify number of cores:
        numCores <- detectCores()

       #Make clusters and initialize using DoSNOW as this allows for a progress bar
        cl <- makePSOCKcluster(numCores)
        registerDoSNOW(cl)

       #Split dataset into chuncks of 10,000 pixels (resulting in 52 chunks)
        chunk_size = 10000
        pixelIDs <- unique(df_pixel_ndsi$pixel_ID)
        pixelIDs_split <- split(pixelIDs, ceiling(seq_along(pixelIDs)/chunk_size))

       #Specify function to sequentially remove outliers from the dataset for each pixel (results in better estimate of day of snowmelt).
        
        #We fit a Generalized Additive Model (GAM) through the NDSI data for every pixel in each data subset. We do this using a 
        #sequential outlier-filtering process. We first fit a GAM through the data and calculate model residuals. We then exclude 
        #all rows from the dataframe where the residual >= (0.25 * the range of the data). We then re-fit a GAM to this reduced 
        #dataset, and again calculate model residuals. We then exclude all rows from the reduced dataframe where the residual >= 
        #(0.1 * the range of the data). This gives us a final dataframe in which outliers are thus sequentially removed. This 
        #whole process is executed using the function f_gam_SeqRemOutliers. Note that a sequential step is required because 
        #initially some datapoints might falsely be assigned a large residual because of one extreme outlier. After removal of 
        #this extreme outlier and refitting of the GAM, it can be better assessed which data points truly have a large residual 
        #and can thus be assigned as 'true' outliers.
        f_gam_SeqRemOutliers <- f_gam_SeqRemOutliers #sourced
        
        #Specify sequential outlier thresholds (Note that these thresholds are relative to the range of the data. Thus if the response
        #variable ranges between -1 and +1, then the total range of the data is 2 and a threshold value of 0.4 thus corresponds to an
        #actual residual threshold of 0.4*2=0.8:
         outlier_thresh_1=outlier_thresh_1
         outlier_thresh_2=outlier_thresh_2
         outlier_removal=outlier_removal

       #Load a function that iterates through all data subsets. Within each data subset it calculates day of snowmelt 
       #for every pixel, by fitting a GAM with sequential outlier removal and then linearly approximating at which day of 
       #year the predicted NDSI value of this GAM changes from above outlier_threshold to below. The code employs 
       #parallel processing using foreach and %dopar% on four local computer cores.
        f_detect_threshold_date_parallel <- f_detect_threshold_date_parallel #sourced  
        
       #Run this function over all data subsets, combine the results and save the resulting dataframe
        system.time({
          print("Calculate the date of snowmelt for each pixel:")
          results <- lapply(1:length(pixelIDs_split), FUN=f_detect_threshold_date_parallel,
                            pixelIDs_split=pixelIDs_split, df_pixel_y=df_pixel_ndsi, pixel_ID_column="pixel_ID", 
                            y="NDSI", x="doy", pixel_gam_plots=F, y_threshold=NDSI_threshold)
          
          #Store date of snowmelt per pixel
          df_pixel_snowmelt <- lapply(results, "[[", 1)
          df_pixel_snowmelt <- as.data.frame(do.call(rbind, do.call(c, df_pixel_snowmelt)))
          colnames(df_pixel_snowmelt)[colnames(df_pixel_snowmelt)=="x_threshold"] <- "doy_snowmelt"
          write.csv(df_pixel_snowmelt, file=paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_Snowmelt_Shapefile.csv"), quote = FALSE, row.names=FALSE)
          
          # #Store plots per pixel
          # plot_pixel_snowmelt <- lapply(results, "[[", 2)
          # plots_per_page = 25
          # plot_pixel_snowmelt <- lapply(plot_pixel_snowmelt, function(x){split(x, ceiling(seq_along(plot_pixel_snowmelt[[1]])/plots_per_page))})
          # plot_pixel_snowmelt <- unname(unlist(plot_pixel_snowmelt, recursive = F))
          # pdf(paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Plot_Pixel_Snowmelt_Shapefile.pdf"), width=20, height=16, onefile = TRUE)
          # for (i in seq(length(plot_pixel_snowmelt))) { do.call("grid.arrange", plot_pixel_snowmelt[[i]]) }
          # dev.off()
          })

       #this took 17685 seconds (5 hours). Parallel processing thus decreased computation time from 80+ hours
       #to c.a. 5 hours!

       # #Read dataframe
       #  df_pixel_snowmelt <- read.csv(file=paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_Snowmelt_Shapefile.csv"), header=TRUE)
       #  df_pixel_snowmelt now contains the date of snowmelt for each pixel within the area depicted by aoi_Shapefile.  
        
       #Clean up the cluster after finishing the parallel runs
        stopCluster(cl)

       #Turn parallel processing off and run sequentially again after this point
        registerDoSEQ()
     
       #The dataframe df_pixel_snowmelt now contains the date of snowmelt for each individual pixel_ID (514584 10mx10m pixels).
       #To be able to plot these data we transform this dataframe to a feature collection and then transform this feature
       #collection to an image.
        
   #(21): Transform df_pixel_snowmelt to a feature collection with random geometry
        
       #Generate some random longitude and latitude values as this is required for an sf object:
        df_pixel_snowmelt$lon <- runif(nrow(df_pixel_snowmelt), 0, 74)
        df_pixel_snowmelt$lat <- runif(nrow(df_pixel_snowmelt), -20, 20)
        
       #Prevent NA in doy_snowmelt column (earthengine cannot deal with NA)  
        df_pixel_snowmelt$doy_snowmelt[is.na(df_pixel_snowmelt$doy_snowmelt)] <- -9999
        
       #memory limits prevent us from processing the whole dataframe at once. Therefore we split it up into subsets of 10000 rows
        chunk_size = 10000
        rowIDs <- 1:nrow(df_pixel_snowmelt)
        rowIDs_split <- split(rowIDs, ceiling(seq_along(rowIDs)/chunk_size))
      
       #Iterate through all dataframe subsets, convert each to a feature collection, and append each to the feature collection FC_initial.
       #Note that this step takes approximately XX hours!
        FC_initial <- ee$FeatureCollection(ee$List(list()))
        for(i in 1:length(rowIDs_split)){ 
          
          #Select subset of rows
           rowID_min <- min(rowIDs_split[[i]])
           rowID_max <- max(rowIDs_split[[i]])
          
          #Select subset of dataframe:
           df_tmp <- df_pixel_snowmelt[rowID_min:rowID_max,]
          
          #Change subset-dataframe to a SF object 
           df_sf_tmp <- st_as_sf(x = df_tmp,                         
                                 coords = c("lon", "lat"),
                                 crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
          
          #Change sf object to an earth engine feature collection by uploading it to the asset folder
           FC_tmp <- sf_as_ee(
             x = df_sf_tmp,
             assetId = paste0(path_asset, "/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_df_pixel_snowmelt_", i),
             overwrite = TRUE,
             monitoring = TRUE,
             via = 'getInfo_to_asset')
           #FC_tmp <- sf_as_ee(df_sf_tmp)
          
          #Add feature collection to an expanding feature collection:
           FC_initial <- FC_initial$merge(FC_tmp)
           FC_initial <- ee$FeatureCollection(FC_initial)               
          
          #print progress
           print(paste0("Progress: ", (100*i)/(length(rowIDs_split))))
           #print(paste0("Size of FC_initial: ", FC_initial$size()$getInfo()))
          
        }
        FC_pixels_snowmelt <- ee$FeatureCollection(FC_initial)
        # FC_pixels_snowmelt$first()$getInfo()
        # FC_pixels_snowmelt$size()$getInfo()
        
        #Inspect assets folder (for debugging)
        #  ee_manage_quota()
        #  ee_manage_assetlist(path_asset)
        
    #(22): Add geometry (latitude and longitude) of each pixel_ID to FC_pixels_snowmelt
         
       #The feature collection that we want to construct should contain 'pixel_ID' and 'doy_snowmelt' as properties and should contain
       #the original geometry corresponding to each pixel_ID. So far, FC_pixels_snowmelt contains a separate feature for each pixel, where
       #each feature contains doy_snowmelt and pixel_ID as a property. However, the geometry (lan/lon) of each feature is randomly chosen. 
       #We need to make sure that the actual geometry matching each pixel_ID is added instead of this random geometry.
        
       #To obtain the corresponding geometry (lat/lon) of every pixel_ID within aoi_Shapefile, we sample from a single image on a 'resolution' 
       #resolution using img$sample(). This gives a single distinct feature for each pixel, including their geometry. We can then re-construct
       #the property pixel_ID for every feature in this feature collection. The resulting feature collection is called FC_pixels_distinct. 
       #The final step is then to join FC_pixels_snowmelt to FC_pixels_distinct based on an inner join with pixel_ID.
       
       #This step results in memory errors (on the server side) when the image contains more than 2.5 million pixels. In the latter case
       #we need to split up our study area in several sub-areas (i.e. create separate adjacent shapefiles using QGIS). If such errors are
       #encountered here then script "6-RGEE_TomVersluijs_PixelDateOfSnowmelt_TAY" should be run instead. 
    
        #(A): Select a single image from the image collection
         img <- s2_clouds_filtered$first()$select("NDSI")
        
        #(B): Create a feature collection of points at the center of each pixel within aoi_Shapefile (each point thus corresponds to a pixel in the image)
         FC_pixels_distinct <- img$sample(
          region=aoi_Shapefile, #Sample all pixels within aoi_Shapefile
          geometries=TRUE,  #Store lat/lon in geometry property
          projection=S2Projection, #Set to native projection of S2 image
          scale=resolution, #10 meter native resolution of image (implicitly set using projection above)
          seed=23, #set fixed seed to be able to reproduce the sample (same as for FC_merged above)
          dropNulls=FALSE) #Make sure NULL pixels are not dropped from image collection
      
        #(C): Make sure there is a NDSI value at each feature within the feature collection:
        #Set the band value to a no data value of -9999 for all features where the band value is NULL.
         FC_pixels_distinct <- FC_pixels_distinct$map(function(feature){
          ndsi <- ee$List(list(feature$get('NDSI'), -9999))$reduce(ee$Reducer$firstNonNull())
          return(feature$set("NDSI", ndsi))})  
        
        #(D): Add latitude and longitude of each pixel as a property to each feature
         FC_pixels_distinct <- FC_pixels_distinct$map(function(feature){
          coordinates <- feature$geometry()$coordinates()
          lon <- coordinates$get(0)
          lat <- coordinates$get(1)
          return(feature$set('lon', lon)$set('lat', lat))
          })
        
        #(E): Add a pixel_ID property to each feature
         FC_pixels_distinct <- FC_pixels_distinct$map(function(feature){
          lon_tmp <- feature$get("lon")
          lat_tmp <- feature$get("lat")
          lon_tmp <- ee$Number$format(lon_tmp, '%.5f')
          lat_tmp <- ee$Number$format(lat_tmp, '%.5f')
          string_tmp <- ee$String("_")
          pixel_ID <- lat_tmp$cat(string_tmp)$cat(lon_tmp)
          return(feature$set('pixel_ID', pixel_ID))
          })  
         
        #(F): Reduce feature collection properties to only contain pixel_ID  
         FC_pixels_distinct <- FC_pixels_distinct$map(function(feature) {
          properties_all <- feature$propertyNames()
          properties_selection <- properties_all$filter(ee$Filter$inList('item', list('pixel_ID')))
          return(feature$select(properties_selection))
          })
        
         # FC_pixels_distinct$first()$getInfo()
         # FC_pixels_distinct$size()$getInfo()
         
         #Thus, FC_pixels_distinct is a featurecollection that contains a single distinct feature for every pixel_ID,
         #and contains the original lat/lon geometry of every pixel. We will now add doy_snowmelt as a new property
         #to FC_pixels_distinct by merging it with FC_pixels_snowmelt.
         
        #(G): Add doy_snowmelt as property to FC_pixels_distinct   
        
          #Use an equals filter to specify how the collections match.
           Filter_pixel_ID <- ee$Filter$equals(
              leftField= 'pixel_ID',
              rightField= 'pixel_ID')
         
          #Define the join.
           innerJoin <- ee$Join$inner()
         
          #Apply the join.
           Pixel_ID_Join <- innerJoin$apply(FC_pixels_distinct, FC_pixels_snowmelt,  Filter_pixel_ID)
         
          #Add features of second ('secondary') feature collection to those of the first ('primary') feature collection 
           FC_Combined <- Pixel_ID_Join$map(function(pair) {
             f1 <- ee$Feature(pair$get('primary'))
             f2 <- ee$Feature(pair$get('secondary'))
             return(f1$set(f2$toDictionary()))
             })  
          
        #(H): To speed-up further computations with this feature collection we as an intermediate step upload and re-download it
         
          #Delete FC_pixels_snowmelt_optimized if it already occured in the asset folder:
          #tryCatch(ee_manage_delete(paste0(path_asset, "/FC_pixels_snowmelt_optimized")), error=function(error_message) {message("path_asset does not yet exist")})

          #Upload to asset folder (This step takes 35 minutes):
           assetid2 <- paste0(path_asset, "/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_FC_pixels_snowmelt_optimized")
             task_vector2 <- ee_table_to_asset(
               collection = FC_Combined,
               overwrite = TRUE,
               assetId = assetid2
               )
             task_vector2$start()
             ee_monitoring(task_vector2, task_time=300, max_attempts = 1000000)
         
          #Check assets folder:
           #ee_manage_quota()
           #ee_manage_assetlist(path_asset)
          
          #Save assetid2 for future downloading of FC_pixels_snowmelt_optimized
           saveRDS(object=assetid2, file=paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Variable_AssetID.Rds"))
           
          #Get feature collection from asset folder and create FC_pixels_snowmelt_optimized
           #assetid2=paste0(path_asset, "/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_FC_pixels_snowmelt_optimized")
           FC_pixels_snowmelt_optimized <- ee$FeatureCollection(assetid2) 
           #FC_pixels_snowmelt_optimized$first()$getInfo()
           #FC_pixels_snowmelt_optimized$size()$getInfo()
       
    #(23): Transform the Feature collection FC_pixels_snowmelt_optimized to an image (with doy_snowmelt as an image band)
           
          #(A): Reduce feature collection FC_pixels_snowmelt_optimized to an Image with a 'resolution' resolution:    
           image_snowmelt <- FC_pixels_snowmelt_optimized$
             filter(ee$Filter$notNull(list('doy_snowmelt')))$ #pre-filter data for nulls that cannot be turned into an image
             filter(ee$Filter$neq(name='doy_snowmelt', value=-9999))$ #pre-filter data for -9999 values
             reduceToImage(properties=list('doy_snowmelt'), reducer=ee$Reducer$first()$setOutputs(list('doy_snowmelt')))$
             reproject(crs=crs, crsTransform=NULL, scale=resolution)
         
          #(B): Extract day of snowmelt in year of interest for a single point
           ee_extract(x=image_snowmelt, y=coordinates_point, fun=ee$Reducer$first(), scale=resolution, sf=TRUE)
         
          #(C): Plot snowmelt day of year as a coloured image
           Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
           Map$addLayer(s2_col$filterDate(paste0(year_ID, "-06-30"), end_date)$first(),list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
           Map$addLayer(image_snowmelt,list(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')
     
          #(D): Export original image to Google Drive (takes c.a. 2 minutes):
         
            #Create task to export the original doy_snowmelt image to Google Drive  
             task_vector3 <- ee_image_to_drive(
              image=image_snowmelt,
              description= paste0(timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt'),
              #scale= resolution,
              region=aoi,
              crs=crs,
              maxPixels=1e13,
              dimensions=ee$Number(4096),
              fileFormat='GeoTIFF'
              )
           
            #Start and monitor export task:
             task_vector3$start()
             ee_monitoring(task_vector3, task_time=30, max_attempts = 1000000)
             ee_drive_to_local(task = task_vector3, dsn=paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt'))
             
          #(E): Export RGB image to Google Drive (takes c.a. 2 minutes):
         
            #Convert original image to an RGB image:
             image_snowmelt_RGB <- image_snowmelt$visualize(bands=c('doy_snowmelt'), min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red'))
             #ee_print(image_snowmelt_RGB)
           
            #Create task to export RGB image to Google Drive:  
             task_vector4 <- ee_image_to_drive(
               image=image_snowmelt_RGB,
               description= paste0(timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt_RGB'),
               #scale= resolution,
               region=aoi,
               crs=crs,
               maxPixels=1e13,
               dimensions=ee$Number(4096),
               fileFormat='GeoTIFF'
               )
          
            #Start and monitor export task:
             task_vector4$start()
             ee_monitoring(task_vector4, task_time=30, max_attempts = 1000000)
             ee_drive_to_local(task = task_vector4, dsn=paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_Image_DoySnowmelt_RGB"))  
             
    #(24): Save workspace 
     save.image(paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Backup_Workspace_PixelDoySnowmelt.RData"))


###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################

          