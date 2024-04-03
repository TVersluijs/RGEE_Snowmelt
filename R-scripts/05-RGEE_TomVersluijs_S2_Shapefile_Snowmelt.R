##################################################################################################################################

#Use Sentinel-2 data to extract time series of the average NDVI, NDMI, and NDSI and the fraction of snowcover for a single area of 
#interest (i.e. a single shapefile). The fraction of snowcover is estimated by calculating the fraction of pixels with an NDSI value
#larger than the user specified NDSI-threshold for each timestep. This corresponds to the method='snowfraction' in other scripts. 
#The user can specify  whether clouds and permanent waterbodies need to be masked, and whether a composite image per day of year 
#needs to be generated (merging multiple satellite photos for that day). The current script requires a shapefile of the study area 
#as input. It only works for small areas of c.a. 50 km2 (larger areas might result in computation errors unless the spatial 
#resolution of the analysis is decreased).

#Note that in this script (05), snowmelt is not calculated based on pixel-level snowmelt data (i.e. 'pixel_gam' method is not 
#implemented). This approach is instead implemented in script "08-RGEE_TomVersluijs_S2_Shapefile_Pixel_Snowmelt.R" and involves 
#fitting of GAMS through NDSI data per pixel and extracting the moment this GAM passes a user-defined NDSI-threshold. This results 
#in a pixel-level map of the date of snowmelt. Script "10-RGEE_TomVersluijs_S2_ExtractSnowFraction.R" can then be used to extract 
#timeseries of the fraction of snowcover for points/polygon(s) of interest from this map.

#Copyright Tom Versluijs 2024-04-03. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R". 
#Note that a GoogleDrive is required. Important: make sure to run this script from within the "RGEE_Snowmelt.Rproj" project file.

##################################################################################################################################

#I: Setup workspace

##################################################################################################################################

  #(1): Clear workspace and set python environment
   rm(list=ls())
   utils::install.packages("here")
   library(here)
   if(file.exists(paste0(here::here(), "/Input/rgee_environment_dir.rds"))){
     rgee_environment_dir <- readRDS(paste0(here::here(), "/Input/rgee_environment_dir.rds"))
     reticulate::use_python(rgee_environment_dir, required=T)
     reticulate::py_config()
     }
    
  #(2): Load packages
   #renv::restore() #revert to last version of R-packages used to successfully run this script (optional).
   utils::install.packages("pacman")
   library(pacman)
   p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr, tidyr, foreach, parallel, doSNOW)

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
     resolution=20 # default maximum resolution for Sentinel-2 = 10m
     
   #(b) Area of interest

     #Specify name of study area (used as prefix in output files)
     area_name="ZAC"

     #Name of Shapefile (located in '/Input/Shapefiles' folder)
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
     start_date <- paste0(year_ID, "-03-15") #-03-15
     end_date <- paste0(year_ID, "-09-15") #-09-15

   #(d) Snow detection

     #NDSI threshold for snow detection
     NDSI_threshold_vector = c(0.4, 0.3, 0.5)

     #Define the snowcover fraction within the aoi for which the date of its occurrence will be calculated
     Snowfraction_threshold_vector = c(0.25, 0.5, 0.75) 
    
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
     #First, all cloudy pixels, without cloudy pixels within radius cld_erosion, will be removed.
     #Second, all remaining pixels will be dilated by radius cld_buffer. Make sure that cld_buffer
     #is larger than cld_erosion. Larger values result in a more coarse cloud-filtering process.
     cld_erosion=250 #default = 250 for conservative cloud estimation
     cld_buffer=270 #default = 270 for conservative cloud estimation

     #Define resolution of cloud-filtering (in meters)
     resolution_cldmsk=150 #(default = 150)

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
     outlier_thresh_2=0.2 #second threshold in outlier removal (relative to y-range of data, default=0.2)
       
     #Specify the degree of GAM smoothing by setting the 'k' parameter
     gam_k_outlier=10 #Number of knots when filtering outliers (default=10)
     gam_k=25 #Number of knots when making model predictions (default=25). 
     #Larger values result in a more precise GAM-fit, but might result in overfitting.
 
     
##################################################################################################################################
 
 #III: Define some additional parameters (automated)
 
##################################################################################################################################
 
  #(5): Automatically define some additional parameters 
 
    #Create a unique data_ID
     data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_S2")           
 
    #Create a timestamp variable
     timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
     
    #Specify starting and ending date as day of year 
     start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
     end_date_doy <- as.numeric(strftime(end_date, format = "%j"))   
     
    #Set mask_clouds to TRUE if a composite image needs to be generated         
     if(create_composite==TRUE){mask_clouds=TRUE}      
 
    #Set mask_clouds to TRUE if water masking is TRUE and mask_water_type==water_mask_Manual
     if(mask_water==TRUE & (mask_water_type=="water_mask_Manual" | mask_water_type=="both")){mask_clouds=TRUE}      
     
    #Create output folder
     if(dir.exists(paste0(here(), "/Output/S2/05_Shapefile_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/S2/05_Shapefile_Snowmelt"), recursive = TRUE)}
     
    #Save all parameters and their values in the environment to a text file 
     file_conn <- file(paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Parameters.txt"), "w")
     for (obj in setdiff(ls(), lsf.str())) {cat(paste(obj, "=", get(obj)), file = file_conn) ; cat("\n", file = file_conn)}
     close(file_conn)
     
##################################################################################################################################
       
#IV: Read and display the unfiltered data
       
##################################################################################################################################
    
  #(6): Read study area shapefile and convert to a feature collection.
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
      paste0('Size of study area calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')
    
  #(7): Extract satellite data including cloud probability band:

    #Extract Sentinel-2 Surface Reflectance data for a given area and date range:
    s2_col<-ee$ImageCollection(paste0('COPERNICUS/', s2_dataset))
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
    
    #Clip all images to the area depicted by 'aoi_Shapefile'
    s2_col <- s2_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})

  #(8): Plot an RGB, NDSI and NDVI image for a single extracted satellite image (for debugging)
    
    # #Select a single image for initial plot:
    # image <- s2_col$filterDate("2022-07-15", end_date)$first()
    # 
    # #Normalized difference snow index:
    # ndsi_s2 <-  image$expression('(B3 - B11) / (B3 + B11)',
    #                              list('B11'=image$select('B11'),
    #                                   'B3'=image$select('B3')))
    # 
    # #Normalized difference vegetation index:
    # ndvi_s2 <- image$expression('(B8 - B4) / (B8 + B4)',
    #                             list('B8'=image$select('B8'),
    #                                  'B4'=image$select('B4')))
    # 
    # #Add Normalized difference water index to this image:
    # ndwi_s2 <- image$expression('(B3 - B8) / (B3 + B8)',
    #                             list('B3'=image$select('B3'),
    #                                  'B8'=image$select('B8')))
    # 
    # #Plot all layers
    # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 6)
    # Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
    # Map$addLayer(ndsi_s2,list(min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
    # Map$addLayer(ndvi_s2,list(min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
    # Map$addLayer(ndwi_s2,list(min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')
    
  #(9): Add NDSI and NDVI to the clipped image collection
    
    #Map Normalized difference functions over the image collection
    s2_col <- s2_col$
      map(getNDSI)$
      map(getNDVI)$
      map(getNDMI)$
      map(getNDWI)
  
  #(10): Create a timeseries gif of RGB images for the shapefile area (for debugging)
    
    # #Check number of images in collection
    # s2_col$size()$getInfo()
    # 
    # #Create a timelapse video of RGB band
    # videoArgs <- list(dimensions=400, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
    # tryCatch({browseURL(s2_col$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})

##################################################################################################################################
        
#V: Filter and mask clouds within the image collection
        
##################################################################################################################################
        
   #There are many images that are almost fully covered in dense clouds. We thus first filter and mask clouds in the image collection.
   #We delete all images with a cloud cover fraction > max_cloud_fraction. In the remaining images we mask all pixels that are 
   #covered by clouds.
    if(mask_clouds==TRUE){
    
      #print message   
      print("Cloud masking = TRUE")
      
      #(9): Add cloud information within the aoi_Shapefile to the image collection: 
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
        
        # #Visual check of cloud mask (for debugging)
        # image <- s2_col$filterDate(paste0(year_ID, "-07-15"), end_date)$first()$
        #   clipToCollection(aoi_Shapefile)$
        #   select("clouds_unadjusted", "clouds", "clouds_inv", "clouds_prob", "B4", "B3", "B2")
        # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 6)
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
        #between both measures can for example arise when the part of the aoi_Shapefile was covered by clouds while most of the larger
        #region (tile) was not or vice versa.

        # #Extract cloud fraction of all images in image collection for the shapefile area
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
        #   theme_classic()
        #   #Thus, this indicates that our local cloud measure does not necessarily correspond to that as measured on the whole tile as
        #   #explained above (which is expected as my measure is more representative of the cloudcover within the study area).
       
      #(10): Exclude all images from the image collection that have a CloudFraction value > max_cloud_fraction, and mask all cloud pixels within the remaining images
    
        #Filter and mask the image collection for clouds
         s2_clouds_filtered <- s2_col$
          #Pre-filter to get less cloudy granules.
          filter(ee$Filter$lt('CloudFraction', max_cloud_fraction))$
          #Apply cloudmask for individual pixels
          map(Add_CloudMask)
        
        # #Create timelapse video of the cloud filtered/masked RGB images (for debugging)
        #  videoArgs <- list(dimensions=180, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
        #  tryCatch({browseURL(s2_clouds_filtered$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
         
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

#VI: Mask permanent waterbodies (ponds, lakes, rivers, sea) within the image collection

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
                                           start_date_NDSI=start_date_NDSI, end_date_NDSI=end_date_NDSI, NDSI_threshold=NDSI_threshold_vector[1], NIR_threshold=NIR_threshold)
           
           
       }
         
      #(D): Apply the water masking function to each image in the collection:
         
       #Define the final water_mask:
       if(mask_water_type=="water_mask_ESA"){water_mask <- water_mask1}
       if(mask_water_type=="water_mask_Manual"){water_mask <- water_mask2}
       if(mask_water_type=="both"){water_mask <- water_mask1$Or(water_mask2)}     
         
       # #Plot permanent waterbodies(for debugging)
       # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 6)
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
    
#VII: Create a composite image for each day by mosaicking all images from that day
    
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
    # doy=176
    # s2_img <- s2_clouds_filtered$filterMetadata('year', 'equals', as.numeric(year_ID))$filterMetadata('doy', 'equals', doy)
    # s2_img <- s2_img$select(c('B4', 'B3', 'B2', 'NDSI', 'clouds_inv'))
    # 
    # # #Number of images in s2_clouds_filtered
    # # s2_img$size()$getInfo()
    # # s2_img$first()$bandNames()$getInfo()
    # # s2_img$first()$propertyNames()$getInfo()
    # 
    # # #Inspect whether there is any overlap between the tiles for this day:
    # # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 6)
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
    #   snow <- ndsi$gt(NDSI_threshold_vector[1])$rename('SNOW')
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
    #   #Print snowfraction
    #   SnowFraction$getInfo()
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

#VIII: Count the total number of unmasked pixels per doy within aoi_Shapefile

##################################################################################################################################

  #(18): Count the total number of unmasked pixels and the total number of pixels per doy within aoi_Shapefile

     #(A): Add pixel counts within the area of interest to each separate image by mapping the pixel count functions over the image collection
      s2_col_composite <- s2_col_composite$
        #Count number of unmasked pixels within aoi_Shapefile
        map(AddPixelCount)$
        #Add NULL to those images in which pixel count could not be calculated
        map(AddNULLPixelCount)

     #(B): Extract pixel_counts of all images in image collection for aoi_Shapefile

       #create a current timestamp to prevent identical names on Google Drive
        current_timestamp0 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))

       #We use ee_table_to_drive() to prevent memory limits
        task_vector0 <- ee_table_to_drive(
          collection = s2_col_composite,
          description = paste0(current_timestamp0, "_", data_ID, "_Res", resolution, "_Pixel_Counts_polygon"),
          fileFormat = "CSV",
          selectors = c('doy', 'unmasked', 'total')
          )

       #Monitor the task
        task_vector0$start()
        print("Count the number of unmasked pixels within the shapefile per doy:")
        ee_monitoring(task_vector0, task_time=60, max_attempts=1000000)

       #Export results to local folder
        exported_stats <- ee_drive_to_local(task = task_vector0, dsn=paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Pixel_Counts_polygon"))
        df_pixelcount <- read.csv(exported_stats)
        unlink(exported_stats)

     #(C): Replace -9999 values by NA (for debugging)
      df_pixelcount$unmasked[df_pixelcount$unmasked < -9000] <- NA
      df_pixelcount$total[df_pixelcount$total < -9000] <- NA
      df_pixelcount$doy[df_pixelcount$doy < -9000] <- NA

     #(D): Add missing dates with 0 unmasked pixels to the dataframe
      df_doy_missing <- data.frame(doy=seq(start_date_doy, end_date_doy)[!(seq(start_date_doy, end_date_doy) %in% df_pixelcount$doy)],
                                   unmasked=0,
                                   total=max(df_pixelcount$total))
      df_pixelcount <- rbind(df_pixelcount, df_doy_missing)
      df_pixelcount <- df_pixelcount[order(df_pixelcount$doy),]
      write.csv(df_pixelcount, file=paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Pixel_Counts_polygon.csv"), quote=FALSE, row.names=FALSE)

     #(E): Create barplot with the pixel counts per day of year within aoi_Shapefile
      df_pixelcount$masked <- df_pixelcount$total - df_pixelcount$unmasked
      df_pixelcount <- df_pixelcount[,-which(colnames(df_pixelcount) %in% "total")]
      df_pixelcount <- tidyr::pivot_longer(df_pixelcount, cols=c(unmasked, masked), names_to = "pixels")
      p_pixelcounts <- ggplot()+
       geom_bar(data=df_pixelcount, aes(x=doy, y=value, fill=pixels), stat="identity",colour="black", position="stack", width=1)+
       #geom_rect(aes(xmin=min(df_pixelcount$doy)-1, xmax=max(df_pixelcount$doy)+1, ymin=0, ymax=max(df_pixelcount$value)), alpha=0.5, fill="white")+
       scale_fill_manual(values = c("black", "#FFA52C"))+
       xlab("Time (day of year)")+
       ylab("Pixel count")+
       theme_classic()

     #(F): Save barplot
      pdf(paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Plot_Pixel_Counts_polygon.pdf"), width=12, height=8)
      print(p_pixelcounts)
      dev.off()

##################################################################################################################################
         
#IX: Create timelapse videos of the RGB bands, and the NDSI-, NDVI-, NDMI- and NDMI-bands (Cloud and Water masked)
         
##################################################################################################################################
    
   #Note that these GIFs generally result in memory errors with areas larger than 100km2.  
         
     # #(19): Create a timeseries GIF of RGB images
     #     videoArgs <- list(dimensions=100, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
     #     tryCatch({browseURL(s2_col_composite$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
     #     #Note that missing pixels are actually not recorded by the satellite and are NOT caused by coding errors
     # 
     # #(20): Create a timeseries GIF of NDSI images
     # 
     #    #Map a visualization function over the image collection. This function converts each image to an RGB image (i.e. a red, green and blue
     #    #colour band with min and max values of 0 and 255 respectively). In the function it can be specified which data band needs to be converted,
     #    #what min and max values of the original band should correspond to 0 and 255 in the RGB image, and which colour pallete should be used.
     #    palette=c('black', '0dffff', '0524ff', 'ffffff')
     #    visFun_NDSI <-  function(img) {
     #      return(img$visualize(bands='NDSI', min=-1, max=1.5, palette=palette)$
     #               copyProperties(img, img$propertyNames()))}
     #    S2_snow_masked_RGB <- s2_col_composite$map(visFun_NDSI)
     # 
     #    #Create a timelapse video of NDSI band
     #    videoArgs <- list(dimensions=100, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
     #    tryCatch({browseURL(S2_snow_masked_RGB$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
     # 
     #  #(21): Create a timeseries GIF of NDVI images
     #    palette <- c("#cccccc", "#f46d43", "#fdae61", "#fee08b", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850")
     #    visFun_NDVI <-  function(img) {
     #      return(img$visualize(bands='NDVI', min=-0.25, max=1, palette=palette)$
     #               copyProperties(img, img$propertyNames()))}
     #    S2_ndvi_masked_RGB <- s2_col_composite$map(visFun_NDVI)
     # 
     #    #Create a timelapse video of NDSI band
     #    videoArgs <- list(dimensions=510, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
     #    tryCatch({browseURL(S2_ndvi_masked_RGB$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
     # 
     #  #(22): Create a timeseries GIF of NDMI images
     #    palette <- c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#6ad99e", "#387ad9", "#003dd6")
     #    visFun_NDMI <-  function(img) {
     #      return(img$visualize(bands='NDMI', min=-0.75, max=1, palette=palette)$copyProperties(img, img$propertyNames()))}
     #    S2_ndmi_masked_RGB <- s2_col_composite$map(visFun_NDMI)
     # 
     #    #Create a timelapse video of NDMI band
     #    videoArgs <- list(dimensions=510, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
     #    tryCatch({browseURL(S2_ndmi_masked_RGB$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
     # 
     #  #(23): Create a timeseries GIF of NDWI images
     #    palette <- c('000000', '0dffff', '0524ff', 'ffffff')
     #    visFun_NDWI <-  function(img) {
     #      return(img$visualize(bands='NDWI', min=-0.5, max=1, palette=palette)$copyProperties(img, img$propertyNames()))}
     #    S2_ndwi_masked_RGB <- s2_col_composite$map(visFun_NDWI)
     # 
     #    #Create a timelapse video of NDWI band
     #    videoArgs <- list(dimensions=510, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
     #    tryCatch({browseURL(S2_ndwi_masked_RGB$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
        
##################################################################################################################################
        
#X: Calculate and visualize the change in the fraction of snowcover over time
        
##################################################################################################################################       
        
  #(24): Extract average NDVI, NDSI and NDWI for the area of interest
  
      #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
      FC_initial_BandValues <- ee$FeatureCollection(ee$List(list()))
      
      #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
      #(in this case the current iteration image) and the second element takes the output value from the iteration that preceeded it. The latter
      #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with 
      #has been defined.
      Extract_BandValuesAtPolygon = Extract_BandValuesAtPolygon #sourced
      
      #Iterate over the ImageCollection 
      FC_merged_BandValues <- ee$FeatureCollection(s2_col_composite$select("NDSI", "NDVI", "NDMI")$iterate(Extract_BandValuesAtPolygon, FC_initial_BandValues))
      
      #create a current timestamp to prevent identical names on Google Drive
      current_timestamp1 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
      
      #Setup task
      task_vector1 <- ee_table_to_drive(
        collection= FC_merged_BandValues,
        description = paste0(current_timestamp1, "_", data_ID, "_Res", resolution, "_Data_MeanBandValues"),
        folder="RGEE_tmp",
        fileFormat="CSV",
        selectors=c("NDSI", "NDVI", "NDMI", 'Date')
        )
      
      #Run and monitor task
      print(paste0("calculating Mean BandValues for ", data_ID))
      task_vector1$start() 
      ee_monitoring(task_vector1, task_time=60, max_attempts=1000000)
      
      #Import results
      exported_stats <- ee_drive_to_local(task=task_vector1, dsn=paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Data_MeanBandValues"))
      df_BandValues <- read.csv(exported_stats)
      unlink(exported_stats)
      
      #Add day of year
      df_BandValues$Date <- as.POSIXct(df_BandValues$Date, format="%Y-%m-%d %H:%M:%S")
      df_BandValues$doy <- as.numeric(strftime(df_BandValues$Date, format = "%j"))
      
      #Remove NAs in the NDVI, NDMI and NDSI variables
      df_BandValues$NDVI[df_BandValues$NDVI < -9000] <- NA #replace -9999 by NA
      df_BandValues$NDSI[df_BandValues$NDSI < -9000] <- NA #replace -9999 by NA
      df_BandValues$NDMI[df_BandValues$NDMI < -9000] <- NA #replace -9999 by NA

      #Sort dataframe by doy
      index <- with(df_BandValues, order(doy))
      df_BandValues <- df_BandValues[index,]
      df_BandValues <- df_BandValues[,c("Date", "doy", "NDVI", "NDSI", "NDMI")]
      
      #Save dataframe with Mean BandValues
      write.csv(df_BandValues, paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Data_MeanBandValues.csv"), row.names = FALSE)
    
  #(25): Extract fraction of snowcover within aoi_Shapefile for each image in the image collection  
  
      #Create an empty dataframe for storing the Snowfraction results per NDSI-threshold
      df_SnowCover <- data.frame(Date=character(),
                                 doy=numeric(),
                                 SnowFraction=numeric(),
                                 NDSI_threshold=character())
      
      #Run the analysis for each level of NDSI_threshold_vector
      for(NDSI_threshold in NDSI_threshold_vector){  
        
        #Print message
        print(paste0("  -Extracting Snowfraction timeseries data for NDSI_threshold = ", NDSI_threshold))
        
        #Store NDSI_threshold as a character (used for naming of outputs)
        NDSI_threshold_char <- gsub("\\.", "_", as.character(NDSI_threshold))
        
        #Create a binary NDSI image and calculate the fraction of pixels above NDSI threshold over time
        
          #Map Snow computation functions over the mosaicked image collection
          s2_snow_filtered <- s2_col_composite$
            #Determine which pixels are snow-covered (NDSI > NDSI threshold)
            map(computeSnow)$
            #add the fraction of snow covered pixels as an image property (excluding cloud masked pixels from the calculations)
            map(AddSnowFraction)$
            #Add a NULL value to images for which the snow fraction could not be calculated
            map(AddNULLSNOW)
          
          # #Mask and display the binary layer (for debugging).
          #  img_snow <- s2_snow_filtered$filterDate(paste0(year_ID, "-05-26"), end_date)$first()
          #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
          #  Map$addLayer(img_snow,list(bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
          #  Map$addLayer(img_snow$select('SNOW'))+
          #  Map$addLayer(img_snow$select('NDSI'), list(min=-1, max = 1, palette = c('lightblue', 'darkblue')))
            
        #Extract fraction of snowcover within aoi_Shapefile for each image in the image collection   
        
          #Create an empty feature collection:
          FC_initial <- ee$FeatureCollection(ee$Feature(NULL))
          
          #Define an Iteration function to extract the fraction of snowcover within the area of interest for each image
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
          FC_merged <- ee$FeatureCollection(s2_snow_filtered$iterate(extract_snowcover, FC_initial))
  
          #Export the feature collection as a .csv table 
          #We export the data instead of using aggregate_array() as the latter might fail due to computation timeouts.
          
          #create a current timestamp to prevent identical names on Google Drive
          current_timestamp2 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
          
          #Setup task
          task_vector2 <- ee_table_to_drive(
            collection= FC_merged,
            description = paste0(current_timestamp2, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_SnowFraction"),
            folder="RGEE_tmp",
            fileFormat="CSV",
            selectors=c('SnowFraction', 'Date')
            )
          
          #Run and monitor task
          task_vector2$start() 
          ee_monitoring(task_vector2, task_time=60, max_attempts=1000000) #250s at 100m
          
          #Import results
          exported_stats <- ee_drive_to_local(task=task_vector2, dsn=paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_SnowFraction"))
          df_SnowCover_new <- read.csv(exported_stats)
          unlink(exported_stats)
          
          #Add day of year
          df_SnowCover_new$Date <- as.POSIXct(df_SnowCover_new$Date, format="%Y-%m-%d %H:%M:%S")
          df_SnowCover_new$doy <- as.numeric(strftime(df_SnowCover_new$Date, format = "%j"))
          
          #Remove NAs in the SnowFraction variable  
          df_SnowCover_new$SnowFraction[df_SnowCover_new$SnowFraction < -9000] <- NA #replace -9999 by NA
          df_SnowCover_new <- df_SnowCover_new[!is.na(df_SnowCover_new$SnowFraction), ]
        
          #Sort dataframe by doy
          index <- with(df_SnowCover_new, order(doy))
          df_SnowCover_new <- df_SnowCover_new[index,]
          df_SnowCover_new <- df_SnowCover_new[,c("Date", "doy", "SnowFraction")]
          
          #Add NDSI_threshold as a new column
          df_SnowCover_new$NDSI_threshold <- as.factor(NDSI_threshold) 
          
          #Add dataframe for current NDSI_threshold to dataframe from previous iterations:
          df_SnowCover <- rbind(df_SnowCover, df_SnowCover_new)
          
          # #For debugging
          #  ggplot() + geom_point(data=df_SnowCover, aes(x=doy, y=SnowFraction, col=NDSI_threshold)) + theme_classic()  
          
          # #Save dataframe for current NDSI_threshold
          #  write.csv(df_SnowCover_new, paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Data_SnowFraction.csv"), row.names = FALSE)
        
          }
      
      #Save dataframe with snowfraction data for all NDSI thresholds
      write.csv(df_SnowCover, paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Data_SnowFraction.csv"), row.names = FALSE)
      
  #(26): Fit statistical models (GAMS) through the extracted snowcover data
    
      #Empty dataframe for storing RAW SnowFraction data including detected outliers
      df_SnowFraction_GAM <- data.frame(NDSI_threshold=character(),
                                        SnowFraction=numeric(),
                                        Date=factor(),
                                        doy=numeric(),
                                        outliers=logical())
      
      #Empty dataframe for storing GAM predictions for SnowFraction data:
      df_SnowFraction_GAM_predictions <- data.frame(NDSI_threshold=character(),
                                                    doy=numeric(), 
                                                    SnowFraction_gam_predict=numeric(), 
                                                    stringsAsFactors=FALSE)
      
      #Loop through all NDSI_thresholds and fit GAMS through the snowfraction data
      for(j in unique(df_SnowCover$NDSI_threshold)){
        
        #For debugging
        #j=unique(df_SnowCover$NDSI_threshold)[1]
        
        #Select NDSI-threshold-specific subset of data:
        df_SnowFraction_GAM_new <- df_SnowCover[df_SnowCover$NDSI_threshold==j & !is.na(df_SnowCover$SnowFraction),
                                                c("NDSI_threshold", "SnowFraction", "Date", "doy")]
        
        #Fit a gam through the NDSI-threshold-specific SnowFraction ~ doy data and employ sequential outlier removal
        df_SnowFraction_GAM_new <- f_gam_SeqRemOutliers(data=df_SnowFraction_GAM_new, y="SnowFraction", x="doy", outlier_removal=outlier_removal, 
                                                        outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                        default_k=gam_k_outlier)
        
        #Sort df_SnowFraction_GAM_new by doy:
        df_SnowFraction_GAM_new <- df_SnowFraction_GAM_new[order(df_SnowFraction_GAM_new$doy),]
        
        #Bind the NDSI-threshold-specific dataframe detected outliers to the dataframe containing all dataframes from previous iterations
        df_SnowFraction_GAM <- rbind(df_SnowFraction_GAM, df_SnowFraction_GAM_new)
        
        #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot smoother curves
        
          #Refit GAM through data without outliers
          index <- which(df_SnowFraction_GAM_new$outliers==FALSE)
          mod_gam <- with(df_SnowFraction_GAM_new[index,], mgcv::gam(SnowFraction ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))
          
          #Use gam to make predictions on a more detailed (1-day) day of year interval
          aoi_SnowFraction_predicted_new <- data.frame(NDSI_threshold=j, doy=seq(min(df_SnowFraction_GAM_new$doy), max(df_SnowFraction_GAM_new$doy), 1))
          aoi_SnowFraction_predicted_new$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=aoi_SnowFraction_predicted_new, type="response")
          aoi_SnowFraction_predicted_new <- aoi_SnowFraction_predicted_new[order(aoi_SnowFraction_predicted_new$doy),]
          aoi_SnowFraction_predicted_new$year <- year_ID
          
          #Add predictions to df_SnowFraction_GAM_predictions dataframe:  
          df_SnowFraction_GAM_predictions <- rbind(df_SnowFraction_GAM_predictions, aoi_SnowFraction_predicted_new)
        
          }
      
      #Save the SnowFraction data and GAM predictions per NDSI-threshold as a .csv file
      #write.csv(df_SnowFraction_GAM, paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Data_SnowFraction_GAM.csv"), row.names = FALSE)
      write.csv(df_SnowFraction_GAM_predictions, paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_GAM_Predictions_SnowFraction.csv"), row.names = FALSE)
      
      #Plot the raw SnowFraction datapoints and gam predictions for each NDSI-threshold:
      
        #Plot SnowFraction and model predictions for all NDSI-thresholds in a single plot
        p_SnowFraction = ggplot()+
         geom_point(data=df_SnowFraction_GAM, aes(x=doy, y=SnowFraction, fill=NDSI_threshold, col=NDSI_threshold))+
         geom_line(data=df_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict, col=NDSI_threshold)) +
         xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
         ylab(paste0("Snowcover fraction within study area in ", year_ID)) +
         theme_classic()
          
        #Plot SnowFraction and model predictions in a separate plot per NDSI-threshold
        p_SnowFraction_grid = ggplot()+ 
          geom_point(data=df_SnowFraction_GAM[df_SnowFraction_GAM$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
          geom_point(data=df_SnowFraction_GAM[df_SnowFraction_GAM$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
          geom_line(data=df_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict), col="#1620de", lwd=1)+
          geom_vline(xintercept = 150, colour="grey", lty=2)+
          geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
          facet_wrap(~NDSI_threshold, ncol=ceiling(length(unique(df_SnowFraction_GAM$NDSI_threshold))^0.5))+
          xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
          ylab(paste0("Snowcover fraction within study area in ", year_ID)) +
          theme_classic()
       
   #(27): Calculate at which day of year the predicted snowfraction is equal to the levels of Snowfraction_threshold_vector
        
        #Setup parallel processing
        numCores <- detectCores()
        cl <- makePSOCKcluster(numCores)
        registerDoSNOW(cl)
        
        #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector
        results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                    pixelIDs_split = list(NDSI_threshold_vector), #levels of NDSI_threshold (input needs to be a list)
                                                    df_pixel_y = df_SnowFraction_GAM_predictions, #dataframe containing GAM predictions
                                                    pixel_ID_column="NDSI_threshold", #Grouping column
                                                    y="SnowFraction_gam_predict", #response variable in GAM
                                                    x="doy", #predictor variable in GAM
                                                    pixel_gam_plots = FALSE, #Should GAM plots be created
                                                    y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated 
        
        #Turn parallel processing off and run sequentially again after this point
        stopCluster(cl)
        registerDoSEQ()
          
        #Extract dates of snowmelt per NDSI-threshold and per SnowFraction-threshold
        df_SnowFraction <- results[[1]]
        df_SnowFraction <- as.data.frame(do.call(rbind, df_SnowFraction))
        colnames(df_SnowFraction)[colnames(df_SnowFraction)=="pixel_ID"] <- "NDSI_threshold"
        colnames(df_SnowFraction)[colnames(df_SnowFraction)=="x_threshold"] <- "doy"
        colnames(df_SnowFraction)[colnames(df_SnowFraction)=="y_threshold"] <- "Snowfraction_threshold"
        
        #Save dates of snowmelt per NDSI-threshold per SnowFraction-threshold as a .csv file
        write.csv(df_SnowFraction, paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Snowmelt_Snowfraction.csv"), row.names = FALSE)
        
        #Add dates of snowmelt to the plot 'p_SnowFraction_grid'   
        p_SnowFraction_Snowmelt_grid <- p_SnowFraction_grid +
          geom_point(data=df_SnowFraction[!is.na(df_SnowFraction$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)
        
        ggsave(plot=p_SnowFraction_Snowmelt_grid, paste0(here(), "/Output/S2/05_Shapefile_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Plot_Snowfraction.pdf"), width=12, height = 10)
        

##################################################################################################################################       
##################################################################################################################################       
##################################################################################################################################       
             
      