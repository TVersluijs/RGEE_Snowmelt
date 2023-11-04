##################################################################################################################################

#Use Sentinel-2 data to extract time series of the average NDVI, NDMI, and NDSI and the fraction of snowcover for all sub areas 
#located within a shapefile (these sub areas can be specified by creating a multipolygon in e.g. QGIS). The fraction of snowcover 
#is estimated by calculating the fraction of pixels per subarea with an NDSI value larger than the user specified NDSI-threshold. 
#This corresponds to the method='snowfraction' in other scripts.

#Copyright Tom Versluijs 2023-11-01. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R". 
#Note that a GoogleDrive is required.

##################################################################################################################################

#I: Setup workspace

##################################################################################################################################

      #(0) Clear workspace and set python environment
       rm(list=ls())
       utils::install.packages("here")
       library(here)
       if(file.exists(paste0(here::here(), "/Input/rgee_environment_dir.rds"))){
         rgee_environment_dir <- readRDS(paste0(here::here(), "/Input/rgee_environment_dir.rds"))
         reticulate::use_python(rgee_environment_dir, required=T)
         reticulate::py_config()
         }
        
      #(1) Load packages
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
       
 #Manually specify parameters of interest

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

     #Name of Shapefiles for each subArea located in Input folder
	 #Follow the guide "Manual_CreateShapefilePolygons.docx" when creating this shapefile.
     shapefile_subareas <- "ZAC_TenEqualSizeVoronoiPolygons_EPSG4326.shp"

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

     #NDSI threshold for snow detection (specify multiple using c())
     NDSI_threshold_vector = c(0.42, 0.3, 0.5)
     
     #Define the snowcover fraction for which the date of its occurrence will be calculated (specify multiple using c())
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
     #First, all cloudy pixels without cloudy pixels within radius cld_erosion will be removed.
     #Second, all remaining pixels will be dilated by radius cld_buffer. Make sure that cld_buffer
     #is larger than cld_erosion. Larger values result in a more course cloud-filtering process.
     cld_erosion=250 #250 for conservative cloud estimation
     cld_buffer=270 #270 for conservative cloud estimation

     #Define resolution of cloud-filtering (in meters)
     resolution_cldmsk=150 #(default = 150)

     #Maximum fraction of cloud cover allowed in each image
     max_cloud_fraction=0.75

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
     
     #Create a timestamp variable
     timestamp <- format(Sys.time(), "%Y%m%d%H%m%S")
     
     #Set mask_clouds to TRUE if a composite image needs to be generated         
     if(create_composite==TRUE){mask_clouds=TRUE}      
     
     #Set mask_clouds to TRUE if water masking is TRUE and mask_water_type==water_mask_Manual
     if(mask_water==TRUE & (mask_water_type=="water_mask_Manual" | mask_water_type=="both")){mask_clouds=TRUE}      
     
     #Create output folder
     if(dir.exists(paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt"), recursive = TRUE)}
     
       
##################################################################################################################################
        
#III: Read and display the unfiltered data
        
##################################################################################################################################

     #(6): Read study area shapefile and convert to a feature collection.
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
       paste0('Size of study area calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')
     
    #(7) Extract satellite data including cloud probability band:
      
      #Extract Sentinel-2 Surface Reflectance data for a given area and date range:
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
       
      #Clip all images to the area depicted by 'aoi_Shapefile'
       s2_col <- s2_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))}) 
     
   #(8): Plot an RGB, NDSI and NDVI image for a single extracted satellite image (for debugging)   
         
      # #Select first images for initial plot:
      #  image <- s2_col$filterDate(paste0(year_ID, "-07-15"), end_date)$first()
      # 
      # #Normalized difference snow index:
      #  ndsi_s2 <-  image$expression('(B3 - B11) / (B3 + B11)',
      #                               list('B11'=image$select('B11'),
      #                                    'B3'=image$select('B3')))
      # 
      # #Normalized difference vegetation index:
      #  ndvi_s2 <- image$expression('(B8 - B4) / (B8 + B4)',
      #                              list('B8'=image$select('B8'),
      #                                   'B4'=image$select('B4')))
      # 
      # #Add Normalized difference water index to this image:
      #  ndwi_s2 <- image$expression('(B3 - B8) / (B3 + B8)',
      #                              list('B3'=image$select('B3'),
      #                                   'B8'=image$select('B8')))
      #  
      # #Plot all layers
      #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
      #  Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
      #  Map$addLayer(ndsi_s2,list(min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
      #  Map$addLayer(ndvi_s2,list(min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
      #  Map$addLayer(ndwi_s2,list(min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')
      
    #(9): Load shapefile containing 10 subareas shown as voronoi polygons around centroids of 10-mean centering of 10000 
    #random points within the shapefile area
      
      #Read shapefiles
       root_fldr <- here()
       aoi_SubAreas <- st_read(paste0(root_fldr, "/Input/Shapefiles/", shapefile_subareas), quiet=T)
      
      #Plot the 10 Sub areas within the shapefile area:
       p1 <- ggplot() + 
         geom_sf(data = aoi_SubAreas, fill=sf.colors(10), col = "black")+
         theme_tom()+
         geom_sf_label(data = aoi_SubAreas, aes(label=Cluster_ID), colour="black", 
                       nudge_x =c(0.001,0,0.004,0,-0.001,0,-0.002,0,0.001,0), 
                       nudge_y =c(0.003,0.002,0,0,0,0.005,0,-0.001,-0.001,-0.001))
       
         ggsave(plot=p1, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_SubAreas.pdf"), width=10, height=8)
      
      #Convert shapefile to an earthengine feature collection:  
       aoi_SubAreas <- st_transform(aoi_SubAreas, crs="EPSG:4326")
       aoi_SubAreas <- sf_as_ee(aoi_SubAreas)
       aoi_SubAreas <- ee$FeatureCollection(aoi_SubAreas)
       #Note that the feature$property "Cluster_ID" has a unique number (0-9) for each of the subareas within the study area 
      
      #Plot featurecollection
       image <- s2_col$first()
       Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
       Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
       Map$addLayer(aoi_SubAreas, list(color="red"), 'SubAreas')  
      
    #(10) Add NDSI and NDVI to the clipped image collection

      #Map Normalized difference functions over the image collection
      s2_col <- s2_col$
        map(getNDSI)$
        map(getNDVI)$
        map(getNDMI)$
        map(getNDWI)
    
    #(11) Create a timeseries gif of RGB images for the shapefile area
      
      # #Check number of images in collection
      # s2_col$size()$getInfo()
      # 
      # #Create a timelapse video of RGB band
      # videoArgs <- list(dimensions=380, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7))
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
         
    #(12) Add cloud information within the shapefile area to the image collection: 
 
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
      
      # #Check if CloudFraction property has been added to each image
      # s2_col$first()$propertyNames()$getInfo()
      
      # #Visual check of cloud mask (for debugging)
      # image <- s2_col$filterDate(paste0(year_ID, "-05-31"), end_date)$first()$
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
      #area of interest only called 'CloudFraction'. Below we compare both measures, as we expect they should be correlated. Differences
      #between both measures can for example arise when the shapefile area was covered by clouds while most of the larger region (tile)
      #was not or vice versa.
      
      # #Extract cloud fraction of all images in image collection for the shapefile area
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
      #   #Thus, this indicates that our local cloud measure does not necessarily correspond to that as measured on the whole tile as
      #   #explained above.
    
    #(13): Exclude all images from the image collection that have a CloudFraction value > max_cloud_fraction, and mask all cloud pixels within the remaining images
  
      #In the analysis for the whole study area we exclude images with more than 20% snowcover within our shapefile
      #area. However, now we're interested in individual pixels. So even when 90% of the shapefile area is 
      #snowcovered, the pixels we're interested in might be located in the 10% cloud free area. We thus only filter
      #images with a 100% cloudcover (cloudFraction=1.0). We apply a cloud mask to mask out cloudy pixels within each 
      #image.
 
      #Mask cloudy pixels within each image: 
      s2_clouds_filtered <- s2_col$
        #Only remove images that are fully cloud-covered
        filter(ee$Filter$lt('CloudFraction', max_cloud_fraction))$
        #Apply cloudmask for individual pixels
        map(Add_CloudMask)
      
      # #Create timelapse video of the cloud filtered/masked RGB images
      # videoArgs <- list(dimensions=350, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("B4", "B3", "B2"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
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
                                          start_date_NDSI=start_date_NDSI, end_date_NDSI=end_date_NDSI, NDSI_threshold=NDSI_threshold_vector[1], NIR_threshold=NIR_threshold)
            
            
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

#VI: Create a composite image for each day by mosaicking all images from that day

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
     
#VI: Extract Average NDSI, NDVI, NDMI and SnowFraction per SubArea      
      
##################################################################################################################################

  #Run the code for all levels of NDSI_threshold
      
    #Specify the NDSI threshold(s) at which an area is perceived as snow-free
     NDSI_threshold_vector = NDSI_threshold_vector
    
    #Create an empty dataframe for storing output per SubArea
     df_SubAreas_BandValues <- data.frame(NDSI=numeric(),
                                          NDVI=numeric(),
                                          NDMI=numeric(),
                                          SnowFraction=numeric(),
                                          Date=character(),
                                          doy=numeric(),
                                          SubArea=character(),
                                          NDSI_threshold=character())
      
    #Run the analysis for each level of NDSI_threshold_vector
     for(NDSI_threshold in NDSI_threshold_vector){  
        
        #Print message
        print(paste0("  -START ANALYSIS FOR NDSI_threshold = ", NDSI_threshold))
        
        #Store NDSI_threshold as a character (used for naming of outputs)
        NDSI_threshold_char <- gsub("\\.", "_", as.character(NDSI_threshold))
      
        #(20): Add binary SNOW-band to the image collection that identifies for each pixel whether it is snowcovered (NDSI > NDSI_threshold)
        
            #Map computeSnow function over the image collection  
            s2_snow_masked <- s2_col_composite$
                map(computeSnow)  
      
            # #Visualize the binary layer (for debugging).
            #  image <- s2_snow_masked$filterDate(paste0(year_ID, "-06-15"), end_date)$first()
            #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
            #  Map$addLayer(image,list(bands=c("B4", "B3", "B2"), min=0, max=10000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
            #  Map$addLayer(image$select('SNOW'))
          
        #(21): Extract Sentinel-2 band values (NDSI, NDVI, NDMI and Snow Fraction) for each sub area within the study area over time
          
            #Create an iteration function that we will use to iterate through all images of the image collection. For each image, the
            #value of certain bands is extracted for each feature (i.e. SubArea) within the feature collection aoi_SubAreas. The 
            #resulting band values (and the datetime of the image) are added as a property to each feature (i.e. SubArea). This results
            #in an updated feature collection specific for the current image. This feature collection is then appended to a list
            #of feature collections from previous image iterations. After iterating through n-images, the result is thus a feature
            #collection list of length n * the number of features within the feature collection.
            
            #Thus, at each iteration we extract band values of interest for all SubAreas within the current image, store this as a  
            #feature collection and add this to an expanding list of feature collections from previous iterations.   
            
            #Store default Sentinel-2 image projection 
             S2Projection <- s2_snow_masked$first()$select("NDSI")$projection()
             #S2Projection$getInfo() 
              
            #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
             FC_initial <- ee$FeatureCollection(ee$List(list()))
            
            #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
            #(in this case the current iteration image) and the second element takes the output value from the iteration that preceeded it. The latter
            #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with 
            #has been defined.
             Extract_BandValuesAtSubAreas = Extract_BandValuesAtSubAreas #sourced
            
            #Iterate over the ImageCollection 
             FC_merged <- ee$FeatureCollection(s2_snow_masked$select("NDSI", "NDVI", "NDMI", "SNOW")$iterate(Extract_BandValuesAtSubAreas, FC_initial))
             #FC_merged$getInfo() #for debugging
             #FC_merged$first()$getInfo() #for debugging
        
        #(22): Transform feature collection with Sentinel-2 Band values for each sub area to a dataframe     
        
            #We export the data instead of using aggregate_array() as the latter might fail due to computation timeouts.
             
            #Setup task
             task_vector <- ee_table_to_drive(
               collection= FC_merged,
               description = paste0(timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_MeanBandValues_SubAreas"),
               folder="RGEE_tmp",
               fileFormat="CSV",
               selectors=c('NDSI', 'NDVI', 'NDMI', 'SNOW', 'Date', 'Cluster_ID')
               )
             
            #Run and monitor task
             task_vector$start() 
             ee_monitoring(task_vector, task_time=60, max_attempts=1000000)
             
            #Import results
             exported_stats <- ee_drive_to_local(task=task_vector, dsn=paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_MeanBandValues_SubAreas"))
             df_SubAreas_BandValues_new <- read.csv(exported_stats)
             unlink(exported_stats)
             
            #Restructure dataframe
             df_SubAreas_BandValues_new$doy <- as.numeric(strftime(df_SubAreas_BandValues_new$Date, format = "%j"))
             colnames(df_SubAreas_BandValues_new) <- c('NDSI', 'NDVI', 'NDMI', 'SnowFraction', 'Date', 'SubArea', 'doy')         
             df_SubAreas_BandValues_new <- df_SubAreas_BandValues_new[ ,c('NDSI', 'NDVI', 'NDMI', 'SnowFraction', 'Date', 'doy', 'SubArea')]
    
            #Change -9999 to NA
             df_SubAreas_BandValues_new$NDSI[df_SubAreas_BandValues_new$NDSI < -9000] <- NA
             df_SubAreas_BandValues_new$NDVI[df_SubAreas_BandValues_new$NDVI < -9000] <- NA
             df_SubAreas_BandValues_new$NDMI[df_SubAreas_BandValues_new$NDMI < -9000] <- NA
             df_SubAreas_BandValues_new$SnowFraction[df_SubAreas_BandValues_new$SnowFraction < -9000] <- NA
             df_SubAreas_BandValues_new$doy[df_SubAreas_BandValues_new$doy < -9000] <- NA
             df_SubAreas_BandValues_new$SubArea[df_SubAreas_BandValues_new$SubArea < -9000] <- NA
    
            #Sort dataframe by ClusterID and doy
             index <- with(df_SubAreas_BandValues_new, order(SubArea, doy))
             df_SubAreas_BandValues_new <- df_SubAreas_BandValues_new[index,]
            
            #Add NDSI_threshold as a new column
             df_SubAreas_BandValues_new$NDSI_threshold <- as.factor(NDSI_threshold)
             
            #Add dataframe for current subarea to dataframe from previous iterations:
             df_SubAreas_BandValues <- rbind(df_SubAreas_BandValues, df_SubAreas_BandValues_new)
             
             # #For debugging
             #  ggplot() + geom_point(data=df_SubAreas_BandValues, aes(x=doy, y=SnowFraction, col=NDSI_threshold)) + theme_classic()  
             
             # #Save dataframe for current subarea
             #  write.csv(df_SubAreas_BandValues_new, paste0("Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, _Resolution", resolution, "_NDSI", NDSI_threshold_char, "_Data_MeanBandValues_SnowFraction_SubAreas.csv"), row.names = FALSE)
           
          }
     
     #Store dataframe with average bandValues and Snowfraction for all subareas for all levels of NDSI_threshold
     write.csv(df_SubAreas_BandValues, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_MeanBandValues_SnowFraction.csv"), row.names = FALSE)
     #df_SubAreas_BandValues <- read.csv(paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_MeanBandValues_SnowFraction.csv"), header=T)
     
##################################################################################################################################

#VII: Fit GAMS through the Average NDSI, NDVI, NDMI and SnowFraction data per SubArea      

##################################################################################################################################

   #(23): Fit GAMS 
      
      #We fit a Generalized Additive Model (GAM) through the data within each sub-area. We do this using a sequential
      #process. We first fit a GAM through the data, calculating model predictions and residuals. We then exclude all 
      #rows from the dataframe where the residual >= (0.25 * the range of the data). We then re-fit a GAM to this reduced 
      #dataset, make predictions and calculate residuals. We then exclude all rows from the reduced dataframe where the
      #residual >= >= (0.1 * the range of the data). This gives us a final dataframe through which we fit a third
      #GAM and store its model predictions given the dataset in which outliers were thus removed through two subsequent 
      #steps. This whole process is executed using the function f_gam_SeqRemOutliers. Note that a sequential step is
      #required because initially some datapoints might falsely be assigned a large residual because of one extreme
      #outlier. After removal of this extreme outlier and refitting of the GAM, it can be better assessed which data
      #points truly have a large residual and can thus be assigned as 'true' outliers.
        f_gam_SeqRemOutliers <- f_gam_SeqRemOutliers #sourced
       
        #Specify sequential outlier thresholds:
         outlier_thresh_1=outlier_thresh_1
         outlier_thresh_2=outlier_thresh_2
         outlier_removal=outlier_removal
         
        #Specify at which NDSI value, the SubArea is considered snow-free
         NDSI_threshold_vector = NDSI_threshold_vector
         
########################################################################################################################################################################################
         
  #(A): SNOWFRACTION - Fit a Generalized Additive Model (GAM) through the snow fraction data within each sub-area  
  
      #(A.1) Create an empty dataframe
        df_SubAreas_SnowFraction_GAM <- data.frame(NDSI_threshold=character(),
                                                   SnowFraction=numeric(),
                                                   Date=factor(),
                                                   doy=numeric(),
                                                   SubArea=factor(),
                                                   outliers=logical())
       
        df_SubAreas_SnowFraction_GAM_predictions <- data.frame(SubArea=character(),
                                                               NDSI_threshold=character(),
                                                               doy=numeric(), 
                                                               SnowFraction_gam_predict=numeric(), 
                                                               stringsAsFactors=FALSE)   
    
      #(A.2) Loop through all SubAreas and fit a separate gam (with sequential outlier removal) to the SubArea-specific SnowFraction data
        for(i in unique(df_SubAreas_BandValues$SubArea)){
          
          #For debugging  
          #i=unique(df_SubAreas_BandValues$SubArea)[1]
         
          #Loop through all NDSI_thresholds
          for(j in unique(df_SubAreas_BandValues$NDSI_threshold)){
            
            #For debugging
            #j=unique(df_SubAreas_BandValues$NDSI_threshold)[1]
          
            #Select SubArea-specific subset of data:
            df_SubArea_SnowFraction_GAM_new <- df_SubAreas_BandValues[df_SubAreas_BandValues$SubArea==i & 
                                                               df_SubAreas_BandValues$NDSI_threshold==j &
                                                               !is.na(df_SubAreas_BandValues$SnowFraction),
                                                               c("NDSI_threshold", "SnowFraction", "Date", "doy", "SubArea")]
                       
            #Fit a gam through the SubArea-specific SnowFraction ~ doy data and emply sequential outlier removal
            df_SubArea_SnowFraction_GAM_new <- f_gam_SeqRemOutliers(data=df_SubArea_SnowFraction_GAM_new, y="SnowFraction", x="doy", outlier_removal=outlier_removal, 
                                                             outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2, 
                                                             default_k=gam_k_outlier)
          
            #Sort df_SubArea_SnowFraction_GAM_new by doy:
            df_SubArea_SnowFraction_GAM_new <- df_SubArea_SnowFraction_GAM_new[order(df_SubArea_SnowFraction_GAM_new$doy),]
             
            #Bind the SubArea-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
             df_SubAreas_SnowFraction_GAM <- rbind(df_SubAreas_SnowFraction_GAM, df_SubArea_SnowFraction_GAM_new)
          
            #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves
           
              #Refit GAM through data
               index <- which(df_SubArea_SnowFraction_GAM_new$outliers==FALSE)
               mod_gam <- with(df_SubArea_SnowFraction_GAM_new[index,], mgcv::gam(SnowFraction ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))
           
              #Use gam to make predictions on a more detailed (1-day) day of year interval
               df_SubArea_SnowFraction_GAM_predictions_new <- data.frame(SubArea=i, NDSI_threshold=j, doy=seq(min(df_SubArea_SnowFraction_GAM_new$doy), max(df_SubArea_SnowFraction_GAM_new$doy), 1))
               df_SubArea_SnowFraction_GAM_predictions_new$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=df_SubArea_SnowFraction_GAM_predictions_new, type="response")
               df_SubArea_SnowFraction_GAM_predictions_new <- df_SubArea_SnowFraction_GAM_predictions_new[order(df_SubArea_SnowFraction_GAM_predictions_new$doy),]
              
              #Add predictions to df_SubAreas_SnowFraction_GAM_predictions dataframe:  
               df_SubAreas_SnowFraction_GAM_predictions <- rbind(df_SubAreas_SnowFraction_GAM_predictions, df_SubArea_SnowFraction_GAM_predictions_new)
           
             }
          
          }
       
        #Change subArea column to factor
        df_SubAreas_SnowFraction_GAM$SubArea <- as.factor(as.character(df_SubAreas_SnowFraction_GAM$SubArea))
        df_SubAreas_SnowFraction_GAM_predictions$SubArea <- as.factor(as.character(df_SubAreas_SnowFraction_GAM_predictions$SubArea))
        
        #Save dataframe with GAM fits for SnowFraction
        #write.csv(df_SubAreas_SnowFraction_GAM, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_SnowFraction_GAM.csv"), row.names = FALSE)
        write.csv(df_SubAreas_SnowFraction_GAM_predictions, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_SnowFraction_Predictions_GAM.csv"), row.names = FALSE)
          
      #(A.3) Plot the raw Snowfraction datapoints and gam predictions for each SubArea:
    
         #Plot Snowfraction and model predictions for all SubAreas in a single plot
          p_SubArea_SnowFraction = ggplot()+
            geom_point(data=df_SubAreas_SnowFraction_GAM, aes(x=doy, y=SnowFraction, fill=SubArea, col=SubArea))+
            geom_line(data=df_SubAreas_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict, col=SubArea)) +
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("Fraction of snow-covered pixels") +
            facet_wrap(~NDSI_threshold, ncol=3)+
            theme_tom()
                
         #Plot SnowFraction and model predictions in a separate plot per SubArea
          p_SubArea_SnowFraction_grid = ggplot()+ 
            geom_point(data=df_SubAreas_SnowFraction_GAM[df_SubAreas_SnowFraction_GAM$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
            geom_point(data=df_SubAreas_SnowFraction_GAM[df_SubAreas_SnowFraction_GAM$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
            geom_line(data=df_SubAreas_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict), col="#1620de", lwd=1.25)+
            geom_vline(xintercept = 150, colour="grey", lty=2)+
            geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
            facet_wrap(~SubArea + NDSI_threshold, ncol=ceiling(length(unique(df_SubAreas_SnowFraction_GAM$SubArea))^0.5))+
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("Fraction of snow-covered pixels") +
            theme_tom()
           
      #(A.4): Calculate at which day of year the SnowFraction value is equal to snowFraction_threshold for each SubArea using predictions from mod_gam
          
          #Create an empty dataframe to store snowfraction dates
          df_SubArea_Snowfraction <- data.frame(NDSI_threshold=character(),
                                                Snowfraction_threshold=character(),
                                                doy=numeric(), 
                                                SubArea=character(), 
                                                stringsAsFactors=FALSE)   
          
          #Setup parallel processing
          numCores <- detectCores()
          cl <- makePSOCKcluster(numCores)
          registerDoSNOW(cl)
          
          #Loop through all SubAreas in the df_SubAreas_SnowFraction_GAM_predictions dataframe and extract snowfraction dates for each SubArea:
          for(i in unique(df_SubAreas_SnowFraction_GAM_predictions$SubArea)){
            
            #For debugging  
            #i=unique(df_SubAreas_SnowFraction_GAM_predictions$SubArea)[1] #for debugging
            
            #Select dataset with GAM predictions for current SubArea
            df_SubArea_gam <- df_SubAreas_SnowFraction_GAM_predictions[df_SubAreas_SnowFraction_GAM_predictions$SubArea==i & !is.na(df_SubAreas_SnowFraction_GAM_predictions$SnowFraction_gam_predict),]
            
            #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector
            results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                        pixelIDs_split = list(NDSI_threshold_vector), #levels of NDSI_threshold (input needs to be a list)
                                                        df_pixel_y = df_SubArea_gam, #dataframe containing GAM predictions
                                                        pixel_ID_column="NDSI_threshold", #Grouping column
                                                        y="SnowFraction_gam_predict", #response variable in GAM
                                                        x="doy", #predictor variable in GAM
                                                        pixel_gam_plots = FALSE, #Should GAM plots be created
                                                        y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated 
            
            #Store dates of snowmelt per SubArea
            df_SubArea_SnowFraction_new <- results[[1]]
            df_SubArea_SnowFraction_new <- as.data.frame(do.call(rbind, df_SubArea_SnowFraction_new))
            colnames(df_SubArea_SnowFraction_new)[colnames(df_SubArea_SnowFraction_new)=="pixel_ID"] <- "NDSI_threshold"
            colnames(df_SubArea_SnowFraction_new)[colnames(df_SubArea_SnowFraction_new)=="x_threshold"] <- "doy"
            colnames(df_SubArea_SnowFraction_new)[colnames(df_SubArea_SnowFraction_new)=="y_threshold"] <- "Snowfraction_threshold"
            df_SubArea_SnowFraction_new$SubArea <- i
            
            #Add snowmelt data for current SubArea to general dataframe
            df_SubArea_Snowfraction <- rbind(df_SubArea_Snowfraction, df_SubArea_SnowFraction_new)
            
          }
            
          #Turn parallel processing off and run sequentially again after this point
          stopCluster(cl)
          registerDoSEQ()
            
          #Change column subArea to a factor:
          df_SubArea_Snowfraction$SubArea <- as.factor(as.character(df_SubArea_Snowfraction$SubArea))
          
          #Save dates of snowmelt per subarea per SnowFraction threshold as a .csv file
          write.csv(df_SubArea_Snowfraction, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Snowmelt_Snowfraction.csv"), row.names = FALSE)
          
          #Add dates of snowmelt to the plot 'p_SubArea_SnowFraction_grid'   
          p_SubArea_SnowFraction_Snowmelt_grid <- p_SubArea_SnowFraction_grid +
            geom_point(data=df_SubArea_Snowfraction[!is.na(df_SubArea_Snowfraction$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)
          
          ggsave(plot=p_SubArea_SnowFraction_Snowmelt_grid, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Plot_Snowmelt_Snowfraction_grid.pdf"), width=16, height = 14)
        
########################################################################################################################################################################################
              
  #(B): NDSI - Fit a Generalized Additive Model (GAM) through the NDSI data within each sub-area

      #(B.1) Create an empty dataframe
        df_SubAreas_NDSI_GAM <- data.frame(NDSI=numeric(),
                                           Date=factor(),
                                           doy=numeric(),
                                           SubArea=factor(),
                                           outliers=logical())

        df_SubAreas_NDSI_GAM_predictions <- data.frame(SubArea=character(),
                                                       doy=numeric(),
                                                       NDSI_gam_predict=numeric(),
                                                       stringsAsFactors=FALSE)

      #(B.2) Loop through all SubAreas and fit a separate gam (with sequential outlier removal) to the SubArea-specific NDSI data
        for(i in unique(df_SubAreas_BandValues$SubArea)){

          #For debugging
          #i=unique(df_SubAreas_BandValues$SubArea)[1]

          #Run code for the first NDSI_threshold (as NDSI-data is collected independently of thresholds)
          j=unique(df_SubAreas_BandValues$NDSI_threshold)[1]

          #Select SubArea-specific subset of data:
          df_SubArea_NDSI_GAM_new <- df_SubAreas_BandValues[df_SubAreas_BandValues$SubArea==i &
                                                            df_SubAreas_BandValues$NDSI_threshold==j &
                                                            !is.na(df_SubAreas_BandValues$NDSI),
                                                            c("NDSI", "Date", "doy", "SubArea")]

          #Fit a gam through the SubArea-specific NDSI ~ doy data and emply sequential outlier removal
          df_SubArea_NDSI_GAM_new <- f_gam_SeqRemOutliers(data=df_SubArea_NDSI_GAM_new, y="NDSI", x="doy", outlier_removal=outlier_removal,
                                                           outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                           default_k=gam_k_outlier)

          #Sort df_SubArea_NDSI_GAM_new by doy:
          df_SubArea_NDSI_GAM_new <- df_SubArea_NDSI_GAM_new[order(df_SubArea_NDSI_GAM_new$doy),]

          #Bind the SubArea-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
           df_SubAreas_NDSI_GAM <- rbind(df_SubAreas_NDSI_GAM, df_SubArea_NDSI_GAM_new)

          #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

            #Refit GAM through data
             index <- which(df_SubArea_NDSI_GAM_new$outliers==FALSE)
             mod_gam <- with(df_SubArea_NDSI_GAM_new[index,], mgcv::gam(NDSI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

            #Use gam to make predictions on a more detailed (1-day) day of year interval
             df_SubArea_NDSI_GAM_predictions_new <- data.frame(SubArea=i, doy=seq(min(df_SubArea_NDSI_GAM_new$doy), max(df_SubArea_NDSI_GAM_new$doy), 1))
             df_SubArea_NDSI_GAM_predictions_new$NDSI_gam_predict <- stats::predict(mod_gam, newdata=df_SubArea_NDSI_GAM_predictions_new, type="response")
             df_SubArea_NDSI_GAM_predictions_new <- df_SubArea_NDSI_GAM_predictions_new[order(df_SubArea_NDSI_GAM_predictions_new$doy),]

            #Add predictions to df_SubAreas_NDSI_GAM_predictions dataframe:
             df_SubAreas_NDSI_GAM_predictions <- rbind(df_SubAreas_NDSI_GAM_predictions, df_SubArea_NDSI_GAM_predictions_new)

          }

        #Change subArea column to factor
        df_SubAreas_NDSI_GAM$SubArea <- as.factor(as.character(df_SubAreas_NDSI_GAM$SubArea))
        df_SubAreas_NDSI_GAM_predictions$SubArea <- as.factor(as.character(df_SubAreas_NDSI_GAM_predictions$SubArea))

        #Save dataframe with GAM fits for NDSI
        #write.csv(df_SubAreas_NDSI_GAM, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_NDSI_GAM.csv"), row.names = FALSE)
        write.csv(df_SubAreas_NDSI_GAM_predictions, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_NDSI_Predictions_GAM.csv"), row.names = FALSE)

      #(B.3) Plot the raw NDSI datapoints and gam predictions for each SubArea:

         #Plot NDSI and model predictions for all SubAreas in a single plot
          p_SubArea_NDSI = ggplot()+
            geom_point(data=df_SubAreas_NDSI_GAM, aes(x=doy, y=NDSI, fill=SubArea, col=SubArea))+
            geom_line(data=df_SubAreas_NDSI_GAM_predictions, aes(x=doy, y=NDSI_gam_predict, col=SubArea)) +
            geom_hline(yintercept = NDSI_threshold_vector, colour="grey", lty=2)+
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("NDSI") +
            theme_tom()
            
         #Plot NDSI and model predictions in a separate plot per SubArea
          p_SubArea_NDSI_grid = ggplot()+
            geom_point(data=df_SubAreas_NDSI_GAM[df_SubAreas_NDSI_GAM$outliers==FALSE,], aes(x=doy, y=NDSI))+
            geom_point(data=df_SubAreas_NDSI_GAM[df_SubAreas_NDSI_GAM$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
            geom_line(data=df_SubAreas_NDSI_GAM_predictions, aes(x=doy, y=NDSI_gam_predict), col="#1620de", lwd=1.25)+
            geom_vline(xintercept = 150, colour="grey", lty=2)+
            geom_hline(yintercept = NDSI_threshold_vector, colour="grey", lty=2)+
            facet_wrap(~SubArea, ncol=ceiling(length(unique(df_SubAreas_NDSI_GAM$SubArea))^0.5))+
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("NDSI") +
            theme_tom()
           
      #(B.4): Calculate at which day of year the NDSI value is equal to NDSI_threshold for each SubArea using predictions from mod_gam

          #Setup parallel processing
          numCores <- detectCores()
          cl <- makePSOCKcluster(numCores)
          registerDoSNOW(cl)
          
          #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in NDSI_threshold_vector
          results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                      pixelIDs_split = list(unique(df_SubAreas_NDSI_GAM_predictions$SubArea)), #levels of NDSI_threshold (input needs to be a list)
                                                      df_pixel_y = df_SubAreas_NDSI_GAM_predictions, #dataframe containing GAM predictions
                                                      pixel_ID_column="SubArea", #Grouping column
                                                      y="NDSI_gam_predict", #response variable in GAM
                                                      x="doy", #predictor variable in GAM
                                                      pixel_gam_plots = FALSE, #Should GAM plots be created
                                                      y_threshold = NDSI_threshold_vector) #Which threshold values for 'y' should be calculated

          #Turn parallel processing off and run sequentially again after this point
          stopCluster(cl)
          registerDoSEQ()  
          
          #Store dates of snowmelt per SubArea
          df_SubArea_NDSI <- results[[1]]
          df_SubArea_NDSI <- as.data.frame(do.call(rbind, df_SubArea_NDSI))
          colnames(df_SubArea_NDSI)[colnames(df_SubArea_NDSI)=="pixel_ID"] <- "SubArea"
          colnames(df_SubArea_NDSI)[colnames(df_SubArea_NDSI)=="x_threshold"] <- "doy"
          colnames(df_SubArea_NDSI)[colnames(df_SubArea_NDSI)=="y_threshold"] <- "NDSI_threshold"
         
          #Change column subArea to a factor:
          df_SubArea_NDSI$SubArea <- as.factor(as.character(df_SubArea_NDSI$SubArea))

          #Save dates of snowmelt per subarea per NDSI threshold as a .csv file
          write.csv(df_SubArea_NDSI, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Snowmelt_NDSI.csv"), row.names = FALSE)

          #Add dates of snowmelt to the plot 'p_SubArea_NDSI_grid'
          p_SubArea_NDSI_Snowmelt_grid <- p_SubArea_NDSI_grid +
            geom_point(data=df_SubArea_NDSI[!is.na(df_SubArea_NDSI$doy),], aes(x=doy, y=NDSI_threshold), col="red", size=3)

          ggsave(plot=p_SubArea_NDSI_Snowmelt_grid, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Plot_Snowmelt_NDSI_grid.pdf"), width=16, height = 14)

########################################################################################################################################################################################
              
  #(C): NDVI - Fit a Generalized Additive Model (GAM) through the NDVI data within each sub-area

      #(C.1) Create an empty dataframe
        df_SubAreas_NDVI_GAM <- data.frame(NDVI=numeric(),
                                           Date=factor(),
                                           doy=numeric(),
                                           SubArea=factor(),
                                           outliers=logical())

        df_SubAreas_NDVI_GAM_predictions <- data.frame(SubArea=character(),
                                                       doy=numeric(),
                                                       NDVI_gam_predict=numeric(),
                                                       stringsAsFactors=FALSE)

      #(C.2) Loop through all SubAreas and fit a separate gam (with sequential outlier removal) to the SubArea-specific NDVI data
        for(i in unique(df_SubAreas_BandValues$SubArea)){

          #For debugging
          #i=unique(df_SubAreas_BandValues$SubArea)[1]

          #Run code for the first NDSI_threshold (as NDVI-data is collected independently of thresholds)
          j=unique(df_SubAreas_BandValues$NDSI_threshold)[1]

          #Select SubArea-specific subset of data:
          df_SubArea_NDVI_GAM_new <- df_SubAreas_BandValues[df_SubAreas_BandValues$SubArea==i &
                                                            df_SubAreas_BandValues$NDSI_threshold==j &
                                                            !is.na(df_SubAreas_BandValues$NDVI),
                                                            c("NDVI", "Date", "doy", "SubArea")]

          #Fit a gam through the SubArea-specific NDVI ~ doy data and emply sequential outlier removal
          df_SubArea_NDVI_GAM_new <- f_gam_SeqRemOutliers(data=df_SubArea_NDVI_GAM_new, y="NDVI", x="doy", outlier_removal=outlier_removal,
                                                           outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                           default_k=gam_k_outlier)

          #Sort df_SubArea_NDVI_GAM_new by doy:
          df_SubArea_NDVI_GAM_new <- df_SubArea_NDVI_GAM_new[order(df_SubArea_NDVI_GAM_new$doy),]

          #Bind the SubArea-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
           df_SubAreas_NDVI_GAM <- rbind(df_SubAreas_NDVI_GAM, df_SubArea_NDVI_GAM_new)

          #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

            #Refit GAM through data
             index <- which(df_SubArea_NDVI_GAM_new$outliers==FALSE)
             mod_gam <- with(df_SubArea_NDVI_GAM_new[index,], mgcv::gam(NDVI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

            #Use gam to make predictions on a more detailed (1-day) day of year interval
             df_SubArea_NDVI_GAM_predictions_new <- data.frame(SubArea=i, doy=seq(min(df_SubArea_NDVI_GAM_new$doy), max(df_SubArea_NDVI_GAM_new$doy), 1))
             df_SubArea_NDVI_GAM_predictions_new$NDVI_gam_predict <- stats::predict(mod_gam, newdata=df_SubArea_NDVI_GAM_predictions_new, type="response")
             df_SubArea_NDVI_GAM_predictions_new <- df_SubArea_NDVI_GAM_predictions_new[order(df_SubArea_NDVI_GAM_predictions_new$doy),]

            #Add predictions to df_SubAreas_NDVI_GAM_predictions dataframe:
             df_SubAreas_NDVI_GAM_predictions <- rbind(df_SubAreas_NDVI_GAM_predictions, df_SubArea_NDVI_GAM_predictions_new)

          }

        #Change subArea column to factor
        df_SubAreas_NDVI_GAM$SubArea <- as.factor(as.character(df_SubAreas_NDVI_GAM$SubArea))
        df_SubAreas_NDVI_GAM_predictions$SubArea <- as.factor(as.character(df_SubAreas_NDVI_GAM_predictions$SubArea))

        #Save dataframe with GAM fits for NDVI
        #write.csv(df_SubAreas_NDVI_GAM, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_NDVI_GAM.csv"), row.names = FALSE)
        write.csv(df_SubAreas_NDVI_GAM_predictions, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_NDVI_Predictions_GAM.csv"), row.names = FALSE)

      #(C.3) Plot the raw NDVI datapoints and gam predictions for each SubArea:

         #Plot NDVI and model predictions for all SubAreas in a single plot
          p_SubArea_NDVI = ggplot()+
            geom_point(data=df_SubAreas_NDVI_GAM, aes(x=doy, y=NDVI, fill=SubArea, col=SubArea))+
            geom_line(data=df_SubAreas_NDVI_GAM_predictions, aes(x=doy, y=NDVI_gam_predict, col=SubArea)) +
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("NDVI") +
            theme_tom()
            
         #Plot NDVI and model predictions in a separate plot per SubArea
          p_SubArea_NDVI_grid = ggplot()+
            geom_point(data=df_SubAreas_NDVI_GAM[df_SubAreas_NDVI_GAM$outliers==FALSE,], aes(x=doy, y=NDVI))+
            geom_point(data=df_SubAreas_NDVI_GAM[df_SubAreas_NDVI_GAM$outliers==TRUE,], aes(x=doy, y=NDVI), col="black", pch=16, alpha=0.2)+
            geom_line(data=df_SubAreas_NDVI_GAM_predictions, aes(x=doy, y=NDVI_gam_predict), col="#1620de", lwd=1.25)+
            geom_vline(xintercept = 150, colour="grey", lty=2)+
            facet_wrap(~SubArea, ncol=ceiling(length(unique(df_SubAreas_NDVI_GAM$SubArea))^0.5))+
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("NDVI") +
            theme_tom()

            ggsave(plot=p_SubArea_NDVI_grid, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Plot_Snowmelt_NDVI_grid.pdf"), width=16, height = 14)
            
########################################################################################################################################################################################

  #(D): NDMI - Fit a Generalized Additive Model (GAM) through the NDMI data within each sub-area

      #(D.1) Create an empty dataframe
        df_SubAreas_NDMI_GAM <- data.frame(NDMI=numeric(),
                                           Date=factor(),
                                           doy=numeric(),
                                           SubArea=factor(),
                                           outliers=logical())

        df_SubAreas_NDMI_GAM_predictions <- data.frame(SubArea=character(),
                                                       doy=numeric(),
                                                       NDMI_gam_predict=numeric(),
                                                       stringsAsFactors=FALSE)

      #(D.2) Loop through all SubAreas and fit a separate gam (with sequential outlier removal) to the SubArea-specific NDMI data
        for(i in unique(df_SubAreas_BandValues$SubArea)){

          #For debugging
          #i=unique(df_SubAreas_BandValues$SubArea)[1]

          #Run code for the first NDSI_threshold (as NDMI-data is collected independently of thresholds)
          j=unique(df_SubAreas_BandValues$NDSI_threshold)[1]

          #Select SubArea-specific subset of data:
          df_SubArea_NDMI_GAM_new <- df_SubAreas_BandValues[df_SubAreas_BandValues$SubArea==i &
                                                            df_SubAreas_BandValues$NDSI_threshold==j &
                                                            !is.na(df_SubAreas_BandValues$NDMI),
                                                            c("NDMI", "Date", "doy", "SubArea")]

          #Fit a gam through the SubArea-specific NDMI ~ doy data and emply sequential outlier removal
          df_SubArea_NDMI_GAM_new <- f_gam_SeqRemOutliers(data=df_SubArea_NDMI_GAM_new, y="NDMI", x="doy", outlier_removal=outlier_removal,
                                                           outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                           default_k=gam_k_outlier)

          #Sort df_SubArea_NDMI_GAM_new by doy:
          df_SubArea_NDMI_GAM_new <- df_SubArea_NDMI_GAM_new[order(df_SubArea_NDMI_GAM_new$doy),]

          #Bind the SubArea-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
           df_SubAreas_NDMI_GAM <- rbind(df_SubAreas_NDMI_GAM, df_SubArea_NDMI_GAM_new)

          #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

            #Refit GAM through data
             index <- which(df_SubArea_NDMI_GAM_new$outliers==FALSE)
             mod_gam <- with(df_SubArea_NDMI_GAM_new[index,], mgcv::gam(NDMI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

            #Use gam to make predictions on a more detailed (1-day) day of year interval
             df_SubArea_NDMI_GAM_predictions_new <- data.frame(SubArea=i, doy=seq(min(df_SubArea_NDMI_GAM_new$doy), max(df_SubArea_NDMI_GAM_new$doy), 1))
             df_SubArea_NDMI_GAM_predictions_new$NDMI_gam_predict <- stats::predict(mod_gam, newdata=df_SubArea_NDMI_GAM_predictions_new, type="response")
             df_SubArea_NDMI_GAM_predictions_new <- df_SubArea_NDMI_GAM_predictions_new[order(df_SubArea_NDMI_GAM_predictions_new$doy),]

            #Add predictions to df_SubAreas_NDMI_GAM_predictions dataframe:
             df_SubAreas_NDMI_GAM_predictions <- rbind(df_SubAreas_NDMI_GAM_predictions, df_SubArea_NDMI_GAM_predictions_new)

          }

        #Change subArea column to factor
        df_SubAreas_NDMI_GAM$SubArea <- as.factor(as.character(df_SubAreas_NDMI_GAM$SubArea))
        df_SubAreas_NDMI_GAM_predictions$SubArea <- as.factor(as.character(df_SubAreas_NDMI_GAM_predictions$SubArea))

        #Save dataframe with GAM fits for NDMI
        #write.csv(df_SubAreas_NDMI_GAM, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_NDMI_GAM.csv"), row.names = FALSE)
        write.csv(df_SubAreas_NDMI_GAM_predictions, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Data_NDMI_Predictions_GAM.csv"), row.names = FALSE)

      #(D.3) Plot the raw NDMI datapoints and gam predictions for each SubArea:

         #Plot NDMI and model predictions for all SubAreas in a single plot
          p_SubArea_NDMI = ggplot()+
            geom_point(data=df_SubAreas_NDMI_GAM, aes(x=doy, y=NDMI, fill=SubArea, col=SubArea))+
            geom_line(data=df_SubAreas_NDMI_GAM_predictions, aes(x=doy, y=NDMI_gam_predict, col=SubArea)) +
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("NDMI") +
            theme_tom()
            
         #Plot NDMI and model predictions in a separate plot per SubArea
          p_SubArea_NDMI_grid = ggplot()+
            geom_point(data=df_SubAreas_NDMI_GAM[df_SubAreas_NDMI_GAM$outliers==FALSE,], aes(x=doy, y=NDMI))+
            geom_point(data=df_SubAreas_NDMI_GAM[df_SubAreas_NDMI_GAM$outliers==TRUE,], aes(x=doy, y=NDMI), col="black", pch=16, alpha=0.2)+
            geom_line(data=df_SubAreas_NDMI_GAM_predictions, aes(x=doy, y=NDMI_gam_predict), col="#1620de", lwd=1.25)+
            geom_vline(xintercept = 150, colour="grey", lty=2)+
            facet_wrap(~SubArea, ncol=ceiling(length(unique(df_SubAreas_NDMI_GAM$SubArea))^0.5))+
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("NDMI") +
            theme_tom()

            ggsave(plot=p_SubArea_NDMI_grid, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Plot_Snowmelt_NDMI_grid.pdf"), width=16, height = 14)
           
########################################################################################################################################################################################
                
    #(E) Plot NDSI, NDVI and NDMI together in a plot per SubArea
    
      #Plot:   
       p_SubArea_BANDS_grid <- ggplot()+ 
         geom_point(data=df_SubAreas_NDSI_GAM[df_SubAreas_NDSI_GAM$outliers==FALSE,], aes(x=doy, y=NDSI))+
         geom_point(data=df_SubAreas_NDSI_GAM[df_SubAreas_NDSI_GAM$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
         geom_point(data=df_SubAreas_NDVI_GAM[df_SubAreas_NDVI_GAM$outliers==FALSE,], aes(x=doy, y=NDVI), col="black")+
         geom_point(data=df_SubAreas_NDVI_GAM[df_SubAreas_NDVI_GAM$outliers==TRUE,], aes(x=doy, y=NDVI), col="black", pch=16, alpha=0.2)+
         geom_point(data=df_SubAreas_NDMI_GAM[df_SubAreas_NDMI_GAM$outliers==FALSE,], aes(x=doy, y=NDMI), col="black")+
         geom_point(data=df_SubAreas_NDMI_GAM[df_SubAreas_NDMI_GAM$outliers==TRUE,], aes(x=doy, y=NDMI), col="black", pch=16, alpha=0.2)+
         geom_line(data=df_SubAreas_NDSI_GAM_predictions, aes(x=doy, y=NDSI_gam_predict), col="#1620de", lwd=1.25)+
         geom_line(data=df_SubAreas_NDVI_GAM_predictions, aes(x=doy, y=NDVI_gam_predict), col="#08a31a", lwd=1.25)+
         geom_line(data=df_SubAreas_NDMI_GAM_predictions, aes(x=doy, y=NDMI_gam_predict), col="#16acde", lwd=1.25)+
         geom_vline(xintercept = 150, colour="grey", lty=2)+
         facet_wrap(~SubArea, ncol=ceiling(length(unique(df_SubAreas_NDMI_GAM$SubArea))^0.5))+
         xlab("Day of year") + 
         ylab("Normalized Difference Band Index") +
         theme_tom()
       
         ggsave(plot=p_SubArea_BANDS_grid, paste0(here(), "/Output/S2/06_Shapefile_SubAreas_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_SubAreas_Plot_AllBands_grid.pdf"), width=14, height = 12)
         
             
##################################################################################################################################       
##################################################################################################################################       
           
        