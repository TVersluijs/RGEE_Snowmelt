##################################################################################################################################

#Extract MODIS satellite data and Calculate the date of snowmelt for every 500mx500m pixel in an area of interest (shapefile)

#Copyright Tom Versluijs 2023-02-08. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R". 
#Note that a GoogleDrive is required.

##################################################################################################################################

#I: Setup workspace

##################################################################################################################################

      #(0): Clear workspace and set python environment
       rm(list=ls())
       library(here)
       rgee_environment_dir <- readRDS(paste0(here::here(), "/Input/rgee_environment_dir.rds"))
       reticulate::use_python(rgee_environment_dir, required=T)
       reticulate::py_config()

      #(1): Load packages
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
       source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "MODIS_AuxiliaryFunctions")
       sapply(source_files, source, chdir = TRUE) ; rm(source_files)
       
      #(4): Initialize earth engine and google drive
       rgee::ee_Initialize(user="tom.versluijs@gmail.com", drive = T)
       
       
##################################################################################################################################
       
#II: Specify parameters of interest
       
##################################################################################################################################       
       
 #(5): Specify parameters used in the analysis

   #(a): MODIS satellite

     #Specify resolution of images in meters
     resolution=500 # native resolution of satellite image (MODIS=500m)

   #(b) Area of interest

     #Specify name of study area
     area_name <- "ZAC"

     #Name of Shapefile located in Input folder (specific area of interest)
     #Follow the guide "Manual_CreateShapefilePolygons.docx" when creating this shapefile.
     shapefile <- "ZAC_Outline_EPSG4326.shp"

     #Coordinate reference system used for calculations
     #EPSG:4326 is recommended for areas spanning multiple UTM zones, but increased computation time (i.e. spherical coordinate system).
     #EPSG:326XX is results in reduced computation time for areas located within a single UTM zone (i.e. planar coordinate system).
     crs <- "EPSG:4326"

   #(c) Dates

     #Specify the year of interest:
     year_ID <- "2022"
     
     #Date range for all images considered in analysis
     start_date <- paste0(year_ID, "-03-15") #choose date (well) before the first snowmelt occurs within the study site
     end_date <- paste0(year_ID, "-09-15") #choose date (well) after last date of tracking

   #(d) Snow detection
     
     #NDSI threshold above which a pixel is perceived as snow
     NDSI_threshold=0.42

   #(e): Cloud masking

     #Should clouds be masked from the analysis (default=TRUE).
     mask_clouds=TRUE
     
     #Define which MODIS cloud masking algorithm has been used ("PGE11", "MOD35", or "Combined)
     MODIS_cloud_masking_algorithm = 'PGE11' #default is PGE11

     #Maximum fraction of cloud cover allowed in each image (should be set to 1.0 for pixel-level analysis)
     max_cloud_fraction=1.0
     
     #Specify the resolution of the QABand   
     resolution_qaband=1000 #resolution of QA_band (MODIS=1000m)   
     
   #(f): Water masking

     #Should permanent waterbodies be masked from the analysis (default=TRUE).
     mask_water=TRUE
     
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
         
      #(5) Create a unique dataID and asset folder for storing the generated datafiles
       
        #Create a unique data_ID
         data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
         data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)
          
        #Create a unique Asset folder (delete this folder if already present) 
         path_asset <- paste0(ee_get_assethome(), "/", data_ID)
         #tryCatch(ee_manage_assetlist(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
         tryCatch(ee_manage_delete(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
         ee_manage_create(path_asset=path_asset, asset_type="Folder")
         ee_manage_assetlist()
      
         
##################################################################################################################################
         
#IV: Read and display the unfiltered data
         
##################################################################################################################################

      #(7) Read study area shapefile as a feature collection and clip image collection to area of interest
         
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
           crs= crs,
           scale= 10,
           maxPixels= 1E13)
         paste0('Size of study area calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')

      #(8) Extract MODIS satellite Surface Reflectance images for specified daterange and general area of interest
       
       #Note that specifying a region of interest using a polygon does not function properly for the MODIS data. 
       #Instead, providing an initial point and then afterwards clipping the image to a desired area of interest works better.
       
        #Specify starting and ending date  
         start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
         end_date_doy <- as.numeric(strftime(end_date, format = "%j"))
       
        #Extract Sentinel-2 Surface Reflectance satellite data
         MODIS_col <- ee$ImageCollection('MODIS/006/MOD09GA')
         MODIS_col <- MODIS_col$
           filterBounds(coordinates_point)$
           filterDate(start_date, end_date)
         
        #Clip all images to area of interest (aoi). Clipping to shapefile will occur at step 16B of this script
         MODIS_col <- MODIS_col$map(function(img){return(img$clip(aoi))})
        
        #Plot a single clipped image: 
         
          #Select a single image for initial plot, clipped to aoi_Shapefile:
           image <- MODIS_col$filterDate(paste0(year_ID, "-06-10"), end_date)$first()$clipToCollection(aoi_Shapefile)
  
          #Add Normalized difference snow index to this image:
           ndsi_MODIS <-  image$expression('(B4 - B6) / (B4 + B6)',
                                        list('B4'=image$select('sur_refl_b04'),
                                             'B6'=image$select('sur_refl_b06')))
        
          #Add Normalized difference vegetation index to this image (CHECKED):
           ndvi_MODIS <- image$expression('(B2 - B1) / (B2 + B1)',
                                       list('B2'=image$select('sur_refl_b02'),
                                            'B1'=image$select('sur_refl_b01')))
           
          #Add Normalized difference Moisture index to this image:
           ndmi_MODIS <- image$expression('(B2 - B6) / (B2 + B6)',
                                       list('B2'=image$select('sur_refl_b02'),
                                            'B6'=image$select('sur_refl_b06'))) 
           
          #Add Normalized difference water index to this image:
           ndwi_MODIS <- image$expression('(B4 - B2) / (B4 + B2)',
                                       list('B4'=image$select('sur_refl_b04'),
                                            'B2'=image$select('sur_refl_b02')))
    
          #Plot all image Bands
           Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
           Map$addLayer(image,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
           Map$addLayer(ndsi_MODIS,list(min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
           Map$addLayer(ndvi_MODIS,list(min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
           Map$addLayer(ndwi_MODIS,list(min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')
       
      #(9) Add NDSI, NDVI, NDMI and NDWI to the clipped image collection
        
        #Map normalized Difference functions over image collection
        MODIS_col <- MODIS_col$
          map(getNDSI)$
          map(getNDVI)$
          map(getNDMI)$
          map(getNDWI)
      
      #(10) Create a timeseries gif of RGB images for the aoi_Shapefile shapefile (for debugging)
        
        # #Check number of images in collection
        #  MODIS_col$size()$getInfo()

        # #Create a timelapse video of RGB band
        #  videoArgs <- list(dimensions=380, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), 
        #                    min=0, max=12000, gamma=c(1.9, 1.7, 1.7))
        #  browseURL(MODIS_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$getVideoThumbURL(videoArgs)) 
        
        # #Create a timelapse video of NDSI band  
        #  palette=c('black', '0dffff', '0524ff', 'ffffff')
        #  visFun_NDSI <-  function(img) {
        #    return(img$visualize(bands='NDSI', min=-1, max=1.5, palette=palette)$copyProperties(img, img$propertyNames()))
        #    }
        #  MODIS_snow_RGB <- MODIS_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$map(visFun_NDSI)
        #  videoArgs <- list(dimensions=380, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
        #  browseURL(MODIS_snow_RGB$getVideoThumbURL(videoArgs)) 
      
         #Note that MODIS_col is only clipped by 'aoi' and not yet by aoi_Shapefile.
         
##################################################################################################################################
        
#V: Filter and mask clouds within the image collection
        
##################################################################################################################################

   #There are many images that are almost fully covered in dense clouds. We thus first filter and mask clouds in the image collection.
   #We delete all images with a cloud cover fraction >= max_cloud_fraction. In the remaining images we mask all pixels that are 
   #covered by clouds. 
   if(mask_clouds==TRUE){      
         
     #print message   
     print("Cloud masking = TRUE")
     
     #The MODIS quality band 'state_1km' provides several cloud flag algorithms (See Wilson 2014, MODIS user guide). The cloud flags from 
     #the 'MOD35'  algorithm are contained in bits 0-1 and have four categories: (0) clear, (1) cloudy, (2) mixed, (3) not set, assumed clear.
     #The cloud flags from the MODIS Land Team 'PGE11' internal cloud masking algorithm are contained in bit 10 and have two categories: 
     #(0) no cloud, (1) cloud. Cirrus cloud flags are stored in bits 8-9 and have four categories: (0) none, (1) small, (2) average, 
     #(3) high. The MOD35 cloud mask uses 22 MODIS bands, ecosystem/landcover type and other environmental data to identify clouds
     #(Wilson et al 2014), while the PGE11 internal cloud mask detects clouds based on two reflective testes and a thermal test. Usually
     #clouds are filtered if one or both of these two algorithms detects clouds (Wilson et al 2014).

      #(11) Add cloud information within the aoi to the image collection: 

        #Add the cloud fraction within the area of interest to each separate image by mapping the cloud functions over the image collection
        MODIS_col <- MODIS_col$
           #Determine which pixels are either opague or cirrus clouds
           map(computeClouds)$
           #Add the fraction of cloud-covered pixels within the area of interest as image property
           map(AddCloudFraction)$
           #Add NULL to those images in which cloudfraction could not be calculated
           map(AddNULLCloudFraction)$
           #Add date and time characteristics to each image
           map(add_Date)
         
        # #Check if Cloud information has been added to the properties of each image (for debugging)
        # MODIS_col$first()$propertyNames()$getInfo()
      
        # #Extract cloud fraction of all images in image collection for the area of interest (for debugging)
        # Combined_clouds_Fraction <- unlist(MODIS_col$aggregate_array('Combined_clouds_Fraction')$getInfo())
        # PGE11_clouds_Fraction <- unlist(MODIS_col$aggregate_array('PGE11_clouds_Fraction')$getInfo())
        # MOD35_clouds_Fraction <- unlist(MODIS_col$aggregate_array('MOD35_clouds_Fraction')$getInfo())
        # doy <- unlist(MODIS_col$aggregate_array('doy')$getInfo())
        
        # #Replace -9999 values by NA (for debugging)
        # Combined_clouds_Fraction[Combined_clouds_Fraction < -9000] <- NA
        # PGE11_clouds_Fraction[PGE11_clouds_Fraction < -9000] <- NA
        # MOD35_clouds_Fraction[MOD35_clouds_Fraction < -9000] <- NA
        # doy[doy < -9000] <- NA
        
        # #Combine all cloud fraction measures in a single dataframe (for debugging)
        # clouds_MODIS <- data.frame(doy=doy,
        #                            Combined_clouds_Fraction=Combined_clouds_Fraction,
        #                            PGE11_clouds_Fraction=PGE11_clouds_Fraction,
        #                            MOD35_clouds_Fraction=MOD35_clouds_Fraction)
        
        # #Assess relationship between both cloud algorithms (for debugging)
        # ggplot()+
        #   geom_line(aes(x=doy, y=MOD35_clouds_Fraction), col="red", lwd=1)+
        #   geom_line(aes(x=doy, y=PGE11_clouds_Fraction), col="black", lwd=1)+
        #   #geom_line(aes(x=doy, y=Combined_clouds_Fraction), col="blue", lwd=1)+
        #   theme_tom()
        
        # #Correlation between MOD35 and PGE11 (for debugging)
        # ggplot()+
        #   geom_point(aes(x=MOD35_clouds_Fraction, y=PGE11_clouds_Fraction), col="black")+
        #   geom_abline(slope=1, intercept=0, lty=2)+
        #   theme_tom()
        
        #Especially during the period before snowmelt, MOD35 seems to always detect at least 25% cloudcover, while PGE11 detects many 
        #cloud free days. The former cannot be correct. Thus, MOD35 results in many false positives when snow is still present. PGE11
        #is more conservative in detecting clouds when snow is present.
  
        # #(for debugging) Manually compare MODIS RGB images to the different cloud detection algorithms during the most critical period for snowmelt.
        # #Thus, make sure to manually alter the date below to check how MODIS cloudcover algorithms correspond to the underlying
        # #RGB images at different points in time. Note that white colours in each mask indicate clouds.
        #  image <- MODIS_col$
        #     filterDate(paste0(year_ID, "-06-20"), end_date)$
        #     first()$
        #     clipToCollection(aoi_Shapefile)$
        #     select('MOD35_clouds', 'PGE11_clouds', 'Combined_clouds', "sur_refl_b01", "sur_refl_b04", "sur_refl_b03")      
        #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
        #  Map$addLayer(image, list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
        #  Map$addLayer(image, list(bands='MOD35_clouds', min=0, max=1, opacity=1), 'MOD35_clouds')+
        #  Map$addLayer(image, list(bands='PGE11_clouds', min=0, max=1, opacity=1), 'PGE11_clouds')+
        #  Map$addLayer(image, list(bands='Combined_clouds', min=0, max=1, opacity=1), 'MOD35 and PGE11 combined')

        #These images also indicate that the PGE11 algorithm does a more conservative job that better matches clouds in the RGB images than the MOD35 algorithm.

      #(12): Exclude all images from the image collection that have a CloudFraction value >= max_cloud_fraction, and mask all cloud pixels within the remaining images
        
        #In the analysis for the whole study area we exclude images with more than 20% snowcover within our area of 
        #interest. However, now we're interested in individual pixels. So even when 90% of the area of interest is 
        #snowcovered, the pixels we're interested in might be located in the 10% cloud free area. We thus only filter
        #images with a 100% cloudcover (cloudFraction=1.0). For all remaining images we apply a cloud mask to mask 
        #out cloudy pixels within each image.
       
        #Cloud function III: Define a cloudmask function for the selected cloud masking algorithm
         cloud_algorithm <- paste0(MODIS_cloud_masking_algorithm, "_clouds")
         
        #Apply a manually constructed cloudmask to each image in the collection:
         clouds_fraction <- paste0(MODIS_cloud_masking_algorithm, "_clouds_Fraction")
         MODIS_clouds_filtered <- MODIS_col$
           #Filter all images with 100% cloudcover (PGE11_clouds_Fraction=1.0)
           filter(ee$Filter$lt(clouds_fraction, max_cloud_fraction))$
           #Apply cloudmask for individual pixels
           map(AddCloudMask)
         
        # #Create timelapse video of the cloud filtered/masked RGB images (for debugging)
        #  videoArgs <- list(dimensions=380, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
        #  browseURL(MODIS_clouds_filtered$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$getVideoThumbURL(videoArgs)) 
     
        # #Create a timelapse video of the cloud filtered/masked NDSI band (for debugging) 
        #  palette=c('black', '0dffff', '0524ff', 'ffffff')
        #  visFun_NDSI <-  function(img) {
        #    return(img$visualize(bands='NDSI', min=-1, max=1.5, palette=palette)$
        #             copyProperties(img, img$propertyNames()))}
        #  MODIS_snow_masked_RGB <- MODIS_clouds_filtered$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$map(visFun_NDSI)
        #  videoArgs <- list(dimensions=380, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
        #  browseURL(MODIS_snow_masked_RGB$getVideoThumbURL(videoArgs)) 

       }
   if(mask_clouds==FALSE){
     
     #print message   
     print("Cloud masking = FALSE")
     
     #Add Date and Time to MODIS collection
     MODIS_clouds_filtered <- MODIS_col$map(add_Date)

   }       
         
   #Store default MODIS image projection 
   modisProjection <- MODIS_clouds_filtered$first()$select("NDSI")$projection()
   #modisProjection$getInfo()
         
   #Note that MODIS_clouds_filtered is not yet clipped by aoi_Shapefile (only by aoi)!


##################################################################################################################################
         
#VI: Mask permanent waterbodies (ponds, lakes, rivers, sea) within the image collection
         
##################################################################################################################################            
      
   #(13): Mask permanent waterbodies if mask_water==TRUE
   if(mask_water==TRUE){   
         
      #(A): print message   
        print("Water masking = TRUE")
      
      #(B): Extract permanent waterbodies from the ESA WorldCover map
    
        #Load auxilliary function
        compute_Water_MODIS=compute_Water_MODIS #sourced
        
        #Extract pixels corresponding to permanent waterbodies 
        water_mask <- compute_Water_MODIS()
        
      #(C): Apply the water masking function to each image in the collection:
      
        # #Plot permanent waterbodies(for debugging)
        # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
        # Map$addLayer(MODIS_clouds_filtered$filterDate(paste0(year_ID, "-08-04"), end_date)$first(), list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
        # Map$addLayer(water_mask,list(min=0, max = 1, palette = c('ffffff', '000000')), 'Water_mask')
      
        #Load auxilliary function to mask water-pixels
        Add_WaterMask_MODIS=Add_WaterMask_MODIS #sourced
      
        #Apply the final watermask:
        MODIS_clouds_filtered <- MODIS_clouds_filtered$map(Add_WaterMask_MODIS)
      
        # #Create a timeseries GIF of RGB images of the water and cloud filtered image collection (for debugging)
        # videoArgs <- list(dimensions=380, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
        # browseURL(MODIS_clouds_filtered$getVideoThumbURL(videoArgs))
      
      }
   if(mask_water==FALSE){
      
      #(A): print message   
      print("Water masking = FALSE")
      
    }
     
   #Note that MODIS_clouds_filtered is not yet clipped by aoi_Shapefile (only by aoi)! 
         
##################################################################################################################################
        
#VII: Calculate the date of snowmelt for every pixel within the study area by fitting a GAM through the NDSI data
        
##################################################################################################################################            
         
      #(14): Calculate the date of snowmelt (NDSI <= NDSI_threshold) for every pixel within the study area (bounding box!) 
       
          #(A): Transform each image to a feature Collection of NDSI values for all pixels within 'aoi'
         
           #Create an iteration function that we will use to iterate through all images of the image collection. For each image, the
           #value of the NDSI band is extracted for each pixel of the image. The resulting NDSI values (+lat/lon, datetime of the image)
           #of each pixel are stored as feature properties in a feature collection (i.e. where each pixel has a separate feature with
           #properties). This results in a feature collection of all pixels for the current image. This feature collection is then
           #appended to a list of feature collections from previous image iterations.

           #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
            FC_initial <- ee$FeatureCollection(ee$List(list()))

           #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
           #(in this case the current iteration image) and the second element takes the output value from the iteration that preceeded it. The latter
           #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with
           #has been defined. Note that inside $map() functions all processing has to be done in the language of the server (javascript
           #Api of google earth engine). Thus, inside $map() functions, client-side functions such as $getInfo() or print(), paste0()
           #cannot be used.
            Extract_BandValuesAtPixels = Extract_BandValuesAtPixels #note that region=aoi

           #Iterate over the ImageCollection (output is a large feature collection)
            FC_merged <- ee$FeatureCollection(MODIS_clouds_filtered$select("NDSI")$iterate(Extract_BandValuesAtPixels, FC_initial))
            #FC_merged$first()$getInfo() #for debugging
           
            #Note that at this point all pixels within 'aoi' are included (not only those in 'aoi_Shapefile'!)
            
         #(B): Transform feature collection to a dataframe:

           #We use ee_table_to_drive() to prevent memory limits
            a=Sys.time()
            task_vector1 <- ee_table_to_drive(
              collection = FC_merged,
              description = paste0(data_ID, "_Pixel_NDSI_Snowmelt"),
              fileFormat = "CSV",
              selectors = c('NDSI', 'Date', 'doy', 'lat', 'lon')
              )

            task_vector1$start()
            print("Transform each image to a feature Collection of NDSI values for all pixels:")
            ee_monitoring(task_vector1, max_attempts = 1000000)

            exported_stats <- ee_drive_to_local(task = task_vector1, dsn=paste0("Output/", data_ID, "_Pixel_NDSI_Snowmelt"))
            df_pixel_ndsi <- read.csv(exported_stats)
            b=Sys.time()
            print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))

           # #Load dataframe (takes ca 2 minutes):
           #  df_pixel_ndsi <- read.csv(paste0(data_ID, "_Pixel_NDSI_Snowmelt.csv"))

           #Make sure each latitude/longitude combination gets its own pixel_ID (takes ca 1 minute):
            df_pixel_ndsi$pixel_ID <- paste0(format(df_pixel_ndsi$lat, nsmall = 5), "_", format(df_pixel_ndsi$lon, nsmall = 5))

         #(C): Calculate the date of snowmelt for each pixel in the dataframe
            
           #Loop through all pixel_ID's, select dataframe for that pixel containing NDSI values
           #as measured in all images in the image collection, fit gam through the NDSI ~ doy data, 
           #determine date when NDSI<NDSI_threshold and store this date of snowmelt together with 
           #the pixel_ID in a new dataframe.

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

           #Split dataset into chuncks of pixels
            chunk_size = 2000
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
           #year the predicted NDSI value of this GAM changes from above outlier_threshold to below (direction="down") or
           #from below outlier_threshold to above (direction="up"). The code employs parallel processing using foreach
           #and %dopar% on four local computer cores.
            f_detect_threshold_date_parallel_2 <- f_detect_threshold_date_parallel_2 #sourced

           #Run the 'f_detect_threshold_date_parallel_2' function over all data subsets, combine the results and save the resulting dataframe and plots
            print("Calculate the date of snowmelt for each pixel within aoi:")
            results <- lapply(1:length(pixelIDs_split), FUN=f_detect_threshold_date_parallel_2, 
                              df=df_pixel_ndsi, pixel_ID_column="pixel_ID", y="NDSI", 
                              x="doy", direction="down", y_threshold=NDSI_threshold)
              
            #Clean up the cluster after finishing the parallel runs
            stopCluster(cl)
            
            #Turn parallel processing off and run sequentially again after this point
            registerDoSEQ()
            
           #Store date of snowmelt per pixel within aoi (bounding box!) as a dataframe
            df_pixel_snowmelt <- lapply(results, "[[", 1)
            df_pixel_snowmelt <- as.data.frame(do.call(rbind, do.call(c, df_pixel_snowmelt)))
            colnames(df_pixel_snowmelt)[colnames(df_pixel_snowmelt)=="x_threshold"] <- "doy_snowmelt"
            write.csv(df_pixel_snowmelt, file=paste0(here(), "/Output/", data_ID, "_Pixel_Snowmelt_aoi.csv"), quote = FALSE, row.names=FALSE)
            ##Read dataframe
            #df_pixel_snowmelt <- read.csv(file=paste0(here(), "/Output/", data_ID, "_Pixel_Snowmelt_aoi.csv"), header=TRUE)  
            
           #Store plots
            plot_pixel_snowmelt <- lapply(results, "[[", 2)
            plots_per_page = 25
            plot_pixel_snowmelt <- lapply(plot_pixel_snowmelt, function(x){split(x, ceiling(seq_along(plot_pixel_snowmelt[[1]])/plots_per_page))})
            plot_pixel_snowmelt <- unname(unlist(plot_pixel_snowmelt, recursive = F))
            pdf(paste0("Output/", data_ID, "_Pixel_NDSI_Snowmelt_aoi.pdf"), width=20, height=16, onefile = TRUE)
            for (i in seq(length(plot_pixel_snowmelt))) { do.call("grid.arrange", plot_pixel_snowmelt[[i]]) }
            dev.off()
        
           #The dataframe df_pixel_snowmelt now contains the date of snowmelt for each individual pixel_ID (514584 10mx10m pixels)
           #withing the aoi (i.e. the defined bounding box). To be able to plot these data we transform this dataframe to a feature 
           #collection and then transform this feature collection to an image.
            
       #(15): Transform df_pixel_snowmelt to a feature collection with random geometry
            
           #Generate some random longitude and latitude values as this is required for an sf object:
            df_pixel_snowmelt$lon <- runif(nrow(df_pixel_snowmelt), 0, 74)
            df_pixel_snowmelt$lat <- runif(nrow(df_pixel_snowmelt), -20, 20)
            
           #Prevent NA in doy_snowmelt column (earthengine cannot deal with NA)  
            df_pixel_snowmelt$doy_snowmelt[is.na(df_pixel_snowmelt$doy_snowmelt)] <- -9999
            
           #memory limits prevent us from processing the whole dataframe at once. Therefore we split it up into subsets with 'chunk_size' rows
            chunk_size = 2000
            rowIDs <- 1:nrow(df_pixel_snowmelt)
            rowIDs_split <- split(rowIDs, ceiling(seq_along(rowIDs)/chunk_size))
          
           #Iterate through all dataframe subsets, convert each to a feature collection, and append each to the feature collection FC_initial.
            FC_initial <- ee$FeatureCollection(ee$List(list()))
            for(i in 1:length(rowIDs_split)){ 
              
              #Select subset of rows
               rowID_min <- min(rowIDs_split[[i]])
               rowID_max <- max(rowIDs_split[[i]])
              
              #Select subset of dataframe:
               df_tmp <- df_pixel_snowmelt[rowID_min:rowID_max,]
              
              #Change subset-dataframe to a SF object 
               print("Transform df_pixel_snowmelt to a feature collection with random geometry:")
               df_sf_tmp <- st_as_sf(x = df_tmp,                         
                                     coords = c("lon", "lat"),
                                     crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
              
              #Change sf object to an earth engine feature collection by uploading it to the asset folder
               FC_tmp <- sf_as_ee(
                 x = df_sf_tmp,
                 assetId = paste0(path_asset, "/", data_ID, "_df_pixel_snowmelt_", i),
                 overwrite = FALSE,
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
            FC_pixels_snowmelt$first()$getInfo()
            FC_pixels_snowmelt$size()$getInfo()
            
            #Inspect assets folder:
             ee_manage_quota()
             ee_manage_assetlist(path_asset)
            
        #(16): Add geometry (latitude and longitude) of each pixel_ID to FC_pixels_snowmelt
             
           #The feature collection that we want to construct should contain 'pixel_ID' and 'doy_snowmelt' as properties and should contain
           #the original geometry corresponding to each pixel_ID. So far, FC_pixels_snowmelt contains a separate feature for each pixel, where
           #each feature contains doy_snowmelt and pixel_ID as a property. However, the geometry (lan/lon) of each feature is randomly chosen. 
           #We need to make sure that the actual geometry matching each pixel_ID is added instead of this random geometry.
            
           #To obtain the corresponding geometry (lat/lon) of every pixel_ID within aoi, we sample from a single image on a 'resolution' 
           #resolution using img$sample(). This gives a single distinct feature for each pixel, including their geometry. We can then re-construct
           #the property pixel_ID for every feature in this feature collection. The resulting feature collection is called FC_pixels_distinct. 
           #The final step is then to join FC_pixels_snowmelt to FC_pixels_distinct based on an inner join with pixel_ID.
           
           #This step results in memory errors (on the server side) when the image contains more than 2.5 million pixels. In the latter case
           #we need to split up our study area in several sub-areas (i.e. create separate adjacent shapefiles using QGIS). If such errors are
           #encountered here then script "6-RGEE_TomVersluijs_PixelDateOfSnowmelt_TAY" should be run instead. 
        
            #(A): Select a single image from the image collection
             img <- MODIS_clouds_filtered$first()$select("NDSI")
            
            #(B): Create a feature collection of points within aoi (where each point corresponds to the center of a pixel in the image)
             FC_pixels_distinct <- img$sample(
              region=aoi, #Sample all pixels within aoi
              geometries=TRUE,  #Store lat/lon in geometry property
              projection=modisProjection, #set to native projection of MODIS image
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
            
             FC_pixels_distinct$first()$getInfo()
             FC_pixels_distinct$size()$getInfo()
             
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
               tryCatch({ee_manage_delete(paste0(path_asset, "/", data_ID, "_FC_pixels_snowmelt_optimized"))}, 
                        error = function(cond){return("Path did not yet exist - no folder deleted")})
               
              #Upload to asset folder:
               assetid2 <- paste0(path_asset, "/", data_ID, "_FC_pixels_snowmelt_optimized")
               task_vector2 <- ee_table_to_asset(
                   collection = FC_Combined,
                   overwrite = FALSE,
                   assetId = assetid2
                   )
               task_vector2$start()
               print("Optimize further calculations with FC_Combined by uploading it to the asset folder:")
               ee_monitoring(task_vector2, max_attempts = 1000000)
             
              #Check assets folder:
               #ee_manage_quota()
               #ee_manage_assetlist(path_asset)
              
              #Get feature collection from asset folder and create FC_pixels_snowmelt_optimized
               #assetid2=paste0("users/escape/", data_ID, "/", data_ID, "_FC_pixels_snowmelt_optimized")
               FC_pixels_snowmelt_optimized <- ee$FeatureCollection(assetid2) 
               #FC_pixels_snowmelt_optimized$first()$getInfo()
               #FC_pixels_snowmelt_optimized$size()$getInfo()
           
        #(17): Transform the Feature collection FC_pixels_snowmelt_optimized to an image (with doy_snowmelt as an image band)
               
              #(A): Reduce feature collection FC_pixels_snowmelt_optimized to an Image with a 500m resolution:    
               image_snowmelt <- FC_pixels_snowmelt_optimized$
                 filter(ee$Filter$notNull(list('doy_snowmelt')))$ #pre-filter data for nulls that cannot be turned into an image
                 filter(ee$Filter$neq(name='doy_snowmelt', value=-9999))$ #pre-filter data for -9999 values
                 reduceToImage(properties=list('doy_snowmelt'), reducer=ee$Reducer$first()$setOutputs(list('doy_snowmelt')))$
                 #reproject(crs=crs, crsTransform=NULL, scale=resolution) #ensures sure that reduceToImage above is done in the crs projection
                 reproject(crs=modisProjection, crsTransform=NULL, scale=resolution) #ensures sure that reduceToImage above is done in the MODIS projection
               
               #image_snowmelt$projection()$getInfo()
               
              #(B): Clip image along The specified shapefile.
               image_snowmelt <- image_snowmelt$clipToCollection(aoi_Shapefile)
               
              #(C): Extract day of snowmelt in year of interest for a single point
               ee_extract(x=image_snowmelt, y=coordinates_point, fun=ee$Reducer$first(), scale=resolution, sf=TRUE)
             
              #(D): Plot snowmelt day of year as a coloured image
               MODIS_image <- MODIS_clouds_filtered$filterDate(paste0(year_ID, "-07-20"), end_date)$first()#$reproject(crs=crs, crsTransform=NULL, scale=resolution)
               Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
               Map$addLayer(MODIS_image,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
               Map$addLayer(image_snowmelt,list(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')
          
              #(E): Export original image to Google Drive (takes c.a. 2 minutes):
             
                #Create task to export the original doy_snowmelt image to Google Drive  
                 task_vector3 <- ee_image_to_drive(
                  fileFormat='GeoTIFF',
                  image=image_snowmelt,
                  description=paste0(data_ID, '_PixelSnowmeltDoy_Image_DoySnowmelt'),
                  region=aoi,
                  #scale=ee$Number(resolution), #defaults to native resolution of image asset.
                  crs="EPSG:3857", #Coordinate reference system of projection of exported image
                  maxPixels=1e9, #maximum allowed number of pixels in exported image
                  dimensions=ee$Number(1024) #maximum dimension
                  )
               
                #Start and monitor export task:
                 task_vector3$start()
                 print("Export original image to Google Drive:")
                 ee_monitoring(task_vector3, max_attempts = 1000000)
                 ee_drive_to_local(task = task_vector3, dsn=paste0("Output/", data_ID, "_PixelSnowmeltDoy_Image_DoySnowmelt"))
                 
              #(F): Export RGB image to Google Drive (takes c.a. 2 minutes):
             
                #Convert original image to an RGB image:
                 image_snowmelt_RGB <- image_snowmelt$visualize(bands=c('doy_snowmelt'), min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red'))
                 #ee_print(image_snowmelt_RGB)
                 #image_snowmelt_RGB$projection()$getInfo()

                #Create task to export RGB image to Google Drive:
                 task_vector4 <- ee_image_to_drive(
                   fileFormat='GeoTIFF',
                   image=image_snowmelt_RGB,
                   description=paste0(data_ID, '_PixelSnowmeltDoy_Image_RGB'),
                   region=aoi,
                   #scale=ee$Number(resolution), #defaults to native resolution of image asset.
                   crs="EPSG:3857", #Coordinate reference system of projection of exported image
                   maxPixels=1e9, #maximum allowed number of pixels in exported image
                   dimensions=ee$Number(1024) #maximum dimension
                   )

                #Start and monitor export task:
                 print("Export RGB image to Google Drive:")
                 task_vector4$start()
                 ee_monitoring(task_vector4, max_attempts = 1000000)
                 ee_drive_to_local(task = task_vector4, dsn=paste0("Output/", data_ID, "_PixelSnowmeltDoy_Image_RGB"))
                 
        #18: Extract the date of snowmelt for all pixels within image_snowmelt (i.e. clipped by aoi_Shapefile)         
                 
              #(A): Extract the date of snowmelt at each pixel within image_snowmelt and store each pixel value as a separate feature.
              #The resulting output is a feature collection of all features (pixels) for the current image
               FC_snowmelt <- image_snowmelt$sample( #sampling is done automatically for all Bands of the input image
                 region=aoi_Shapefile, #All pixels within aoi_Shapefile will be stored as a separate feature
                 geometries=TRUE,  #if TRUE, add center of sampled pixel as the geometry property of the output feature
                 projection=modisProjection, #Set to native projection of MODIS image
                 seed=23, #Create reproducible results using the same random seed
                 dropNulls=FALSE) #If TRUE, the result is post-filtered to drop features that have a NULL value for NDSI
                 
              #(B): Make sure there is a doy_snowmelt value at each feature within the feature collection:
              #Set the band value to a no data value of -9999 for all features where the band value is NULL.
               FC_snowmelt <- FC_snowmelt$map(function(feature){
                  doy_snowmelt <- ee$List(list(feature$get('doy_snowmelt'), -9999))$reduce(ee$Reducer$firstNonNull())
                  return(feature$set("doy_snowmelt", doy_snowmelt))})
                 
              #(C): Add latitude and longitude of each pixel as a property to each feature
               FC_snowmelt <- FC_snowmelt$map(function(feature){
                 coordinates <- feature$geometry()$coordinates()
                 lon <- coordinates$get(0)
                 lat <- coordinates$get(1)
                 return(feature$set('lon', lon)$set('lat', lat))
                  })
                 
              #(D): Transform feature collection FC_snowmelt to a dataframe:
                 
                #We use ee_table_to_drive() to prevent memory limits
                 a=Sys.time()
                 task_vector5 <- ee_table_to_drive(
                   collection = FC_snowmelt,
                   description = paste0(data_ID, "_Pixel_Snowmelt_shapefile"),
                   fileFormat = "CSV",
                   selectors = c('doy_snowmelt', 'lat', 'lon')
                   )
               
                #Execute task  
                 task_vector5$start()
                 print("Transform image_snowmelt to a feature Collection of doy_snowmelt values for all pixels:")
                 ee_monitoring(task_vector5, max_attempts = 1000000)
                 exported_stats <- ee_drive_to_local(task = task_vector5, dsn=paste0("Output/", data_ID, "_Pixel_Snowmelt_shapefile"))
                 df_pixel_snowmelt_shapefile <- read.csv(exported_stats)
                 b=Sys.time()
                 print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))
      
                #Make sure each latitude/longitude combination gets its own pixel_ID (takes ca 1 minute):
                 df_pixel_snowmelt_shapefile$pixel_ID <- paste0(format(df_pixel_snowmelt_shapefile$lat, nsmall = 5), "_", format(df_pixel_snowmelt_shapefile$lon, nsmall = 5))
                 
                #Only select columns "pixel_ID" and "doy_snowmelt"
                 df_pixel_snowmelt_shapefile <- df_pixel_snowmelt_shapefile[,c("pixel_ID", "doy_snowmelt")]
                 
              #(E): Save dataframe
                write.csv(df_pixel_snowmelt_shapefile, file=paste0(here(), "/Output/", data_ID, "_Pixel_Snowmelt_shapefile.csv"), quote = FALSE, row.names=FALSE)
                 
        #(19): Save workspace
         #save.image(paste0(here(), "/Output/", data_ID, "_Backup_Workspace_PixelDateOfSnowmelt.RData"))
         
    #The snowmelt image is now completed and can be downloaded as .tif file from the RGEE_backup folder on the specified Google Drive.
    #The optional code below uses the generated snowmelt image to extract the snowmelt day of year at specific points of interest.

         
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################

  # #OPTIONAL 1: Extract snowmelt data from image_snowmelt at points of interest
  # 
  #   #(0): Clear workspace
  #    rm(list=ls())
  # 
  #   #(1): Load packages
  #    library(pacman)
  #    library(rgee)
  #    p_load(sf,
  #           ggplot2,
  #           mgcv,
  #           dplyr)
  # 
  #    #(2): Define ggplot2 plotting theme
  #     theme_tom <- function(){
  #       theme_classic() %+replace%
  #         theme(axis.title = element_text(size=18),
  #               axis.text = element_text(size=16),
  #               legend.position = "none",
  #               strip.background = element_rect(fill = "white", colour = "black", size = rel(2)),
  #               complete = TRUE)}
  # 
  #   #(3): Initialize earth engine and google drive
  #    rgee::ee_Initialize(user="tom.versluijs@gmail.com", drive = T)
  #    #ee_check()
  # 
  #   #(4): Specify parameters used in the analysis
  # 
  #     #(a): Unique data identifiers:
  #       area_name="ZAC"
  #       year_ID="2019"
  #       data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
  # 
  #     #(b): Approximate central point of study area
  #       coordinates_point <- c(-20.49, 74.49)
  # 
  #     #(c): Coordinate reference system used for calculations
  #       #EPSG:4326 is recommended for areas spanning multiple UTM zones, but increased computation time (i.e. spherical coordinate system).
  #       #EPSG:326XX is results in reduced computation time for areas located within a single UTM zone (i.e. planar coordinate system).
  #       crs <- "EPSG:4326"
  # 
  #     #(d): Date range for all images considered in analysis
  #       start_date <- paste0(year_ID, "-03-15") #choose date (well) before the first snowmelt occurs within the study site
  #       end_date <- paste0(year_ID, "-08-15") #choose date (well) after last date of tracking
  # 
  #     #(e): Buffer radius around each point of interest
  #       Buffer_radius_m=0
  # 
  #     #(f): Resolution of sampling in meters
  #       resolution=500 # native resolution of satellite image (sentinel-2=10m, MODIS=500m)
  # 
  #     #(g): Define which MODIS cloud masking algorithm has been used ("PGE11", "MOD35", or "Combined)
  #       MODIS_cloud_masking_algorithm = 'PGE11'
  #       data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)
  # 
  #     #(h): Name of shapefile
  #       shapefile <- "ZAC_Outline_EPSG4326.shp"
  # 
  #   #(5): Create a unique Asset folder
  #      path_asset <- paste0(ee_get_assethome(), "/", data_ID)
  # 
  #   #(6): Calculate date range in day of year format
  #      start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
  #      end_date_doy <- as.numeric(strftime(end_date, format = "%j"))
  # 
  #   #(7): Read aoi_Shapefile shapefile in root folder
  #      root_fldr <- here()
  #      aoi_Shapefile <- st_read(paste0(root_fldr, "/Input/Shapefiles/", shapefile), quiet=T)
  #      aoi_Shapefile <- st_transform(aoi_Shapefile, crs="EPSG:4326")
  #      aoi_Shapefile <- sf_as_ee(aoi_Shapefile)
  #      aoi_Shapefile <- ee$FeatureCollection(aoi_Shapefile)
  # 
  #   #(8): Reload MODIS sinusoidal projection
  #      MODIS_col<-ee$ImageCollection('MODIS/006/MOD09GA')
  #      point <- ee$Geometry$Point(coordinates_point[1], coordinates_point[2])
  #      MODIS_col <- MODIS_col$filterBounds(point)$filterDate(start_date, end_date)
  #      modisProjection <-MODIS_col$first()$select("sur_refl_b01")$projection()
  #      #modisProjection$getInfo()
  # 
  #   #(9): Re-construct generated snowmelt image ('image_snowmelt')
  # 
  #      #Load feature collection FC_pixels_snowmelt_optimized
  #      assetid2=paste0(path_asset, "/", data_ID, "_FC_pixels_snowmelt_optimized")
  #      FC_pixels_snowmelt_optimized <- ee$FeatureCollection(assetid2)
  # 
  #      #Reduce feature collection FC_pixels_snowmelt_optimized to an Image with MODIS sinusoidal projection
  #       image_snowmelt <- FC_pixels_snowmelt_optimized$
  #         filter(ee$Filter$notNull(list('doy_snowmelt')))$ #pre-filter data for nulls that cannot be turned into an image
  #         filter(ee$Filter$neq(name='doy_snowmelt', value=-9999))$ #pre-filter data for -9999 values
  #         reduceToImage(properties=list('doy_snowmelt'), reducer=ee$Reducer$first()$setOutputs(list('doy_snowmelt')))$
  #         reproject(crs=modisProjection, crsTransform=NULL) #ensures sure that reduceToImage above is done in the MODIS projection
  #         #reproject(crs=crs, crsTransform=NULL, scale=resolution) #ensures sure that reduceToImage above is done in the crs projection
  # 
  #       #Clip image along The loaded shapefile.
  #        image_snowmelt <- image_snowmelt$clipToCollection(aoi_Shapefile)
  # 
  #      #Plot snowmelt image on a map. Note that 'Map' automatically converts to mercator projection!
  #       Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
  #       Map$addLayer(image_snowmelt,list(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')
  # 
  #   #(10): Define a featurecollection of points corresponding to movement trajectories of chicks in Zackenberg in 2019
  #     ZAC_Paths <- read.csv(paste0("Input/ZAC19_ChickMovement.csv"), header=T)[,c("ChickID", "LON_x", "LAT_y", "DateTime", "Hatch_Date")]
  #     colnames(ZAC_Paths) <- c("ChickID", "LON_x", "LAT_y", "DateTime", "Hatch_Date")
  # 
  #     #Add a unique LocationID to every unique lat/lon combination in ZAC_Paths
  #      Locations <- unique(ZAC_Paths[,c("LON_x", "LAT_y")])
  #      Locations$LocationID <- paste0("Location_", 1:nrow(Locations))
  #      ZAC_Paths <- left_join(ZAC_Paths, Locations, by=c("LON_x", "LAT_y"))
  # 
  #     #Transform Locations to an ee-object with MODIS sinusoidal projection
  #      Locations <- st_as_sf(Locations, coords = c("LON_x", "LAT_y"), crs=crs)
  #      Locations <- st_transform(Locations, crs=crs) #as an intermediate step
  #      Locations <- sf_as_ee(Locations[,c("LocationID", "geometry")], proj='epsg:4326') #as an intermediate step
  # 
  #     #Apply transformation of all locations to MODIS sinusoidal projection
  #      transform_to_modis <- function(feature){
  #        transformed_feature <- feature$transform(proj=modisProjection)
  #        return(transformed_feature)
  #        }
  #      Locations <- Locations$map(transform_to_modis)
  # 
  #     #Add buffer around all point locations
  #      Buffer_radius_m=Buffer_radius_m
  #      if(Buffer_radius_m>0){
  #          bufferBy <- function(Buffer_radius_m) {
  #            return(function(feature) {
  #              return(feature$buffer(Buffer_radius_m))
  #              })
  #            }
  #          Locations <- Locations$map(bufferBy(Buffer_radius_m))
  #          }
  # 
  #      #Locations$first()$getInfo()
  # 
  #      #Locations is now a feature collection. A feature in our analysis comprises the feature type (point in case no buffer
  #      #was added, or a circular polygon with radius Buffer_radius_m with a center the original point location), Location_ID
  #      #and the feature's geometry as lat and lon coordinates. A feature collection corresponds to a collection of such features
  #      #(i.e. a collection of different Locations, with their respective feature types, LocationIDs, DateTimes and coordinates).
  # 
  #     #Plot feature collection on a map (i.e. locations as points). Note that 'this'Map' automatically converts to Mercator projection!
  #      Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
  #      Map$addLayer(image_snowmelt,list(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')+
  #      Map$addLayer(Locations, list(color="red"), paste0("MovementPositions_", year_ID))
  # 
  #   #(11): Extract Sentinel-2 snowmelt_doy values for points of interest from image_snowmelt (output a feature collection)
  # 
  #     #Extract snowmelt_doy values within the current image at each point in the Locations feature collection. We add this doy
  #     #as a property to each feature (i.e. movement position). The resulting output is a feature collection for the current img
  #     #where the snowmelt_doy values are added as a property for each movement position.
  #      FC_image <- image_snowmelt$reduceRegions(collection = Locations,
  #                                               reducer = ee$Reducer$mean()$setOutputs(list('Snowmelt_doy')),
  #                                               #scale=resolution,
  #                                               crs=modisProjection,
  #                                               crsTransform=NULL)
  # 
  #     #FC_image$getInfo()
  # 
  #     #reduceRegion does not return any values when the point is masked by clouds or does not fall within
  #     #the boundaries of the image. Therefore we manually have to add an empty snowmelt_doy property to those features.
  #     #We therefore set the snowmelt_doy value to a no data value of -9999 for all features where the snowmelt_doy value is NULL.
  #      FC_image <- FC_image$map(function(feature){
  #          snow <- ee$List(list(feature$get('Snowmelt_doy'), -9999))$reduce(ee$Reducer$firstNonNull())
  #          return(feature$set("Snowmelt_doy", snow))
  #          })
  # 
  #    #(12): Transform feature collection with Band values to a dataframe
  # 
  #      #Extract the data for all features (locations) within all images of the image collection (summer year_ID) for the area of interest (aoi)
  #       snow_aoi <- unlist(FC_image$aggregate_array('Snowmelt_doy')$getInfo())
  #       location_aoi <- unlist(FC_image$aggregate_array('LocationID')$getInfo())
  # 
  #      #Store all data in a dataframe
  #       Locations_BandValues <- data.frame(Snowmelt_doy=round(as.numeric(snow_aoi), 5),
  #                                           LocationID=location_aoi)
  # 
  #      #Change -9999 to NA
  #       Locations_BandValues$Snowmelt_doy[Locations_BandValues$Snowmelt_doy < -9000] <- NA
  # 
  #      #Sort dataframe by LocationID
  #       index <- with(Locations_BandValues, order(LocationID))
  #       Locations_BandValues <- Locations_BandValues[index,]
  # 
  #    #(13) Add Snowmelt_doy values per location to the corresponding date-location combination in ZAC_Paths
  #        ZAC_Paths$Date_doy <- as.numeric(strftime(ZAC_Paths$DateTime , format = "%j"))
  #        ZAC_Paths <- left_join(ZAC_Paths, Locations_BandValues, by=c("LocationID"))
  #        write.csv(ZAC_Paths, paste0(here(), "/Output/", data_ID, "_ChickMovement_Snowmelt_doy_Buffer", Buffer_radius_m, ".csv"), row.names = FALSE)
  # 
  #    #(14) Plot each individuals path through a landscape of Snowmelt doy values
  #        p1 <- ggplot()+
  #          geom_line(data=ZAC_Paths, aes(x=Date_doy, y=Snowmelt_doy ), col="#1620de", lwd=1.25)+
  #          geom_abline()+
  #          facet_wrap(~ChickID, ncol=ceiling(length(unique(ZAC_Paths$ChickID))^0.5))+
  #          geom_smooth(data=ZAC_Paths, aes(x=Date_doy, y=Snowmelt_doy), method="lm", formula=y~x, col="red", size=1.25, se=T)+
  #          geom_point(data=ZAC_Paths, aes(x=Date_doy, y=Snowmelt_doy ), col="black")+
  #          xlab("Day of year") +
  #          ylab("Day of snowmelt") +
  #          theme_tom()
  #
  #          ggsave(plot=p1, paste0(here(), "/Output/", data_ID, "_MovementPaths_Snowmelt_doy_Grid_Buffer", Buffer_radius_m, ".pdf"), width=10, height = 8)
  # 
  #     #(15): There do not seem to be any clear trajectories in individuals moving up to areas with a later snowmelt (except maybe 8218556)
  #     #We therefore calculate the average encountered Snowmelt_doy values along each individual's foraging trajectory
  #        ZAC_Chicks_MeanBandvalue <- ZAC_Paths[!is.na(ZAC_Paths$Snowmelt_doy),] %>% group_by(ChickID) %>% summarise(mean_Snowmelt_doy=mean(Snowmelt_doy),
  #                                                                                                                   mean_doy=mean(Date_doy))
  # 
  #        #Load estimated growth parameters per chick
  #        ZAC_ChickGrowth <- read.csv(paste0("Input/ZAC19_ChickGrowth.csv"), header=T)
  #        ZAC_ChickGrowth <- ZAC_ChickGrowth[ZAC_ChickGrowth$max_age>=5,]
  # 
  #        #Combine growth parameters with mean Band value during foraging
  #        ZAC_Chicks_MeanBandvalue <- left_join(ZAC_Chicks_MeanBandvalue, ZAC_ChickGrowth, by="ChickID")
  # 
  #        #Combine growth parameters with mean Band value during foraging
  #        p2 <- ggplot()+
  #          geom_smooth(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_Snowmelt_doy, y=d, size=max_age), method="lm", formula=y~x, col="red", size=1.25)+
  #          geom_point(data=ZAC_Chicks_MeanBandvalue, aes(x=mean_Snowmelt_doy, y=d, size=max_age), col="black")+
  #          xlab("Mean day of snowmelt across foraging locations") +
  #          ylab("Upper asymptote bodymass (gram)") +
  #          theme_tom()
  #
  #          ggsave(plot=p2, paste0(here(), "/Output/", data_ID, "_GrowthParameter_Mean_Snowmelt_doy_Buffer", Buffer_radius_m, ".pdf"), width=10, height = 8)
  # 
  #        #Construct linear mixed models:
  #        lmm0 <- lme(d ~ 1, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #        lmm1 <- lme(d ~ mean_Snowmelt_doy, random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #        lmm2 <- lme(d ~ mean_Snowmelt_doy + I(mean_Snowmelt_doy^2), random=~1|Brood, data=ZAC_Chicks_MeanBandvalue, method="ML")
  #        MuMIn::AICc(lmm0, lmm1, lmm2)
  #        rm(lmm0, lmm1, lmm2)
  # 
  #        #Calculate mean visiting date vs mean date of snowmelt for all chicks (55 days!)
  #        mean(ZAC_Chicks_MeanBandvalue$mean_doy - ZAC_Chicks_MeanBandvalue$mean_Snowmelt_doy)
    
#####################################################################################################################################          
       
# #OPTIONAL 2: Add movement trajectories to Snowmelt map:
# 
#   #(0): Clear workspace
#    rm(list=ls())
# 
#   #(1): Load packages
#    library(pacman)
#    library(rgee)
#    p_load(stars,
#           ggplot2,
#           colorspace,
#           ggnewscale,
#           future)
# 
#    #Define ggplot2 plotting theme
#     theme_tom <- function(){
#        theme_classic() %+replace%
#              theme(axis.title = element_text(size=18),
#                    axis.text = element_text(size=16),
#                    legend.position = "none",
#                    strip.background = element_rect(fill = "white", colour = "black", size = rel(2)),
#                    complete = TRUE)}
# 
#   #(2): Initialize earth engine and google drive
#     rgee::ee_Initialize(user="tom.versluijs@gmail.com", drive = T)
#     #ee_check()
# 
#   #(3) Specify a unique data identifier
#     area_name="ZAC"
#     year_ID="2019"
#     data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
# 
#     #Define which MODIS cloud masking algorithm has been used ("PGE11", "MOD35", or "Combined)
#     MODIS_cloud_masking_algorithm = 'PGE11'
#     data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)
# 
#   #(4): Load the .tif file for the snowmelt image. This .tif file has to be downloaded manually from my google drive
#   #and placed in the local working directory (/Input) of this script.
#    tif=paste0("Input/", data_ID, "_PixelSnowmeltDoy_Image_RGB.tif")
#    x=read_stars(tif)
#    x_RGB <- st_rgb(x,3)
# 
#   #(5): Plot snowmelt map
#    p3 <- ggplot() +
#            geom_stars(data=x_RGB)+
#            scale_fill_identity()+
#            labs(x = 'Longitude', y = 'Latitude') +
#            theme_tom() +
#            geom_rect(aes(xmin=508920, xmax=519800, ymin=8261350,ymax=8272600), alpha=0.25, fill="black")
# 
#   #(6): Plot trajectories of chicks with Age>=5 (same as used in chick growth analysis):
#     ChickID=c("8218553", "8218556", "8218557","8218558","8218560","8218571","8218572","8218580","8218581", "8218582")
#     df <- read.csv(paste0("Input/ZAC19_ChickMovement.csv"), header=TRUE)
#     df <- df[order(df$DateTime),]
#     ChickID_all <- data.frame(ChickID=unique(df$ChickID))
#     ChickID_all$index <- 1:nrow(ChickID_all)
#     ChickID_path <- data.frame(ChickID=ChickID, col_start=NA, col_end=NA)
#     n <- nrow(ChickID_path)
#     ChickID_path$col_end <- c("red","green","yellow","blue","orchid","orange", "turquoise1", "bisque", "seagreen3","lightskyblue")[1:length(ChickID)]
#     ChickID_path$col_start <- lighten(ChickID_path$col_end, amount=0.4)
#     ChickID_all <- merge(ChickID_all, ChickID_path, by="ChickID", all.x=T)
#     ChickID_all$col_start[is.na(ChickID_all$col_start)] <- "#808080"
#     ChickID_all$col_end[is.na(ChickID_all$col_end)] <- "#808080"
#     ChickID_all <- ChickID_all[order(ChickID_all$index),c(1,3,4)]
# 
#     #Plot movement data on top of snowmelt map
#      for(i in 1:length(ChickID)){
#        p3 <- p3+
#          geom_point(data=df[df$ChickID==ChickID[i],], aes(x=LON_x, y=LAT_y, col=Age), size = 2)+
#          geom_path(data=df[df$ChickID==ChickID[i],], aes(x=LON_x, y=LAT_y, color=Age), size=1.25)+
#          scale_color_gradient(low=ChickID_all$col_start[ChickID_all$ChickID==ChickID[i]],
#                               high=ChickID_all$col_end[ChickID_all$ChickID==ChickID[i]],
#                               breaks = seq(0, max(na.omit(df$Age)), by = 1))+
#          new_scale_color()
#          }
#      ggsave(plot=p3, paste0(here(), "/Output/", data_ID, "_Chick_Movement_Snowmelt.pdf"), width=13 , height=12 , units="in")
         
         
##############################################################################################################################################################         
             
          
