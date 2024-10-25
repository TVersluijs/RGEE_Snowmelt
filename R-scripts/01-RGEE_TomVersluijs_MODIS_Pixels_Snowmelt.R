##################################################################################################################################

#Extract MODIS satellite data and calculate the date of snow melt for every 500mx500m pixel in an area of interest (shapefile). 
#The user can specify whether clouds and permanent waterbodies need to be masked. Snow melt is calculated per pixel by fitting 
#a GAM through the average NDSI data and by extracting the moment this GAM crosses a user specified NDSI threshold. After 
#creating the snow melt map, script "10-RGEE_TomVersluijs_S2_ExtractSnowFraction.R" can be used to extract time series of the 
#fraction of snow cover for points/polygon(s) of interest from this map.

#Copyright Tom Versluijs 2024-10-25. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R". 
#Note that a GoogleDrive is required. Important: make sure to run this script from within the "RGEE_Snowmelt.Rproj" project file.

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
       p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr, tidyr, foreach, parallel, doSNOW, gridExtra)       

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

     #MODIS dataset
     MODIS_dataset <- "MODIS/061/MOD09GA"
       
     #Specify resolution of images in meters
     resolution=500 #default maximum resolution for MODIS = 500m

   #(b) Area of interest

     #Specify name of study area
     area_name <- "ZAC" #max length three characters

     #Name of Shapefile (located in '/Input/Shapefiles' folder)
     #Follow the guide "Manual_CreateShapefilePolygons.docx" when creating this shapefile.
     shapefile <- "ZAC_Outline_EPSG4326.shp"

     #Coordinate reference system used for calculations
     #EPSG:4326 is recommended
     #EPSG:326XX might result in reduced computation time for areas located within a single UTM zone (i.e. planar coordinate system).
     crs <- "EPSG:4326"

   #(c) Dates

     #Specify the year of interest:
     year_ID <- "2023"
     
     #Date range for all images considered in analysis
     start_date <- paste0(year_ID, "-03-15") #choose date (well) before the first snow melt occurs within the study area
     end_date <- paste0(year_ID, "-09-15") #choose date (well) after last snow melt occurs within the study area

   #(d) Snow detection
     
     #NDSI threshold above which a pixel is perceived as snow (only a single threshold allowed)
     NDSI_threshold=0.4

   #(e): Cloud masking

     #Should clouds be masked from the analysis (default=TRUE).
     mask_clouds=TRUE
     
     #Define which MODIS cloud masking algorithm has been used ("PGE11", "MOD35", or "Combined)
     MODIS_cloud_masking_algorithm = 'PGE11' #default is PGE11

     #Maximum fraction of cloud cover allowed in each image 
     max_cloud_fraction=1.0 #1.0 equals 100% cloud cover
     
     #Specify the spatial resolution of the cloud masking algorithm
     resolution_qaband=1000 #default maximum resolution for MODIS = 1000m
     
   #(f): Water masking

     #Should permanent waterbodies be masked from the analysis using a 250m resolution Terra Land Water Mask dataset (default=TRUE).
     mask_water=TRUE
     
   #(g): GAM sequential outlier filtering
       
     #parameters for auxiliary function to filter outliers from GAM using a two-step approach
     outlier_removal=TRUE #should sequential outlier removal be employed when fitting GAMs to the data
     outlier_thresh_1=0.4 #first threshold in outlier removal (relative to y-range of data, default=0.4)
     outlier_thresh_2=0.2 #second threshold in outlier removal (relative to y-range of data, default=0.2)
     
     #Specify the degree of GAM smoothing by setting the 'k' parameter
     gam_k_outlier=10 #Number of knots when filtering outliers (default=10)
     gam_k=25 #Number of knots when making model predictions (default=25). 
     #Larger values result in a more precise GAM-fit, but at a cost of computation time
       
     #Specify whether plots of the GAM fit per pixel should be generated as a pdf file
     pixel_gam_plots=FALSE #significantly increases computation time
     
   #(h): Counts of unmasked pixels
     
     #Should counts of the number of unmasked pixels per day of year within the shapefile area be conducted (increases computation time)
     pixel_counts=TRUE
     
     
##################################################################################################################################
         
#III: Define some additional parameters (automated)
         
##################################################################################################################################
         
      #(5) Create a unique dataID and asset folder for storing the generated datafiles
       
        #Create a unique data_ID
         if(nchar(area_name)>3){area_name <- substr(area_name, start = 1, stop = 3)}
         data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
         data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)
          
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
      
        #Create output folder
         if(dir.exists(paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt"), recursive = TRUE)}
         
        #Save all parameters and their values in the environment to a text file 
         file_conn <- file(paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Parameters.txt"), "w")
         for (obj in setdiff(ls(), lsf.str())) {cat(paste(obj, "=", get(obj)), file = file_conn) ; cat("\n", file = file_conn)}
         close(file_conn)
         
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
           scale= resolution,
           maxPixels= 1E13)
         paste0('Size of study area calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')

      #(8) Extract MODIS satellite Surface Reflectance images for specified daterange and general area of interest
       
       #Note that specifying a region of interest using a polygon does not function properly for the MODIS data. 
       #Instead, providing an initial point and then afterwards clipping the image to a desired area of interest works better.
       
        #Specify starting and ending date  
         start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
         end_date_doy <- as.numeric(strftime(end_date, format = "%j"))
       
        #Extract MODIS satellite data
         MODIS_col <- ee$ImageCollection(MODIS_dataset)
         MODIS_col <- MODIS_col$
           filterBounds(aoi)$
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
        #  tryCatch({browseURL(MODIS_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$getVideoThumbURL(videoArgs)) }, error = function(cond){return("Too many pixels. Reduce dimensions.")})
        
        # #Create a timelapse video of NDSI band  
        #  palette=c('black', '0dffff', '0524ff', 'ffffff')
        #  visFun_NDSI <-  function(img) {
        #    return(img$visualize(bands='NDSI', min=-1, max=1.5, palette=palette)$copyProperties(img, img$propertyNames()))
        #    }
        #  MODIS_snow_RGB <- MODIS_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$map(visFun_NDSI)
        #  videoArgs <- list(dimensions=380, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
        #  tryCatch({browseURL(MODIS_snow_RGB$getVideoThumbURL(videoArgs)) }, error = function(cond){return("Too many pixels. Reduce dimensions.")})
      
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
        #  tryCatch({browseURL(MODIS_clouds_filtered$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$getVideoThumbURL(videoArgs)) }, error = function(cond){return("Too many pixels. Reduce dimensions.")})
     
        # #Create a timelapse video of the cloud filtered/masked NDSI band (for debugging) 
        #  palette=c('black', '0dffff', '0524ff', 'ffffff')
        #  visFun_NDSI <-  function(img) {
        #    return(img$visualize(bands='NDSI', min=-1, max=1.5, palette=palette)$
        #             copyProperties(img, img$propertyNames()))}
        #  MODIS_snow_masked_RGB <- MODIS_clouds_filtered$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$map(visFun_NDSI)
        #  videoArgs <- list(dimensions=380, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
        #  tryCatch({browseURL(MODIS_snow_masked_RGB$getVideoThumbURL(videoArgs)) }, error = function(cond){return("Too many pixels. Reduce dimensions.")})

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
      
      #(B): Extract permanent waterbodies from a 250m resolution Terra Land Water Mask dataset
    
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
        # tryCatch({browseURL(MODIS_clouds_filtered$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
      
      }
   if(mask_water==FALSE){
      
      #(A): print message   
      print("Water masking = FALSE")
      
    }
     
   #Note that MODIS_clouds_filtered is not yet clipped by aoi_Shapefile (only by aoi)! 
         
##################################################################################################################################

#VII: Count the total number of unmasked pixels per doy within aoi_Shapefile

##################################################################################################################################            
   
   #(14): Count the total number of unmasked pixels and the total number of pixels per doy within aoi_Shapefile
   if(pixel_counts==TRUE){
   
     #(A): Add pixel counts within the area of interest to each separate image by mapping the pixel count functions over the image collection
      MODIS_clouds_filtered <- MODIS_clouds_filtered$
        #Count number of unmasked pixels within aoi_Shapefile
        map(AddPixelCount)$
        #Add NULL to those images in which pixel count could not be calculated
        map(AddNULLPixelCount)
  
     #(B): Extract pixel_counts of all images in image collection for aoi_Shapefile
    
       #create a current timestamp to prevent identical names on Google Drive
        current_timestamp0 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
       
       #We use ee_table_to_drive() to prevent memory limits
        task_vector0 <- ee_table_to_drive(
          collection = MODIS_clouds_filtered,
          description = paste0(current_timestamp0, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char,  "_Data_Pixel_Counts_polygon"),
          fileFormat = "CSV",
          selectors = c('doy', 'unmasked', 'total')
          )
       
       #Monitor the task
        task_vector0$start()
        print("Count the number of unmasked pixels within the shapefile per doy:")
        ee_monitoring(task_vector0, quiet=T, max_attempts=1000000)
       
       #Export results to local folder
        exported_stats <- ee_drive_to_local(task = task_vector0, dsn=paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Counts_polygon"))
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
      write.csv(df_pixelcount, file=paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Counts_polygon.csv"), quote=FALSE, row.names=FALSE)
      
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
      pdf(paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Plot_Pixel_Counts_polygon.pdf"), width=12, height=8)
      print(p_pixelcounts)
      dev.off()
      
   }
     
##################################################################################################################################
        
#VIII: Calculate the date of snowmelt for every pixel within the study area by fitting a GAM through the NDSI data
        
##################################################################################################################################            
         
      #(15): Calculate the date of snowmelt (NDSI <= NDSI_threshold) for every pixel within the study area (bounding box!) 
       
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

           #create a current timestamp to prevent identical names on Google Drive
            current_timestamp1 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
            
           #We use ee_table_to_drive() to prevent memory limits
            a=Sys.time()
            task_vector1 <- ee_table_to_drive(
              collection = FC_merged,
              description = paste0(current_timestamp1, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char,  "_Pixel_NDSI_bbox"),
              fileFormat = "CSV",
              selectors = c('NDSI', 'Date', 'lat', 'lon')
              )

            task_vector1$start()
            print("Transform each image to a feature Collection of NDSI values for all pixels:")
            ee_monitoring(task_vector1, max_attempts = 1000000)

            exported_stats <- ee_drive_to_local(task = task_vector1, dsn=paste0("Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_NDSI_bbox"))
            df_pixel_ndsi <- read.csv(exported_stats)
            b=Sys.time()
            print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))

           # #Load dataframe (takes ca 2 minutes):
           #  df_pixel_ndsi <- read.csv(paste0(timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_NDSI_bbox.csv"))

           #Add day of year
            df_pixel_ndsi$doy <- as.numeric(format(as.POSIXct(df_pixel_ndsi$Date, format = "%Y-%m-%d %H:%M:%S"), "%j"))          
            
           #Make sure each latitude/longitude combination gets its own pixel_ID (takes ca 1 minute):
            df_pixel_ndsi$pixel_ID <- paste0(format(round(df_pixel_ndsi$lat, 5), nsmall = 5), "_", format(round(df_pixel_ndsi$lon, 5), nsmall = 5))

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
            f_detect_threshold_date_parallel <- f_detect_threshold_date_parallel #sourced

           #Run the 'f_detect_threshold_date_parallel' function over all data subsets, combine the results and save the resulting dataframe and plots
            print("Calculate the date of snowmelt for each pixel within aoi:")
            results <- lapply(1:length(pixelIDs_split), FUN=f_detect_threshold_date_parallel, 
                              pixelIDs_split=pixelIDs_split, df_pixel_y=df_pixel_ndsi, pixel_ID_column="pixel_ID",
                              y="NDSI", x="doy", pixel_gam_plots=pixel_gam_plots, y_threshold=NDSI_threshold)
                            
            #Clean up the cluster after finishing the parallel runs
            stopCluster(cl)
            
            #Turn parallel processing off and run sequentially again after this point
            registerDoSEQ()
            
           #Store date of snowmelt per pixel within aoi (bounding box!) as a dataframe
            df_pixel_snowmelt <- lapply(results, "[[", 1)
            df_pixel_snowmelt <- as.data.frame(do.call(rbind, do.call(c, df_pixel_snowmelt)))
            colnames(df_pixel_snowmelt)[colnames(df_pixel_snowmelt)=="x_threshold"] <- "doy_snowmelt"
            #write.csv(df_pixel_snowmelt, file=paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_Snowmelt_bbox.csv"), quote = FALSE, row.names=FALSE)
            #df_pixel_snowmelt <- read.csv(file=paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_Snowmelt_bbox.csv"), header=TRUE)  
            
           #Store GAM plots per pixel within aoi (bounding box!) as a list
            if(pixel_gam_plots==TRUE){
              plot_pixel_snowmelt <- lapply(results, "[[", 2)
              plot_pixel_snowmelt <- do.call(c, plot_pixel_snowmelt)
              }
        
           #The dataframe df_pixel_snowmelt now contains the date of snowmelt for each individual pixel_ID (514584 10mx10m pixels)
           #withing the aoi (i.e. the defined bounding box). To be able to plot these data we transform this dataframe to a feature 
           #collection and then transform this feature collection to an image.
            
       #(16): Transform df_pixel_snowmelt to a feature collection with random geometry
            
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
                 assetId = paste0(path_asset, "/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_df_pixel_snowmelt_", i),
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
            
        #(17): Add geometry (latitude and longitude) of each pixel_ID to FC_pixels_snowmelt
             
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
              scale=resolution, #sampling resolution in meters
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
             
              #create a current timestamp to prevent identical names on Google Drive
               current_timestamp2 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
               
              #Delete FC_pixels_snowmelt_optimized if it already occured in the asset folder:
               tryCatch({ee_manage_delete(paste0(path_asset, "/", current_timestamp2, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_FC_pixels_snowmelt_optimized"))}, 
                        error = function(cond){return("Path did not yet exist - no folder deleted")})
               
              #Upload to asset folder:
               assetid2 <- paste0(path_asset, "/", current_timestamp2, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_FC_pixels_snowmelt_optimized")
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
              
              #Save assetid2 for future downloading of FC_pixels_snowmelt_optimized
               saveRDS(object=assetid2, file=paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Variable_AssetID.Rds"))
              
              #Get feature collection from asset folder and create FC_pixels_snowmelt_optimized
               #assetid2=paste0(path_asset, "/", current_timestamp2, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_FC_pixels_snowmelt_optimized")
               FC_pixels_snowmelt_optimized <- ee$FeatureCollection(assetid2) 
               #FC_pixels_snowmelt_optimized$first()$getInfo()
               #FC_pixels_snowmelt_optimized$size()$getInfo()
           
        #(18): Transform the Feature collection FC_pixels_snowmelt_optimized to an image (with doy_snowmelt as an image band)
               
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
             
                #create a current timestamp to prevent identical names on Google Drive
                 current_timestamp3 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
               
                #Create task to export the original doy_snowmelt image to Google Drive  
                 task_vector3 <- ee_image_to_drive(
                  fileFormat='GeoTIFF',
                  image=image_snowmelt,
                  description=paste0(current_timestamp3, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt'),
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
                 ee_drive_to_local(task = task_vector3, dsn=paste0("Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_Image_DoySnowmelt"))
                 
              #(F): Export RGB image to Google Drive (takes c.a. 2 minutes):
             
                #Convert original image to an RGB image:
                 image_snowmelt_RGB <- image_snowmelt$visualize(bands=c('doy_snowmelt'), min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red'))
                 #ee_print(image_snowmelt_RGB)
                 #image_snowmelt_RGB$projection()$getInfo()

                #create a current timestamp to prevent identical names on Google Drive
                 current_timestamp4 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
                 
                #Create task to export RGB image to Google Drive:
                 task_vector4 <- ee_image_to_drive(
                   fileFormat='GeoTIFF',
                   image=image_snowmelt_RGB,
                   description=paste0(current_timestamp4, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt_RGB'),
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
                 ee_drive_to_local(task = task_vector4, dsn=paste0("Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_Image_DoySnowmelt_RGB"))
                 
        #(19): Extract the date of snowmelt for all pixels within image_snowmelt (i.e. clipped by aoi_Shapefile)         
                 
              #(A): Extract the date of snowmelt at each pixel within image_snowmelt and store each pixel value as a separate feature.
              #The resulting output is a feature collection of all features (pixels) for the current image
               FC_snowmelt <- image_snowmelt$sample( #sampling is done automatically for all Bands of the input image
                 region=aoi_Shapefile, #All pixels within aoi_Shapefile will be stored as a separate feature
                 geometries=TRUE,  #if TRUE, add center of sampled pixel as the geometry property of the output feature
                 projection=modisProjection, #Set to native projection of MODIS image
                 scale=resolution, #sampling resolution in meters
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
                 
                #create a current timestamp to prevent identical names on Google Drive
                 current_timestamp5 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
               
                #We use ee_table_to_drive() to prevent memory limits
                 a=Sys.time()
                 task_vector5 <- ee_table_to_drive(
                   collection = FC_snowmelt,
                   description = paste0(current_timestamp5, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Snowmelt_polygon"),
                   fileFormat = "CSV",
                   selectors = c('doy_snowmelt', 'lat', 'lon')
                   )
               
                #Execute task  
                 task_vector5$start()
                 print("Transform image_snowmelt to a feature Collection of doy_snowmelt values for all pixels:")
                 ee_monitoring(task_vector5, max_attempts = 1000000)
                 exported_stats <- ee_drive_to_local(task = task_vector5, dsn=paste0("Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Snowmelt_polygon"))
                 df_pixel_snowmelt_shapefile <- read.csv(exported_stats)
                 b=Sys.time()
                 print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))
      
                #Make sure each latitude/longitude combination gets its own pixel_ID (takes ca 1 minute):
                 df_pixel_snowmelt_shapefile$pixel_ID <- paste0(format(round(df_pixel_snowmelt_shapefile$lat, 5), nsmall = 5), "_", format(round(df_pixel_snowmelt_shapefile$lon, 5), nsmall = 5))
                 
                #Only select columns "pixel_ID" and "doy_snowmelt"
                 df_pixel_snowmelt_shapefile <- df_pixel_snowmelt_shapefile[,c("pixel_ID", "doy_snowmelt")]
                 
                #Store the pixelID of all pixels within aoi_Shapefile        
                 pixelIDs_shapefile <- unique(df_pixel_snowmelt_shapefile$pixel_ID) 
                 
              #(E): Save dataframe with the date of snowmelt for all pixels within aoi_Shapefile (i.e. buffer zone)
                write.csv(df_pixel_snowmelt_shapefile, file=paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Snowmelt_polygon.csv"), quote = FALSE, row.names=FALSE)
                 
              #(F): Store NDSI timeseries for all pixels within aoi Shapefile (i.e. buffer zone)
                df_pixel_ndsi_shapefile <- df_pixel_ndsi[df_pixel_ndsi$pixel_ID %in% pixelIDs_shapefile,]
                write.csv(df_pixel_ndsi_shapefile, file=paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_Pixel_NDSI_polygon.csv"), quote = FALSE, row.names=FALSE)
                unlink(paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Pixel_NDSI_bbox.csv"))
                
              #(G): Store plots of NDSI timeseries for all pixels within aoi_Shapefile (i.e. buffer zone)
                if(pixel_gam_plots==TRUE){
                  plot_pixel_snowmelt_shapefile <- plot_pixel_snowmelt[which(pixelIDs %in% pixelIDs_shapefile)]
                  plots_per_page = 25
                  plot_pixel_snowmelt_shapefile <- split(plot_pixel_snowmelt_shapefile, ceiling(seq_along(plot_pixel_snowmelt_shapefile)/plots_per_page))
                  pdf(paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Plot_Pixel_NDSI_Snowmelt_polygon.pdf"), width=20, height=16, onefile = TRUE)
                  for (i in seq(length(plot_pixel_snowmelt_shapefile))) { do.call("grid.arrange", plot_pixel_snowmelt_shapefile[[i]]) }
                  dev.off()
                  }
                
        #(20): Save workspace
         #save.image(paste0(here(), "/Output/MODIS/01_Pixels_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Backup_Workspace_PixelDateOfSnowmelt.RData"))
         
        #(21): Print concluding remarks
         {cat("\n") ; cat("\n")
          print("--------------------------------------------------------------------------------------------------------------------------")
          print(paste0("THE ANALYSIS HAS COMPLETED"))
          cat("\n")
          print(paste0("-Generated data and plots can be found at ", here(), "/Output/MODIS/01_Pixels_Snowmelt/"))
          print("--------------------------------------------------------------------------------------------------------------------------")
          cat("\n")}
         
###########################################################################################################################################################################
###########################################################################################################################################################################
###########################################################################################################################################################################
