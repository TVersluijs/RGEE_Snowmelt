##################################################################################################################################

#This script requires the snow melt maps generated using script "02-RGEE_TomVersluijs_MODIS_Shapefile_Pixel_Snowmelt.R" as input.
#It imports the pixel-level maps for all analyzed years and transforms them into an image with the change in the timing of 
#snowmelt over the years for each pixel (i.e. slope of linear regression) and another image with the average timing of snowmelt 
#over the years for each pixel (i.e. intercept of linear regression).

#Copyright Tom Versluijs 2023-11-01. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

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
      source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "0-RGEE_MODIS_AuxiliaryFunctions_SnowmeltChangeAndIntercept.R")
      sapply(source_files, source, chdir = TRUE) ; rm(source_files)
      
      #(4): Initialize earth engine and google drive
      rgee::ee_Initialize(user='tom.versluijs@gmail.com', drive=TRUE)
      #ee_check()
      

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
     area_name <- "ZAC"

     #Name of Shapefile used to create MODIS annual snowmelt maps
     #Follow the guide "Manual_CreateShapefilePolygons.docx" when creating this shapefile.
     shapefile <- "ZAC_Outline_EPSG4326.shp"
     
     #Coordinate reference system used for calculations
     #EPSG:4326 is recommended for areas spanning multiple UTM zones, but increased computation time (i.e. spherical coordinate system).
     #EPSG:326XX is results in reduced computation time for areas located within a single UTM zone (i.e. planar coordinate system).
     crs <- "EPSG:4326"

   #(c) Dates

     #Specify date and year ranges of all images
     years <- 2000:2023
     start_mmdd <- "-03-15"
     end_mmdd <- "-09-15"

   #(d): Cloud masking

     #Define which MODIS cloud masking algorithm has been used for each year in the analysis ("PGE11", "MOD35", or "Combined)
     MODIS_cloud_masking_algorithm = 'PGE11' #default is PGE11


#################################################################################################################
   
    #Create output folder
    if(dir.exists(paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt"))==FALSE){dir.create(paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt"), recursive = TRUE)}
     
    #Create a timestamp variable
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
     
    #Save all parameters and their values in the environment to a text file 
    file_conn <- file(paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_Parameters.txt"), "w")
    for (obj in setdiff(ls(), lsf.str())) {cat(paste(obj, "=", get(obj)), file = file_conn) ; cat("\n", file = file_conn)}
    close(file_conn)
    
   #(5): Read study area shapefile and convert to a feature collection.
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
      
  #(6): Loop through all MODIS snowmelt images and store all as an imagecollection and a dataframe
      
      #Reload MODIS sinusoidal projection
       start_date <- paste0(years[1], start_mmdd) 
       end_date <- paste0(years[1], end_mmdd)
       start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
       end_date_doy <- as.numeric(strftime(end_date, format = "%j"))
       MODIS_col<-ee$ImageCollection(MODIS_dataset)
       MODIS_col <- MODIS_col$filterBounds(aoi)$filterDate(start_date, end_date)
       modisProjection <- MODIS_col$first()$select("sur_refl_b01")$projection()
       #modisProjection$getInfo()
      
      #Load year-specific MODIS feature collection, convert to an image and save in a list of images and feature collection list
       MODIS_images_list <- list()
       MODIS_FeatureCollection <- ee$FeatureCollection(ee$List(list()))
       MODIS_FeatureCollection_clipped <- ee$FeatureCollection(ee$List(list()))
       for(i in years){
      
        #(A): Load previously generated snowmelt FeatureCollection of year 'i'
          
           #Define a data_ID
            year_ID <- as.character(i)
            data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
            data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)

           #Define the corresponding asset folder
            path_asset <- paste0(ee_get_assethome(), "/", data_ID)
            #ee_manage_assetlist()
            
           #Create overview of filenames in path_asset
            path_asset_files <- ee_manage_assetlist(path_asset)$ID

           #Re-load the already generated snowmelt image as a feature collection
            assetid=path_asset_files[which(grepl("_FC_pixels_snowmelt_optimized", path_asset_files))]
            FC_pixels_MODIS_snowmelt <- ee$FeatureCollection(assetid)
            #FC_pixels_MODIS_snowmelt$first()$getInfo()
         
        #(B): Transform the feature collection of year 'i' into an Image:
          
           #Reduce feature collection FC_pixels_MODIS_snowmelt to an Image with a 500m resolution:
            image_snowmelt <- FC_pixels_MODIS_snowmelt$
             filter(ee$Filter$notNull(list('doy_snowmelt')))$ #pre-filter data for nulls that cannot be turned into an image
             filter(ee$Filter$neq(name='doy_snowmelt', value=-9999))$ #pre-filter data for -9999 values
             reduceToImage(properties=list('doy_snowmelt'), reducer=ee$Reducer$first()$setOutputs(list('doy_snowmelt')))$
             #reproject(crs=crs, crsTransform=NULL, scale=resolution)
             reproject(crs=modisProjection, crsTransform=NULL, scale=resolution)
         
           #Add year as a property and band to each image 
            image_snowmelt <- image_snowmelt$
             set('Year_property', year_ID)$
             addBands(ee$Image$constant(i)$rename('Year'))$
             cast(list('Year'='int'))
         
           #Plot image on a map (for debugging)
            #Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
            #Map$addLayer(image_snowmelt,list(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')

           #Store image_snowmelt in a list
            MODIS_images_list[[which(i==years)]] <- image_snowmelt
         
        #(C): Transform image_snowmelt to a FeatureCollection of sampling points (where points correspond to pixels) at a 500 meter resolution
           
            #Extract Pixel_year and doy_snowmelt at each pixel within the current image and store each sample value as a separate feature
            #The resulting output is a feature collection of all features (sampled pixels) for the current image.
            FC_image <- image_snowmelt$sample( #sampling is done automatically for all Bands of the image
              region=aoi, #All pixels within aoi will be stored as a separate feature
              geometries=TRUE,  #if TRUE, add center of sampled pixel as the geometry property of the output feature
              projection=modisProjection, #Set to native resolution of satellite image (10m)
              scale=resolution, #sampling resolution in meters
              seed=23, #Create reproducable results using the same random seed
              dropNulls=FALSE) #If TRUE, the result is post-filtered to drop features that have a NULL value for all bands
            
            #Make sure there is a NDSI value at each feature within the feature collection (now redundant due to dropNulls=T above):
            #Set the band value to a no data value of -9999 for all features where the band value is NULL.
             FC_image <- FC_image$map(function(feature){
               Year <- ee$List(list(feature$get('Year'), -9999))$reduce(ee$Reducer$firstNonNull())
               doy_snowmelt <- ee$List(list(feature$get('doy_snowmelt'), -9999))$reduce(ee$Reducer$firstNonNull())
               return(feature$set("Year", Year)$set("doy_snowmelt", doy_snowmelt))})
     
            #Add latitude and longitude of each pixel as a property to each feature
            FC_image <- FC_image$map(function(feature){
              coordinates <- feature$geometry()$coordinates()
              lon <- coordinates$get(0)
              lat <- coordinates$get(1)
              return(feature$set('lon', lon)$set('lat', lat))
              })
            
            #Add a pixel_ID property to each feature
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
            
            #Merge the feature collection of the current image (FC_image) onto the feature collection MODIS_FeatureCollection
            MODIS_FeatureCollection <- ee$FeatureCollection(MODIS_FeatureCollection)$merge(FC_image)
           
        #(D): Transform a clipped image_snowmelt to a FeatureCollection of sampling points (where points correspond to pixels) at a 500 meter resolution
            
            #Extract Pixel_year and doy_snowmelt at each pixel within the current image and store each sample value as a separate feature
            #The resulting output is a feature collection of all features (sampled pixels) for the current image.
            FC_image_clipped <- image_snowmelt$clipToCollection(aoi_Shapefile)$sample( #sampling is done automatically for all Bands of the image
              region=aoi_Shapefile, #All pixels within aoi_Shapefile will be stored as a separate feature
              geometries=TRUE,  #if TRUE, add center of sampled pixel as the geometry property of the output feature
              projection=modisProjection, #Set to native resolution of satellite image (10m)
              scale=resolution, #sampling resolution in meters.
              seed=23, #Create reproducable results using the same random seed
              dropNulls=FALSE) #If TRUE, the result is post-filtered to drop features that have a NULL value for all bands
            
            #Make sure there is a NDSI value at each feature within the feature collection (now redundant due to dropNulls=T above):
            #Set the band value to a no data value of -9999 for all features where the band value is NULL.
            FC_image_clipped <- FC_image_clipped$map(function(feature){
              Year <- ee$List(list(feature$get('Year'), -9999))$reduce(ee$Reducer$firstNonNull())
              doy_snowmelt <- ee$List(list(feature$get('doy_snowmelt'), -9999))$reduce(ee$Reducer$firstNonNull())
              return(feature$set("Year", Year)$set("doy_snowmelt", doy_snowmelt))})
            
            #Add latitude and longitude of each pixel as a property to each feature
            FC_image_clipped <- FC_image_clipped$map(function(feature){
              coordinates <- feature$geometry()$coordinates()
              lon <- coordinates$get(0)
              lat <- coordinates$get(1)
              return(feature$set('lon', lon)$set('lat', lat))
            })
            
            #Add a pixel_ID property to each feature
            FC_image_clipped <- FC_image_clipped$map(function(feature){
              #feature=FC_image_clipped$first() #for debugging
              lon_tmp <- feature$get("lon")
              lat_tmp <- feature$get("lat")
              lon_tmp <- ee$Number$format(lon_tmp, '%.5f')
              lat_tmp <- ee$Number$format(lat_tmp, '%.5f')
              string_tmp <- ee$String("_")
              pixel_ID <- lat_tmp$cat(string_tmp)$cat(lon_tmp)
              #feature <- feature$set('pixel_ID', pixel_ID) #for debugging
              return(feature$set('pixel_ID', pixel_ID))
            })  
            
            #Merge the feature collection of the current image (FC_image_clipped) onto the feature collection MODIS_FeatureCollection_clipped
            MODIS_FeatureCollection_clipped <- ee$FeatureCollection(MODIS_FeatureCollection_clipped)$merge(FC_image_clipped)
           
        #(E): Print progress
          print(paste0("Progress: year ", year_ID,  " finished"))
         
      }
      
      #Change data_ID and asset folder to a general MODIS folder (i.e. not specific to a certain year)
       data_ID <- paste0(area_name, "_MODIS_", MODIS_cloud_masking_algorithm)
       path_asset <- paste0(ee_get_assethome(), "/", data_ID)
       tryCatch(ee_manage_delete(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
       ee_manage_create(path_asset=path_asset, asset_type="Folder")
       
     #(A): Store all images of the region 'aoi' as an ImageCollection
       MODIS_images_snowmelt <- ee$ImageCollection$fromImages(MODIS_images_list)
   
       #Specify a function to clip image collection to aoi_Shapefile 
       clip_to_aoi_Outline <-  function(img) {
         return(img$clipToCollection(aoi_Shapefile))
         }
       
       #Create a RGB-timelapse video of (clipped) MODIS snowmelt images
       visFun_Snowmelt <-  function(img) {
        return(img$visualize(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red'))$copyProperties(img, img$propertyNames()))
        }
       MODIS_images_snowmelt_RGB <- MODIS_images_snowmelt$map(clip_to_aoi_Outline)$map(visFun_Snowmelt)
       videoArgs <- list(dimensions=600, region=aoi, framesPerSecond=4, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
       tryCatch({browseURL(MODIS_images_snowmelt_RGB$getVideoThumbURL(videoArgs)) }, error = function(cond){return("Too many pixels. Reduce dimensions.")})
      
     #(B): Transform the feature collection for region 'aoi' to a dataframe where each row contains a pixel_ID, year and date of snow melt:
       
       #create a current timestamp to prevent identical names on Google Drive
       current_timestamp1 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
       
       a=Sys.time()
       task_vector1 <- ee_table_to_drive(
         collection = MODIS_FeatureCollection,
         description = paste0(current_timestamp1, "_", data_ID, "_Pixel_Snowmelt_MODIS_bbox"),
         fileFormat = "CSV",
         selectors = c('doy_snowmelt', 'Year', 'lat', 'lon', 'pixel_ID')
         )
      
       task_vector1$start()
       ee_monitoring(task_vector1, max_attempts = 1000000)
      
       exported_stats <- ee_drive_to_local(task = task_vector1, dsn=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, "_Pixel_Snowmelt_bbox"))
       MODIS_pixel_snowmelt <- read.csv(exported_stats)
       unlink(exported_stats)
       b=Sys.time()
       print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))
      
       #Change -9999 to NA
       MODIS_pixel_snowmelt$doy_snowmelt[MODIS_pixel_snowmelt$doy_snowmelt < -9000] <- NA

       #Sort dataframe by pixel_ID
       MODIS_pixel_snowmelt <- MODIS_pixel_snowmelt[order(MODIS_pixel_snowmelt$pixel_ID, MODIS_pixel_snowmelt$Year),]
   
    #(C): Transform the feature collection for the region 'aoi_Shapefile' to a dataframe where each row contains a pixel_ID, year and date of snowmelt
       
       #create a current timestamp to prevent identical names on Google Drive
       current_timestamp2 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
       
       a=Sys.time()
       task_vector2 <- ee_table_to_drive(
         collection = MODIS_FeatureCollection_clipped,
         description = paste0(current_timestamp2, "_", data_ID, "_Pixel_Snowmelt_MODIS_polygon"),
         fileFormat = "CSV",
         selectors = c('doy_snowmelt', 'Year', 'lat', 'lon', 'pixel_ID')
         )
       
       task_vector2$start()
       ee_monitoring(task_vector2, max_attempts = 1000000)
       
       exported_stats <- ee_drive_to_local(task = task_vector2, dsn=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, "_Pixel_Snowmelt_polygon"))
       MODIS_pixel_snowmelt_clipped <- read.csv(exported_stats)
       b=Sys.time()
       print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))
       
       #Change -9999 to NA
       MODIS_pixel_snowmelt_clipped$doy_snowmelt[MODIS_pixel_snowmelt_clipped$doy_snowmelt < -9000] <- NA
       
       #Sort dataframe by pixel_ID
       MODIS_pixel_snowmelt_clipped <- MODIS_pixel_snowmelt_clipped[order(MODIS_pixel_snowmelt_clipped$pixel_ID, MODIS_pixel_snowmelt_clipped$Year),]
       
      #At this point, MODIS_images_snowmelt is a collection of images for the region 'aoi' (i.e. unclipped) indicating the
      #date of snowmelt for all pixels for all years of interest based on MODIS satellite imagery. MODIS_pixel_snowmelt is
      #a dataframe indicating the date of snowmelt for all pixels in all years within the 'aoi' region (i.e. unclipped). 
      #MODIS_pixel_snowmelt_clipped is a dataframe indicating the date of snowmelt for all pixels in all years within 
      #aoi_Shapefile (i.e. only the subset of pixels that fall within aoi_Shapefile).
         
######################################################################################################################################         
       
  # #(7): (OPTIONAL) Calculate year specific snowmelt parameters (mean and variance in day of snowmelt) based on all pixels within aoi_Shapefile  
  #     MODIS_snowmelt_param <- MODIS_pixel_snowmelt_clipped %>% 
  #                                           group_by(Year) %>% 
  #                                           summarize(mean_snowmelt=mean(na.omit(doy_snowmelt)),
  #                                                     var_snowmelt=var(na.omit(doy_snowmelt)))
  #      
  #    #(7.1):Variance in snowmelt vs Year 
  #      with(MODIS_snowmelt_param, plot(var_snowmelt ~ Year))
  #      lm0 <- with(MODIS_snowmelt_param, lm(var_snowmelt ~ 1)) 
  #      lm1 <- with(MODIS_snowmelt_param, lm(var_snowmelt ~ Year)) 
  #      lm2 <- with(MODIS_snowmelt_param, lm(var_snowmelt ~ Year + I(Year^2))) 
  #      AIC(lm0, lm1, lm2)
  #      #Linear model is best supported. Heterogeneity in snowmelt dates increases over time
  #      
  #      #Fit predictions of lm0 to data
  #      newdata <- data.frame(Year=seq(min(MODIS_snowmelt_param$Year), max(MODIS_snowmelt_param$Year), 1))
  #      newdata$lm_predict <- predict(lm1, newdata=newdata, type="response")
  #      
  #      #Construct a plot with model predictions
  #      p_snowmelt_year_var <- ggplot2::ggplot()+
  #        geom_point(data=MODIS_snowmelt_param, aes(x=Year, y=var_snowmelt))+
  #        geom_line(data=newdata, aes(x=Year, y=lm_predict), col = "red")+
  #        xlab("Year")+
  #        ylab("Variance in date of snowmelt")+
  #        theme_classic()
  #      
  #    #(7.2): Variance in snowmelt vs mean snowmelt  
  #      with(MODIS_snowmelt_param, plot(var_snowmelt ~ mean_snowmelt))
  #      lm0 <- with(MODIS_snowmelt_param, lm(var_snowmelt ~ 1)) 
  #      lm1 <- with(MODIS_snowmelt_param, lm(var_snowmelt ~ mean_snowmelt)) 
  #      lm2 <- with(MODIS_snowmelt_param, lm(var_snowmelt ~ mean_snowmelt + I(mean_snowmelt^2))) 
  #      AIC(lm0, lm1, lm2)
  #      #quadratic model is best supported. Largest heterogeneity in snowmelt in early or late years...
  #     
  #      #Fit predictions of lm0 to data
  #      newdata <- data.frame(mean_snowmelt=seq(min(MODIS_snowmelt_param$mean_snowmelt), max(MODIS_snowmelt_param$mean_snowmelt), 1))
  #      newdata$lm_predict <- predict(lm2, newdata=newdata, type="response")
  #     
  #      #Construct a plot with model predictions
  #      p_snowmelt_mean_var <- ggplot2::ggplot()+
  #       geom_point(data=MODIS_snowmelt_param, aes(x=mean_snowmelt, y=var_snowmelt))+
  #       geom_line(data=newdata, aes(x=mean_snowmelt, y=lm_predict), col = "red")+
  #       xlab("Mean date of snowmelt")+
  #       ylab("Variance in date of snowmelt")+
  #       theme_classic()

            
######################################################################################################################################         

    #(8): Calculate the change in the timing of snow melt per pixel
       
       #(Method I): Use GEE linearFit function on all images in the ImageCollection 'MODIS_images_snowmelt':
        
         #This method uses the linearFit function from google earth engine to fit a linear model through the doy_snowmelt ~ Year
         #data for all pixels in the image collection. This method is easy, but no very flexible regarding the output or in the
         #model fitting options.
         
         #Calculate slope image.
          slope <- MODIS_images_snowmelt$select(list('Year', 'doy_snowmelt'))$reduce(ee$Reducer$linearFit())$select('scale')
          #slope <- slope$reproject(crs=crs, crsTransform=NULL, scale=resolution)
          slope <- slope$reproject(crs=modisProjection, crsTransform=NULL, scale=resolution)
          slope <- slope$clipToCollection(aoi_Shapefile)
          
         #Define visualization arguments.
          visArgs = list(min=-1, max=1, palette= c('green', 'yellow', 'red'))
          
         #Plot on map (green indicate advancement in date of snowmelt, red indicate delay)
          Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
          Map$addLayer(slope, visArgs, paste0("'", years[1], "-", rev(years)[1], " Change in date of snowmelt"))
          
       #(Method II): Manually fit a linear regression through all pixels in the dataframe 'MODIS_pixel_snowmelt':
          
          #Calculate the change in the date of snowmelt and snowmelt intercept for each pixel in the region 'aoi':
          #Loop through all pixel_ID's, select dataframe for that pixel containing Date of snowmelt values for all
          #years (i.e. for all images in the image collection). Fit a linear model through the day of snowmelt ~ year
          #data, and extract the slope and intercept of this fitted model. 
         
          #Specify number of cores:
          numCores <- detectCores()
          
          #Make clusters and initialize using DoSNOW as this allows for a progress bar
          cl <- makePSOCKcluster(numCores)
          registerDoSNOW(cl)
          
          #Split dataset into chuncks of 10,000 pixels (resulting in 1 chunk in this case)
          chunk_size = 2500
          pixel_IDs <- unique(MODIS_pixel_snowmelt$pixel_ID)
          pixel_IDs_split <- split(pixel_IDs, ceiling(seq_along(pixel_IDs)/chunk_size))
          
          #Load the function that iterates through all data subsets. Within each data subset it calculates the change in date 
          #of snowmelt over years for every pixel by fitting a linear model and extracting the slope. The code employs parallel 
          #processing using foreach and %dopar% on numCores local computer cores.
          f_ChangeInSnowmelt_FitLinearModel_parallel_2 <- f_ChangeInSnowmelt_FitLinearModel_parallel_2 #sourced
          
          #Run this datasubset_parallelfunc over all data subsets, combine the results and save the resulting dataframe and plots
          print("Calculate the change in the date of snowmelt for each pixel:")
          results <- lapply(1:length(pixel_IDs_split), FUN=f_ChangeInSnowmelt_FitLinearModel_parallel_2,
                            MODIS_pixel_y=MODIS_pixel_snowmelt, pixel_ID_column="pixel_ID", 
                            y="doy_snowmelt", x="Year")
          
          #Clean up the cluster after finishing the parallel runs
          stopCluster(cl)
          
          #Turn parallel processing off and run sequentially again after this point
          registerDoSEQ()
          
          #Store change in the date of snowmelt for all pixels within the region 'aoi' in a dataframe
          df_pixel_ChangeInSnowmelt <- lapply(results, "[[", 1)
          df_pixel_ChangeInSnowmelt <- as.data.frame(do.call(rbind, do.call(c, df_pixel_ChangeInSnowmelt)))
          colnames(df_pixel_ChangeInSnowmelt)[colnames(df_pixel_ChangeInSnowmelt)=="y_change"] <- "snowmelt_change"
          colnames(df_pixel_ChangeInSnowmelt)[colnames(df_pixel_ChangeInSnowmelt)=="y_intercept"] <- "snowmelt_intercept"
          #write.csv(df_pixel_ChangeInSnowmelt, file=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, "_Pixel_ChangeInSnowmelt_bbox.csv"), quote = FALSE, row.names=FALSE)
          #df_pixel_ChangeInSnowmelt <- read.csv(file=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, "_Pixel_ChangeInSnowmelt.csv"), header=TRUE) 
          
          #The dataframe df_pixel_ChangeInSnowmelt contains the change in the date of snowmelt for each individual pixel_ID within the
          #region 'aoi'. To be able to plot these data we need to transform this dataframe to a feature collection and then transform this 
          #feature collection to an image.
        
         #Calculate some summarizing statistics:    
          
          #Before transforming df_pixel_ChangeInSnowmelt to an image we can first calculate some summarizing statistics to look at changes 
          #in timing of snowmelt. In this case we're specifically interested in all pixels within aoi_Shapefile (and not all pixels within the
          #region 'aoi'). We thus first define which pixels of the current selection fall within aoi_Shapefile and store this as an index variable.
          index <- which(df_pixel_ChangeInSnowmelt$pixel_ID %in% unique(MODIS_pixel_snowmelt_clipped$pixel_ID))
          write.csv(df_pixel_ChangeInSnowmelt[index,], file=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, "_Pixel_ChangeInSnowmelt_polygon.csv"), quote = FALSE, row.names=FALSE)
          
          #Calculate the change in timing of snowmelt within aoi_Shapefile
          hist(na.omit(df_pixel_ChangeInSnowmelt$snowmelt_change[index]), nclass = 30)
          mean(na.omit(df_pixel_ChangeInSnowmelt$snowmelt_change[index]))
          #An average DELAY in the date of snowmelt of 0.011 days per year! But, this is largely due to the very late year 2018.
          
          #Combine all plots of changes in timing of snowmelt for all pixels within aoi_Shapefile in a single list
          plot_pixel_snowmelt <- lapply(results, "[[", 2)
          plot_pixel_snowmelt <- do.call(c, plot_pixel_snowmelt)
          
          #Restrict to only those pixels wihtin aoi_shapefile
          plot_pixel_snowmelt <- plot_pixel_snowmelt[index]

          #Plot the GAM fits with 25 plots per pdf page
          plots_per_page = 25
          plot_pixel_snowmelt <- split(plot_pixel_snowmelt, ceiling(seq_along(plot_pixel_snowmelt)/plots_per_page))
          pdf(paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, "_Pixel_ChangeInSnowmelt_polygon.pdf"), width=20, height=16, onefile = TRUE)
          for (i in seq(length(plot_pixel_snowmelt))) { do.call("grid.arrange", plot_pixel_snowmelt[[i]]) }
          dev.off()
          
         #Transform df_pixel_ChangeInSnowmelt to a feature collection with random geometry
          
          #Generate some random longitude and latitude values as this is required for an sf object:
          df_pixel_ChangeInSnowmelt$lon <- runif(nrow(df_pixel_ChangeInSnowmelt), 0, 74)
          df_pixel_ChangeInSnowmelt$lat <- runif(nrow(df_pixel_ChangeInSnowmelt), -20, 20)
          
          #Prevent NA in snowmelt_change and snowmelt_intercept columns (earthengine cannot deal with NA)  
          df_pixel_ChangeInSnowmelt$snowmelt_change[is.na(df_pixel_ChangeInSnowmelt$snowmelt_change)] <- -9999
          df_pixel_ChangeInSnowmelt$snowmelt_intercept[is.na(df_pixel_ChangeInSnowmelt$snowmelt_intercept)] <- -9999
          
          #memory limits prevent us from processing the whole dataframe at once. Therefore we split it up into subsets of 2500 rows
          chunk_size = 2500
          rowIDs <- 1:nrow(df_pixel_ChangeInSnowmelt)
          rowIDs_split <- split(rowIDs, ceiling(seq_along(rowIDs)/chunk_size))
          
          #Iterate through all dataframe subsets, convert each to a feature collection, and append each to the feature collection FC_initial.
          FC_initial <- ee$FeatureCollection(ee$List(list()))
          for(i in 1:length(rowIDs_split)){ 
            
            #Select subset of rows
            rowID_min <- min(rowIDs_split[[i]])
            rowID_max <- max(rowIDs_split[[i]])
            
            #Select subset of dataframe:
            df_tmp <- df_pixel_ChangeInSnowmelt[rowID_min:rowID_max,]
            
            #Change subset-dataframe to a SF object 
            print("Transform df_pixel_ChangeInSnowmelt to a feature collection with random geometry:")
            df_sf_tmp <- st_as_sf(x = df_tmp,                         
                                  coords = c("lon", "lat"),
                                  crs=crs)
            
            #Change sf object to an earth engine feature collection by uploading it to the asset folder
            FC_tmp <- sf_as_ee(
              x = df_sf_tmp,
              assetId = paste0(path_asset, "/", timestamp, "_", data_ID, "_df_pixel_ChangeInSnowmelt_", i),
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
          FC_pixels_ChangeInSnowmelt <- ee$FeatureCollection(FC_initial)
          FC_pixels_ChangeInSnowmelt$first()$getInfo()
          FC_pixels_ChangeInSnowmelt$size()$getInfo()
          
          #Inspect assets folder:
          ee_manage_quota()
          ee_manage_assetlist(path_asset)
          
        #Add geometry (latitude and longitude) of each pixel_ID to FC_pixels_ChangeInSnowmelt
          
          #The feature collection that we want to construct should contain 'pixel_ID' and 'snowmelt_change' as properties and should contain
          #the original geometry corresponding to each pixel_ID. So far, FC_pixels_ChangeInSnowmelt contains a separate feature for each pixel, 
          #where each feature contains snowmelt_change and pixel_ID as a property. However, the geometry (lan/lon) of each feature is randomly 
          #chosen. We need to make sure that the actual geometry matching each pixel_ID is added instead of this random geometry.
          
          #To obtain the corresponding geometry (lat/lon) of every pixel_ID within aoi_Shapefile, we sample from a single MODIS image on a 'resolution' 
          #resolution using img$sample(). This gives a single distinct feature for each pixel, including their geometry. We can then re-construct
          #the property pixel_ID for every feature in this feature collection. The resulting feature collection is called FC_pixels_distinct. 
          #The final step is then to join FC_pixels_ChangeInSnowmelt to FC_pixels_distinct based on an inner join with pixel_ID.
          
          #This step results in memory errors (on the server side) when the image contains more than 2.5 million pixels. In the latter case
          #we need to split up our study area in several sub-areas (i.e. create separate adjacent shapefiles using QGIS).
          
          #(A): Select a single image from the image collection
          img <- MODIS_images_snowmelt$first()$select("doy_snowmelt")
          
          #(B): Create a feature collection of points at a 'resolution' resolution within aoi_Shapefile (each point thus corresponds to a pixel in the image)
          FC_pixels_distinct <- img$sample(
            region=aoi, #Sample all pixels within aoi_Shapefile
            geometries=TRUE,  #Store lat/lon in geometry property
            projection=modisProjection, #Set to native projection of MODIS image
            scale=resolution, #sampling resolution in meters
            seed=23, #set fixed seed to be able to reproduce the sample (same as for FC_merged above)
            dropNulls=FALSE) #Make sure NULL pixels are not dropped from image collection
          
          #(C): Make sure there is a doy_snowmelt value at each feature within the feature collection:
          #Set the band value to a no data value of -9999 for all features where the band value is NULL.
          FC_pixels_distinct <- FC_pixels_distinct$map(function(feature){
            doy_snowmelt <- ee$List(list(feature$get('doy_snowmelt'), -9999))$reduce(ee$Reducer$firstNonNull())
            return(feature$set("doy_snowmelt", doy_snowmelt))})
          
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
            #and contains the original lat/lon geometry of every pixel. We will now add snowmelt_change as a new property
            #to FC_pixels_distinct by merging it with FC_pixels_ChangeInSnowmelt.
          
          #(G): Add snowmelt_change as property to FC_pixels_distinct   
          
            #Use an equals filter to specify how the collections match.
            Filter_pixel_ID <- ee$Filter$equals(
              leftField= 'pixel_ID',
              rightField= 'pixel_ID')
          
            #Define the join.
            innerJoin <- ee$Join$inner()
          
            #Apply the join.
            Pixel_ID_Join <- innerJoin$apply(FC_pixels_distinct, FC_pixels_ChangeInSnowmelt,  Filter_pixel_ID)
          
            #Add features of second ('secondary') feature collection to those of the first ('primary') feature collection 
            FC_Combined <- Pixel_ID_Join$map(function(pair) {
              f1 <- ee$Feature(pair$get('primary'))
              f2 <- ee$Feature(pair$get('secondary'))
              return(f1$set(f2$toDictionary()))
            })  
          
          #(H): To speed-up further computations with this feature collection we as an intermediate step upload and re-download it
          
            #create a current timestamp to prevent identical names on Google Drive
            current_timestamp3 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
            
            #Delete FC_pixels_ChangeInSnowmelt_optimized if it already occured in the asset folder:
            tryCatch({ee_manage_delete(paste0(path_asset, "/", current_timestamp3, "_", data_ID, "_FC_pixels_ChangeInSnowmelt_optimized"))}, 
                      error = function(cond){return("Path did not yet exist - no folder deleted")})
          
            #Upload to asset folder:
            assetid2 <- paste0(path_asset, "/", current_timestamp3, "_", data_ID, "_FC_pixels_ChangeInSnowmelt_optimized")
            task_vector3 <- ee_table_to_asset(
              collection = FC_Combined,
              overwrite = TRUE,
              assetId = assetid2
              )
            task_vector3$start()
            print("Optimize further calculations with FC_Combined by uploading it to the asset folder:")
            ee_monitoring(task_vector3, max_attempts = 1000000)
          
          #Check assets folder:
          #ee_manage_quota()
          #ee_manage_assetlist(path_asset)
          
          #Get feature collection from asset folder and create FC_pixels_ChangeInSnowmelt_optimized
          #assetid2=paste0(path_asset, "/", current_timestamp3, "_", data_ID, "_FC_pixels_ChangeInSnowmelt_optimized")
          FC_pixels_ChangeInSnowmelt_optimized <- ee$FeatureCollection(assetid2) 
          #FC_pixels_ChangeInSnowmelt_optimized$first()$getInfo()
          #FC_pixels_ChangeInSnowmelt_optimized$size()$getInfo()
          
############################################################################################################################################################        
       
      #(9.1): Transform the Feature collection FC_pixels_ChangeInSnowmelt_optimized to an image (with snowmelt_change as an image band)
          
          #(A): Reduce feature collection FC_pixels_ChangeInSnowmelt_optimized to an Image with a 500m resolution:    
          image_snowmelt_slope <- FC_pixels_ChangeInSnowmelt_optimized$
            filter(ee$Filter$notNull(list('snowmelt_change')))$ #pre-filter data for nulls that cannot be turned into an image
            filter(ee$Filter$neq(name='snowmelt_change', value=-9999))$ #pre-filter data for -9999 values
            reduceToImage(properties=list('snowmelt_change'), reducer=ee$Reducer$first()$setOutputs(list('snowmelt_change')))$
            #reproject(crs=crs, crsTransform=NULL, scale=resolution)
            reproject(crs=modisProjection, crsTransform=NULL, scale=resolution)
          
            #clip image by aoi_Shapefile
            image_snowmelt_slope <- image_snowmelt_slope$clipToCollection(aoi_Shapefile)
         
          #(B): Extract change in day of snowmelt in year of interest for a single point
          ee_extract(x=image_snowmelt_slope, y=coordinates_point, fun=ee$Reducer$first(), scale=resolution, sf=TRUE)
          
          #(C): Plot snowmelt day of year as a coloured image
          Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
          #Map$addLayer(MODIS_col$filterDate(paste0(year_ID, "-07-30"), end_date)$first(),list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
          Map$addLayer(image_snowmelt_slope,list(bands="snowmelt_change", min=-1, max=1, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')
          
          #(D): Export change in snowmelt image to Google Drive (takes c.a. 2 minutes):
          
            #create a current timestamp to prevent identical names on Google Drive
            current_timestamp4 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
          
            #Create task to export the original snowmelt_change image to Google Drive  
            task_vector4 <- ee_image_to_drive(
              image=image_snowmelt_slope,
              description= paste0(current_timestamp4, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Slope'),
              #scale= resolution,
              region=aoi,
              crs="EPSG:3857",
              maxPixels=1e9,
              fileFormat='GeoTIFF',
              dimensions=ee$Number(1024) #maximum dimension
              )
          
            #Start and monitor export task:
            task_vector4$start()
            print("Export original image to Google Drive:")
            ee_monitoring(task_vector4, max_attempts = 1000000)
            ee_drive_to_local(task = task_vector4, dsn=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Slope'))
            
          #(E): Export RGB image to Google Drive (takes c.a. 2 minutes):
          
            #Convert original image to an RGB image:
            image_snowmelt_slope_RGB <- image_snowmelt_slope$visualize(bands=c('snowmelt_change'), min=-1, max=1, palette=c('green', 'yellow', 'red'))
            #ee_print(image_snowmelt_slope_RGB)
          
            #create a current timestamp to prevent identical names on Google Drive
            current_timestamp5 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
            
            #Create task to export RGB image to Google Drive:
            task_vector5 <- ee_image_to_drive(
              image=image_snowmelt_slope_RGB,
              description= paste0(current_timestamp5, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Slope_RGB'),
              #scale= resolution,
              region=aoi,
              crs="EPSG:3857",
              maxPixels=1e9,
              fileFormat='GeoTIFF',
              dimensions=ee$Number(1024) #maximum dimension
              )
          
            #Start and monitor export task:
            print("Export RGB image to Google Drive:")
            task_vector5$start()
            ee_monitoring(task_vector5, max_attempts = 1000000)
            ee_drive_to_local(task = task_vector5, dsn=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Slope_RGB'))
            
       #(9.2): Transform the Feature collection FC_pixels_ChangeInSnowmelt_optimized to an image (with snowmelt_intercept as an image band)
            
            #(A): Reduce feature collection FC_pixels_ChangeInSnowmelt_optimized to an Image with a 500m resolution:    
            image_snowmelt_intercept <- FC_pixels_ChangeInSnowmelt_optimized$
              filter(ee$Filter$notNull(list('snowmelt_intercept')))$ #pre-filter data for nulls that cannot be turned into an image
              filter(ee$Filter$neq(name='snowmelt_intercept', value=-9999))$ #pre-filter data for -9999 values
              reduceToImage(properties=list('snowmelt_intercept'), reducer=ee$Reducer$first()$setOutputs(list('snowmelt_intercept')))$
              #reproject(crs=crs, crsTransform=NULL, scale=resolution)
              reproject(crs=modisProjection, crsTransform=NULL, scale=resolution)
            
            #clip image by aoi_Shapefile
            image_snowmelt_intercept <- image_snowmelt_intercept$clipToCollection(aoi_Shapefile)
            
            #(B): Extract change in day of snowmelt in year of interest for a single point
            ee_extract(x=image_snowmelt_intercept, y=coordinates_point, fun=ee$Reducer$first(), scale=resolution, sf=TRUE)
            
            #(C): Plot snowmelt day of year as a coloured image
            Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
            #Map$addLayer(MODIS_col$filterDate(paste0(year_ID, "-07-30"), end_date)$first(),list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
            Map$addLayer(image_snowmelt_intercept,list(bands="snowmelt_intercept", min=start_date_doy+50, max=end_date_doy-50, palette=c('green', 'yellow', 'red')), 'snowmelt_intercept')
            
            #(D): Export change in snowmelt image to Google Drive (takes c.a. 2 minutes):
            
              #create a current timestamp to prevent identical names on Google Drive
              current_timestamp6 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
            
             #Create task to export the original snowmelt_change image to Google Drive  
              task_vector6 <- ee_image_to_drive(
                image=image_snowmelt_intercept,
                description= paste0(current_timestamp6, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Intercept'),
                #scale= resolution,
                region=aoi,
                crs="EPSG:3857",
                maxPixels=1e9,
                fileFormat='GeoTIFF',
                dimensions=ee$Number(1024) #maximum dimension
                )
            
             #Start and monitor export task:
              task_vector6$start()
              print("Export original image to Google Drive:")
              ee_monitoring(task_vector6, max_attempts = 1000000)
              ee_drive_to_local(task = task_vector6, dsn=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Intercept'))
              
            #(E): Export RGB image to Google Drive (takes c.a. 2 minutes):
            
             #Convert original image to an RGB image:
              image_snowmelt_intercept_RGB <- image_snowmelt_intercept$visualize(bands=c('snowmelt_intercept'), min=start_date_doy+50, max=end_date_doy-50, palette=c('green', 'yellow', 'red'))
              #ee_print(image_snowmelt_intercept_RGB)
            
             #create a current timestamp to prevent identical names on Google Drive
              current_timestamp7 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
              
             #Create task to export RGB image to Google Drive:
              task_vector7 <- ee_image_to_drive(
                image=image_snowmelt_intercept_RGB,
                description= paste0(current_timestamp7, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Intercept_RGB'),
                #scale= resolution,
                region=aoi,
                crs="EPSG:3857",
                maxPixels=1e9,
                fileFormat='GeoTIFF',
                dimensions=ee$Number(1024) #maximum dimension
                )
            
              #Start and monitor export task:
              print("Export RGB image to Google Drive:")
              task_vector7$start()
              ee_monitoring(task_vector7, max_attempts = 1000000)
              ee_drive_to_local(task=task_vector7, dsn=paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, '_Pixel_Image_SnowmeltDoy_Intercept_RGB'))
              
            #Save workspace
            #save.image(paste0(here(), "/Output/MODIS/03_Shapefile_Pixel_ChangeInSnowmelt/", timestamp, "_", data_ID, "_Backup_Workspace_PixelChangeInSnowmeltDoy.RData"))        
            
               
##################################################################################################################################################
##################################################################################################################################################
            
          