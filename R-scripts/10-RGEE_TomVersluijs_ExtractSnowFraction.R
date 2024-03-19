##################################################################################################################################

#Extract timeseries of the fraction of snowcover from a pixel-level snow melt map for a set of input locations. This corresponds 
#to the method 'pixel_gam' in the other scripts. Input locations can either be point locations with a corresponding buffer zone, 
#or a collection of polygons in a shapefile.

#This script depends on a snowmelt map generated for MODIS using script "02-RGEE_TomVersluijs_#MODIS_Shapefile_Pixel_Snowmelt.R", 
#or for Sentinel-2 using script "08-RGEE_TomVersluijs_S2_Shapefile_Pixel_Snowmelt.R". Please run either script before continuing 
#with the analysis below.

#Copyright Tom Versluijs 2023-11-07. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

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
     p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr, stars, foreach, parallel, doSNOW, gridExtra, rnaturalearth, rnaturalearthdata)
     devtools::install_github("ropensci/rnaturalearthhires")
  
    #(2): Define ggplot2 plotting theme
     theme_tom <- function(){
       theme_classic() %+replace%
       theme(axis.title = element_text(size=18),
           axis.text = element_text(size=16),
           legend.position = "none",
           strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2)),
           complete = TRUE)}
  
    #(3): Initialize earth engine and google drive
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

   #(a): Snowmelt map
     
      #NOTE: Check script 02, or script 08 for the correct parameter settings and values!
       
      #Which satellite was used to generate the snow melt map
       satellite="Sentinel2" #either "MODIS" or "Sentinel2"
       #Note that accurate estimates of timeseries of the fraction of snowcover for MODIS (i.e. 500m pixel resolution) are 
       #only possible when the areas over which this is calculated are at least c.a. > 2km x 2km (i.e. 16 pixels).
     
      #Specify at which resolution the snow melt map was generated (in meters)
       resolution=20
       
      #Specify what NDSI threshold was used to generate the snow melt map 
       NDSI_threshold=0.4 #only a single value allowed 
       
      #What coordinate reference system was used to generate the snow melt map
       crs <- "EPSG:32627"
       
      #Specify the name of the study area
       area_name="ZAC"
       
      #For which year was the snow melt map generated
       year_ID <- "2022"
       
      #What were the min and max date ranges for which snow melt was calculated
       start_date <- paste0(year_ID, "-03-15") 
       end_date <- paste0(year_ID, "-09-15")
  
      #Set which MODIS cloud filtering algorithm was used when generating the snowmelt map (only applicable when satellite="MODIS")
       MODIS_cloud_masking_algorithm = 'PGE11' 
       
  #(b) Specify at which locations the fraction of snowmelt will be extracted from the snow melt map
       
      #Location type
       Location_type="shapefile" #either "points" or "shapefile"
       
      #Name of file with Locations of interest:
       input_locations <- "ZAC_TenEqualSizeVoronoiPolygons_EPSG4326.shp" #"TestLocations_Zackenberg.csv", #"ZAC_TenEqualSizeVoronoiPolygons_EPSG4326.shp", "ZAC_Outline_EPSG4326.shp"
    
       #When Location_type="points":
         #-Make sure that 'input_locations' has the columns "LON_x", "LAT_y"
         #-The "xxxx.csv" file should be placed in the '/Input' folder
      
       #When Location_type="shapefile":
         #-Make sure that each polygon has a separate 'LocationID'
         #-The "xxxx.shp" file should be placed in the '/Input/Shapefiles' folder
       
      #Buffer radius around locations (in meters, only applicable when Location_type="points")
       Buffer_radius_m=250
       
  #(c) Thresholds:
       
      #Define the snowcover fraction for which the date of its occurrence will be calculated (specify multiple using c())
       Snowfraction_threshold_vector = seq(0.25, 0.75, by=0.25)
       
      #GAM sequential outlier filtering
       
       #parameters for auxiliary function to filter outliers from GAM using a two-step approach
        outlier_removal=FALSE #should sequential outlier removal be employed when fitting GAMs to the data
        outlier_thresh_1=0.4 #first threshold in outlier removal (relative to y-range of data, default=0.4)
        outlier_thresh_2=0.2 #second threshold in outlier removal (relative to y-range of data, default=0.2)
       
       #Specify the degree of GAM smoothing by setting the 'k' parameter
        gam_k_outlier=10 #Number of knots when filtering outliers (default=10)
        gam_k=50 #Number of knots when making model predictions (default=50). 
        #Larger values result in a more precise GAM-fit, but might result in overfitting.

      
##################################################################################################################################

#III: Define some additional parameters (automated)

##################################################################################################################################

  #Load auxiliary functions
   source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = paste0(satellite, "_AuxiliaryFunctions"))
   sapply(source_files, source, chdir = TRUE) ; rm(source_files)
   
  #Specify some satellite specific parameters
   if(satellite=="MODIS"){
     
     #Specify the folder where the snowmelt image is stored:
     dir_SnowmeltImage <- paste0(here(), "/Output/MODIS/02_Shapefile_Pixel_Snowmelt/")
     
     #Create a unique data_ID:
     data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
     data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)
     
     #Specify and create output folder
     dir_Output <- paste0(here(), "/Output/MODIS/10_Extract_SnowFraction")
     if(dir.exists(dir_Output)==FALSE){dir.create(dir_Output, recursive = TRUE)}
    
     }
   if(satellite=="Sentinel2"){
     
     #Specify the folder where the snowmelt image is stored:
     dir_SnowmeltImage <- paste0(here(), "/Output/S2/08_Shapefile_Pixel_Snowmelt/")
     
     #Create a unique data_ID:
     data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_S2")
     
     #Specify and create output folder
     dir_Output <- paste0(here(), "/Output/S2/10_Extract_SnowFraction")
     if(dir.exists(dir_Output)==FALSE){dir.create(dir_Output, recursive = TRUE)}
     
     }
   
  #Create a timestamp variable
   timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
   
  #Store NDSI_threshold as a character (used for naming of outputs)
   NDSI_threshold_char <- gsub("\\.", "_", as.character(NDSI_threshold))  
   
  #Specify starting and ending date as day of year 
   start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
   end_date_doy <- as.numeric(strftime(end_date, format = "%j"))      
   
  #Load the path for asset folder that contains the snow melt image
   
     #Specify asset path for 'data_ID'
      path_asset <- paste0(ee_get_assethome(), "/", data_ID)
     
     #Load the name of the assetid .Rds file from the output of script 02 (MODIS) or script 08 (Sentinel2)
      assetid_file <- list.files(path=dir_SnowmeltImage, full.names=F, recursive = T, pattern = paste0(data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Variable_AssetID.Rds"))
      #list.files(path=dir_SnowmeltImage, full.names=F, recursive = T, pattern = paste0("_Variable_AssetID.Rds"))
     
     #Read assetid file and store the assetid path in the variable assetid1
      assetid1 <- readRDS(file=paste0(paste0(dir_SnowmeltImage, assetid_file)))
     
     #Check if this assetid path exists in the 'path_asset' asset folder
      if(!assetid1 %in% ee_manage_assetlist(path_asset)$ID){print("AssetID not found in asset folder. Inspect ee_manage_assetlist(path_asset) to find the correct assetID.")}

   #Check if buffer size is large enough relative to resolution
    if(Location_type=="points" & (Buffer_radius_m < 2 * resolution)){print("ERROR: Buffer size too small relative to specified resolution of image. Increase Buffer_radius_m for better estimates of fraction of snowcover.")}
      
      
##################################################################################################################################

#IV: Read dataframe with Locations of interest

##################################################################################################################################

  #(1): For point locations:    
    if(Location_type=="points"){
      
      #Load dataframe with point locations
      df_locations <- read.csv(paste0(here(), "/Input/", input_locations), header=T)[,c("LON_x", "LAT_y")]
      
      #Add a unique LocationID to every unique lat/lon combination
      df_locations <- unique(df_locations[,c("LON_x", "LAT_y")])
      df_locations$LocationID <- as.factor(1:nrow(df_locations))
      
      #Transform to an ee object
      Locations_sf <- st_as_sf(df_locations, coords = c("LON_x", "LAT_y"), crs = "epsg:4326")
      Locations_ee <- st_transform(Locations_sf, crs="epsg:4326")
      Locations_ee <- sf_as_ee(Locations_ee[,c("LocationID", "geometry")])
      
      #Calculate mean location
      coordinates_point <- Locations_ee$geometry()$centroid()$getInfo()$coordinates
      
      #Add buffer around all point locations
      if(Buffer_radius_m>0){
        bufferBy <- function(Buffer_radius_m) {
          return(function(feature) {
            return(feature$buffer(Buffer_radius_m))   
          })
        }
        Locations_ee <- Locations_ee$map(bufferBy(Buffer_radius_m))
      }
      
      #Locations is now a feature collection. A feature in our analysis comprises the feature type (point in case no buffer
      #was added, or a circular polygon with radius Buffer_radius_m with a center the original point location), Location_ID  
      #and the feature's geometry as lat and lon coordinates. A feature collection corresponds to a collection of such features 
      #(i.e. a collection of different Locations, with their respective feature types, LocationIDs and coordinates).
      
     #Plot the point locations including buffer zone
      
       #(1): Using ggplot
        
         #Extract bounding box around all points
         Locations_sf <- ee_as_sf(Locations_ee)
         Locations_bbox <- st_bbox(Locations_sf) #bounding box surrounding all point locations
         zoom_level=1.5 #specify by how much the map needs to be zoomed out
        
         #Expand boundingbox by zoom_level * data_range
         x_range <- Locations_bbox$xmax - Locations_bbox$xmin
         y_range <- Locations_bbox$ymax - Locations_bbox$ymin
         Locations_bbox$xmin <- Locations_bbox$xmin - zoom_level * x_range
         Locations_bbox$ymin <- Locations_bbox$ymin - zoom_level * y_range
         Locations_bbox$xmax <- Locations_bbox$xmax + zoom_level * x_range
         Locations_bbox$ymax <- Locations_bbox$ymax + zoom_level * y_range
        
         #Download world map
         world_map <- ne_countries(scale = "large", returnclass = "sf")
         
         #Create plot of point locations on world map
         p_points <- ggplot() +
           geom_sf(data=world_map, fill= "antiquewhite")+
           geom_sf(data=Locations_sf , fill="red", col = "black")+ #fill=sf.colors(nrow(Locations_sf))
           geom_sf_text(data=Locations_sf, aes(label=LocationID), colour="white")+
           theme_tom()+ theme(panel.background = element_rect(fill = "aliceblue"))+
           coord_sf(xlim=c(Locations_bbox$xmin, Locations_bbox$xmax) , ylim=c(Locations_bbox$ymin, Locations_bbox$ymax))
        
         ggsave(plot=p_points, paste0(dir_Output, "/", timestamp, "_", data_ID, "_Locations_Points.pdf"), width=10, height=8)
        
       #(2): Using Earth Engine
         
          #Plot feature collection on a map (i.e. locations as points)
          Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
          Map$addLayer(Locations_ee, list(color="red"), paste0("Locations_", year_ID))
          
       #(3): Create an affix for output naming
         output_affix <- paste0("_Buffer", Buffer_radius_m, "_Res", resolution)
        
    }
      
  #(2): For shapefiles:
    if(Location_type=="shapefile"){
      
      #Read shapefiles
      df_locations <- st_read(paste0(here(), "/Input/Shapefiles/", input_locations), quiet=T)
      colnames(df_locations) <- c("LocationID", "geometry")
      
      #Set column LocationID to a factor
      df_locations$LocationID <- as.factor(df_locations$LocationID)
      
      #Plot the locations:
      p_shapefile <- ggplot() + 
        geom_sf(data = df_locations, fill=sf.colors(nrow(df_locations)), col = "black")+
        geom_sf_label(data = df_locations, aes(label=LocationID), colour="black")+
        theme_tom()
      
      ggsave(plot=p_shapefile, paste0(dir_Output, "/", timestamp, "_", data_ID, "_Locations_Shapefile.pdf"), width=10, height=8)
      
      #Convert shapefile to an earthengine feature collection:  
      Locations_ee <- st_transform(df_locations, crs="EPSG:4326")
      Locations_ee <- sf_as_ee(Locations_ee)
      Locations_ee <- ee$FeatureCollection(Locations_ee)
      
      #Calculate mean location
      coordinates_point <- Locations_ee$geometry()$centroid()$getInfo()$coordinates
      
      #Plot featurecollection
      Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
      Map$addLayer(Locations_ee, list(color="red"), paste0("Locations_", year_ID))
      
      #Create an affix for output naming
      output_affix <- paste0("_Res", resolution)
      
    }

   
##################################################################################################################################

#IV: Recreate generated snow melt image

##################################################################################################################################

  #(1): Download 'FC_pixels_snowmelt_optimized' from the asset folder and transform to feature collection
   FC_pixels_snowmelt_optimized <- ee$FeatureCollection(assetid1)

  #(2): Reduce feature collection FC_pixels_snowmelt_optimized to a snow melt Image:
   image_snowmelt <- FC_pixels_snowmelt_optimized$
     filter(ee$Filter$notNull(list('doy_snowmelt')))$ #pre-filter data for nulls that cannot be turned into an image
     filter(ee$Filter$neq(name='doy_snowmelt', value=-9999))$ #pre-filter data for -9999 values
     reduceToImage(properties=list('doy_snowmelt'), reducer=ee$Reducer$first()$setOutputs(list('doy_snowmelt')))$
     reproject(crs=crs, crsTransform=NULL, scale=resolution) #ensures sure that reduceToImage above is done in the S2 projection
  
  #(3): Plot snowmelt day of year as a coloured image
   Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
   Map$addLayer(image_snowmelt,list(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')+
   Map$addLayer(Locations_ee, list(color="red"), paste0("Locations_", year_ID))

  #(4): extract the average day of snowmelt for all locations
   ee_extract(x=image_snowmelt, y=Locations_ee, fun=ee$Reducer$mean(), scale=resolution, sf=TRUE)
   

##################################################################################################################################

#V: Extract the date of snow melt for all pixels for each location

##################################################################################################################################

#(1): Extract the date of snowmelt for all pixels for each location

  #(A): Extract the date of snowmelt at each pixel for all locations in 'Locations_ee':
  FC_snowmelt_locations <- image_snowmelt$sampleRegions( #sampling is done automatically for all Bands of the input image
    collection=Locations_ee, #All pixels within aoi_Shapefile will be stored as a separate feature
    scale=resolution, #sampling resolution in meters
    #projection=S2Projection, #Set to native projection of S2 image
    geometries=TRUE)  #if TRUE, add center of sampled pixel as the geometry property of the output feature
  
  #(B): Make sure there is a doy_snowmelt value at each feature within the feature collection:
  #Set the band value to a no data value of -9999 for all features where the band value is NULL.
  FC_snowmelt_locations <- FC_snowmelt_locations$map(function(feature){
    doy_snowmelt <- ee$List(list(feature$get('doy_snowmelt'), -9999))$reduce(ee$Reducer$firstNonNull())
    return(feature$set("doy_snowmelt", doy_snowmelt))})
  
  #(C): Add latitude and longitude of each pixel as a property to each feature
  FC_snowmelt_locations <- FC_snowmelt_locations$map(function(feature){
    coordinates <- feature$geometry()$coordinates()
    lon <- coordinates$get(0)
    lat <- coordinates$get(1)
    return(feature$set('lon', lon)$set('lat', lat))
    })
  
  #(D): Transform feature collection FC_snowmelt_locations to a dataframe:

    #Use ee_table_to_drive() to prevent memory limits
     task_vector5 <- ee_table_to_drive(
      collection = FC_snowmelt_locations,
      description = paste0(timestamp, "_", data_ID, output_affix, "_Locations_Data_Pixel_Snowmelt"),
      fileFormat = "CSV",
      selectors = c('LocationID', 'lat', 'lon', 'doy_snowmelt')
      )
    
    #Execute task
     task_vector5$start()
     cat("\n")
     print("  -Extract the date of snowmelt from the snowmelt image for all pixels within each Location:")
     ee_monitoring(task_vector5, quiet=FALSE, max_attempts = 1000000)
    
    #Copy file from Google Drive to R-environment
     exported_stats <- ee_drive_to_local(task = task_vector5, dsn=paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Locations_Data_Pixel_Snowmelt"))
     df_locations_pixel_snowmelt <- read.csv(exported_stats)
     unlink(exported_stats)

    #Adjust 'df_locations_pixel_snowmelt'
    
      #Set LocationID to a factor
       df_locations_pixel_snowmelt$LocationID <- as.factor(df_locations_pixel_snowmelt$LocationID)
     
      #Make sure each latitude/longitude combination gets its own pixel_ID
       df_locations_pixel_snowmelt$pixel_ID <- paste0(format(df_locations_pixel_snowmelt$lat, nsmall = 5), "_", format(df_locations_pixel_snowmelt$lon, nsmall = 5))
      
      #Add NDSI_threshold as a new column
       df_locations_pixel_snowmelt$NDSI_threshold <- as.factor(NDSI_threshold)
      
      #Save dataframe
       write.csv(df_locations_pixel_snowmelt, file=paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Locations_Data_Pixel_Snowmelt.csv"), quote = FALSE, row.names=FALSE)


##################################################################################################################################

#VI: Extract timeseries of the fraction of snowcovered pixels for each Location

##################################################################################################################################
       
#(1): Calculate the fraction of snowcovered pixels per location per day

    #Create empty dataframe to store the fraction of snowcover per day of year for each location
    df_Locations_SnowFraction <- data.frame(LocationID=character(),
                                            doy=numeric(),
                                            SnowFraction=numeric(),
                                            NDSI_threshold=numeric(),
                                            LON_x=numeric(),
                                            LAT_x=numeric())     
       
    #Calculate the fraction of snowcover over time per location
    for(Location_i in unique(df_locations_pixel_snowmelt$LocationID)){
        
        ##for debugging
        #Location_i <- unique(df_locations_pixel_snowmelt$LocationID)[1]
        
        #Select dataset for current location
        df_Location_pixel_snowmelt <- df_locations_pixel_snowmelt[df_locations_pixel_snowmelt$LocationID==Location_i,]
        
        #Remove NAs from df_Location_pixel_snowmelt
        df_Location_pixel_snowmelt <- df_Location_pixel_snowmelt[df_Location_pixel_snowmelt$doy_snowmelt>0,]
      
        #Create empty dataframe to store the fraction of snowcover per day of year for a single location
        df_Location_SnowFraction_new <- data.frame(LocationID=character(length(start_date_doy:end_date_doy)),
                                                   doy=start_date_doy:end_date_doy,
                                                   SnowFraction=numeric(length(start_date_doy:end_date_doy)))
        
        #Calculate the fraction of snowcovered pixels for each doy between 'start_date_doy' and 'end_date_doy'
        for(i in 1:nrow(df_Location_SnowFraction_new)){
          
            #Select doy
            doy_i <- df_Location_SnowFraction_new[i,"doy"]
            
            #Count fraction of pixels in df_Location_pixel_snowmelt which are still snowcovered on this doy
            doy_i_SnowcoverFraction <- length(which(df_Location_pixel_snowmelt$doy_snowmelt > doy_i)) / length(df_Location_pixel_snowmelt$doy_snowmelt)
            
            #Store the snowfraction at this doy
            df_Location_SnowFraction_new[i,"SnowFraction"] <- doy_i_SnowcoverFraction
            
            }
        
        #Add LocationID
        df_Location_SnowFraction_new$LocationID <- Location_i
        
        #Add NDSI_threshold as a new column
        df_Location_SnowFraction_new$NDSI_threshold <- as.factor(NDSI_threshold)
        
        #Add coordinates to each location
        if(Location_type=="points"){df_Location_SnowFraction_new <- left_join(df_Location_SnowFraction_new, df_locations, by=c("LocationID"))}
        
        #Add dataframe for current location to dataframe from previous iterations:
        df_Locations_SnowFraction <- rbind(df_Locations_SnowFraction, df_Location_SnowFraction_new)
        
        # #Save dataframe for current location
        # write.csv(df_Location_SnowFraction_new, paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Location_", Location_i, "_Data_SnowFraction.csv"), row.names = FALSE)

      }
       
    #Save snowfraction data per location
    write.csv(df_Locations_SnowFraction, paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Locations_Data_SnowFraction.csv"), row.names = FALSE)
    

##################################################################################################################################

#VII: Fit a Generalized Additive Model (GAM) through the fraction of snowcovered pixels for each Location

##################################################################################################################################

    #(A) Specify for which snowfraction the corresponding date needs to be extracted
    Snowfraction_threshold_vector = Snowfraction_threshold_vector
    
    #(B) Create empty dataframes for storing GAM fits
    df_Locations_SnowFraction_GAM <- data.frame(SnowFraction=numeric(),
                                                Date=factor(),
                                                doy=numeric(),
                                                LocationID=factor(),
                                                outliers=logical())
    
    df_Locations_SnowFraction_GAM_predictions <- data.frame(LocationID=character(),
                                                            doy=numeric(),
                                                            SnowFraction_gam_predict=numeric(),
                                                            stringsAsFactors=FALSE)
    
    #(C) Loop through all Locations and fit a separate gam with sequential outlier removal to the location specific pixel-level SnowFraction data
    for(i in unique(df_Locations_SnowFraction$LocationID)){
      
      #For debugging
      #i=unique(df_Locations_SnowFraction$LocationID)[1]
  
        #Select Location-specific subset of data:
        df_Location_SnowFraction_GAM <- df_Locations_SnowFraction[df_Locations_SnowFraction$LocationID==i & 
                                                              !is.na(df_Locations_SnowFraction$SnowFraction),
                                                              c("SnowFraction", "doy", "LocationID")]
        
        #Fit a gam through the Location-specific mean_NDSI ~ doy data and employ sequential outlier removal
        df_Location_SnowFraction_GAM <- f_gam_SeqRemOutliers(data=df_Location_SnowFraction_GAM, y="SnowFraction", x="doy", outlier_removal=outlier_removal,
                                                               outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                               default_k=gam_k_outlier)
        
        #Sort df_Location_SnowFraction_GAM by doy:
        df_Location_SnowFraction_GAM <- df_Location_SnowFraction_GAM[order(df_Location_SnowFraction_GAM$doy),]
        
        #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
        df_Locations_SnowFraction_GAM <- rbind(df_Locations_SnowFraction_GAM, df_Location_SnowFraction_GAM)
        
        #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves
        
        #Refit GAM through data
        index <- which(df_Location_SnowFraction_GAM$outliers==FALSE)
        mod_gam <- with(df_Location_SnowFraction_GAM[index,], mgcv::gam(SnowFraction ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))
        
        #Use gam to make predictions on a more detailed interval
        df_Locations_SnowFraction_GAM_predictions_new <- data.frame(LocationID=i, doy=seq(min(df_Location_SnowFraction_GAM$doy), max(df_Location_SnowFraction_GAM$doy), 0.01))
        df_Locations_SnowFraction_GAM_predictions_new$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=df_Locations_SnowFraction_GAM_predictions_new, type="response")
        df_Locations_SnowFraction_GAM_predictions_new <- df_Locations_SnowFraction_GAM_predictions_new[order(df_Locations_SnowFraction_GAM_predictions_new$doy),]
        df_Locations_SnowFraction_GAM_predictions_new$year <- year_ID
        
        #Add predictions to df_Locations_SnowFraction_GAM_predictions dataframe:
        df_Locations_SnowFraction_GAM_predictions <- rbind(df_Locations_SnowFraction_GAM_predictions, df_Locations_SnowFraction_GAM_predictions_new)
  
    }
  
    #Change column LocationID to a factor:
    df_Locations_SnowFraction_GAM$LocationID <- as.factor(as.character(df_Locations_SnowFraction_GAM$LocationID))
    df_Locations_SnowFraction_GAM_predictions$LocationID <- as.factor(as.character(df_Locations_SnowFraction_GAM_predictions$LocationID))
    
    #Save dataframe with GAM fits for Location specific Snowfraction data
    #write.csv(df_Locations_SnowFraction_GAM, paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Locations_Data_SnowFraction_GAM.csv"), row.names = FALSE)
    write.csv(df_Locations_SnowFraction_GAM_predictions, paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Locations_GAM_Predictions_SnowFraction.csv"), row.names = FALSE)
    
    #(D) Plot the raw SnowFraction datapoints and gam predictions for each Location:
    
      #Plot SnowFraction and model predictions for all locations in a single plot
        p_Locations_SnowFraction = ggplot()+
         geom_point(data=df_Locations_SnowFraction_GAM, aes(x=doy, y=SnowFraction, fill=LocationID, col=LocationID))+
         geom_line(data=df_Locations_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict, col=LocationID)) +
         xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
         ylab(paste0("Snowcover fraction per location in ", year_ID)) +
         theme_tom()
  
    #(E) Calculate at which day of year the pixel-level SnowFraction is equal to Snowfraction_threshold_vector for each location using predictions from mod_gam
    
        #Setup parallel processing
        numCores <- detectCores()
        cl <- makePSOCKcluster(numCores)
        registerDoSNOW(cl)
        
        #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector
        results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                    pixelIDs_split = list(unique(df_Locations_SnowFraction_GAM_predictions$LocationID)), #levels LocationID (input needs to be a list)
                                                    df_pixel_y = df_Locations_SnowFraction_GAM_predictions, #dataframe containing GAM predictions
                                                    pixel_ID_column="LocationID", #Grouping column
                                                    y="SnowFraction_gam_predict", #response variable in GAM
                                                    x="doy", #predictor variable in GAM
                                                    pixel_gam_plots = FALSE, #Should GAM plots be created
                                                    y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated
            
        #Store dates of snowmelt per Location
        df_Snowmelt_Locations <- results[[1]]
        df_Snowmelt_Locations <- as.data.frame(do.call(rbind, df_Snowmelt_Locations))
        colnames(df_Snowmelt_Locations)[colnames(df_Snowmelt_Locations)=="pixel_ID"] <- "LocationID"
        colnames(df_Snowmelt_Locations)[colnames(df_Snowmelt_Locations)=="x_threshold"] <- "doy"
        colnames(df_Snowmelt_Locations)[colnames(df_Snowmelt_Locations)=="y_threshold"] <- "Snowfraction_threshold"
     
        #Turn parallel processing off and run sequentially again after this point
        stopCluster(cl)
        registerDoSEQ()
        
        #Add coordinates to each location
        if(Location_type=="points"){df_Snowmelt_Locations <- left_join(df_Snowmelt_Locations, df_locations, by=c("LocationID"))}
        
        #Add NDSI threshold
        df_Snowmelt_Locations$NDSI_threshold <- NDSI_threshold
        
        #Save dates of snowmelt per Location per SnowFraction threshold as a .csv file
        write.csv(df_Snowmelt_Locations, paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Locations_Snowmelt_Snowfraction.csv"), row.names = FALSE)
  
        #Create a separate plot with GAM predictions per location:
        
          #Create an empty list to store plots
          list_plots_snowfraction <- list(list())
          
          #Store all levels of LocationID and sort them in ascending order
          levelIDs <- sort(unique(df_Locations_SnowFraction_GAM_predictions$LocationID))
          
          #Loop through all levels of LocationID
          for(i in levelIDs){
            
            #i=levelIDs[1]
            
            #Create an index variable for parameter i
            i_index <- which(levelIDs == i)
            
            #select datasets for current Location:
            df_Location_SnowFraction_GAM <- df_Locations_SnowFraction_GAM[df_Locations_SnowFraction_GAM$LocationID==i,]
            df_Location_SnowFraction_GAM_predictions <- df_Locations_SnowFraction_GAM_predictions[df_Locations_SnowFraction_GAM_predictions$LocationID==i,]
            df_Location_SnowfractionDate <- df_Snowmelt_Locations[df_Snowmelt_Locations$LocationID==i,]
            
            #create plot for current Location and store it in list_plots_snowfraction:
            list_plots_snowfraction[[i_index]] <- ggplot()+ 
              geom_point(data=df_Location_SnowFraction_GAM[df_Location_SnowFraction_GAM$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
              geom_point(data=df_Location_SnowFraction_GAM[df_Location_SnowFraction_GAM$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
              geom_line(data=df_Location_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict), col="#1620de", lwd=1.25)+
              geom_point(data=df_Location_SnowfractionDate[!is.na(df_Location_SnowfractionDate$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)+
              geom_vline(xintercept = 150, colour="grey", lty=2)+
              geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
              xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
              ylab("Fraction of snow-covered pixels") +
              ggtitle(paste0("Polygon: ", i))+
              theme_tom()
            
          }
          
          #Plot SnowFraction and model predictions in a separate plot per Location
          plots_snowfraction <- list(list_plots_snowfraction)
          plots_snowfraction <- do.call(c, plots_snowfraction)
          plots_per_page = 25
          plots_snowfraction <- split(plots_snowfraction, ceiling(seq_along(plots_snowfraction)/plots_per_page))
          pdf(paste0(dir_Output, "/", timestamp, "_", data_ID, output_affix, "_Locations_Plot_Snowfraction.pdf"), width=20, height=16, onefile = TRUE)
          for (k in seq(length(plots_snowfraction))) { do.call("grid.arrange", plots_snowfraction[[k]]) }
          dev.off()

             

#####################################################################################################################################          
#####################################################################################################################################          
#####################################################################################################################################          
