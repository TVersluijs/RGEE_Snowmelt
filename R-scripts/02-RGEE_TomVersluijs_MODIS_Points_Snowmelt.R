#####################################################################################################################################

#The date of snow melt is calculated based on MODIS data for all locations specified in an input file. The user can specify a 
#bufferzone (radius) to depict the area in which snow melt will be analysed per location. No shapefile is required as input for 
#this script, which allows for input locations to be spaced far apart. The user can specify whether clouds and permanent water 
#bodies need to be masked. Snow melt is analysed within each location's buffer zone based on one of the following methods 
#(specified by the user by setting the parameter 'method'): 
# (I):   'avg_NDSI':     Calculate the average NDSI value over time within each point's buffer zone, fit a GAM through these data 
#                        and calculate when this model passes the specified NDSI threshold representing the moment of snow melt.
#                        In addition, time series of the average NDVI and NDMI are extracted within each point's buffer zone.
# (II):  'snowfraction': Calculate the fraction of pixels within each buffer zone over time where NDSI > 'NDSI_threshold', 
#                        fit a GAM through these data and extract the moment when this model passes a user-specified 
#                       'Snowfraction_threshold'.
# (III): 'pixel_gam':    Fit a GAM through the NDSI data for each pixel within each point's buffer zone, and calculate when this 
#                        function passes NDSI_threshold. Then use these pixel-specific dates of snowmelt to calculate a fraction
#                        of snowcovered pixels for each day of year. Then fit a GAM through these pixel-specific snowfraction 
#                        data and extract the moment when this model passes a user-specified 'Snowfraction_threshold'.

#The 'pixel_gam' method is preferred, because fitting a GAM through the NDSI data per pixel ensures that noise is filtered out
#on a pixel level. This clean data can then be used to look at the fraction of snow-covered pixels over time. This contrasts
#with the 'snowfraction' method because this simply calculates the fraction of pixels per timestep where NDSI>NDSI_threshold,
#which thus still includes all unfiltered noise. It is harder to justify the avg_NDSI method, because it is rather unclear
#what this average NDSI value entails.


#Copyright Tom Versluijs 2024-07-19. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R". 
#Note that a GoogleDrive is required. Important: make sure to run this script from within the "RGEE_Snowmelt.Rproj" project file.

#####################################################################################################################################

#I: Setup workspace

#####################################################################################################################################

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
             strip.background = element_rect(fill = "white", colour = "black", linewidth = rel(2)),
             complete = TRUE)}
       
       #(3): Load auxiliary functions
       source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "MODIS_AuxiliaryFunctions")
       sapply(source_files, source, chdir = TRUE) ; rm(source_files)
       
      #(4): Initialize earth engine and google drive
       rgee::ee_Initialize(user='tom.versluijs@gmail.com', drive=TRUE)

##################################################################################################################################
       
#II: Specify parameters of interest
       
##################################################################################################################################       
       
#Manually specify parameters of interest

  #(a): MODIS satellite

    #MODIS dataset
    MODIS_dataset <- "MODIS/061/MOD09GA"

    #Spatial resolution of satellite image
    #Note that this parameter only affects the calculation of the fraction of snowcover using the function
    #'AddSnowFraction' and the calculation of average bandvalues using the function 'Extract_BandValuesAtPoins'.
    #Detection of cloud-pixels occurs at resolution 'resolution_cldmsk' specified below. Detection of water-pixels
    #occurs at the native resolution of the MODIS dataset (500m).
    resolution=500 #default maximum resolution for MODIS = 500m

  #(b) Area of interest

    #Specify name of study area (used as prefix in output files)
    area_name="ZAC" #max length three characters

    #Coordinate reference system used for calculations
    #EPSG:4326 is recommended
    #EPSG:326XX might result in reduced computation time for areas located within a single UTM zone (i.e. planar coordinate system).
    crs <- "EPSG:4326"

  #(c) Point locations

    #Name of file with Locations of interest
    input_locations <- "TestLocations.csv"
    #Make sure it contains the coordinates in decimal degrees in the columns "LON_x" and "LAT_y"

    #Buffer radius around each point location (in meters)
    Buffer_radius_m=5000

  #(c) Dates

    #Specify the year of interest:
    year_ID <- "2022"

    #Date range of all images considered for analysis
    start_date <- paste0(year_ID, "-03-15") #choose date (well) before the first snowmelt occurs within the study site
    end_date <- paste0(year_ID, "-09-15") #choose date (well) after last date of tracking

  #(d) Snow detection

    #NDSI threshold for snow detection (specify multiple using c())
    NDSI_threshold_vector = c(0.4, 0.3, 0.5)

    #Define the snowcover fraction for which the date of its occurrence will be calculated (specify multiple using c())
    Snowfraction_threshold_vector = c(0.25, 0.5, 0.75)

    #Define the preferred method for the analysis of snowmelt
    method=c("avg_NDSI", "snowfraction", "pixel_gam") #either "avg_NDSI", "snowfraction", "pixel_gam", or a combination using c() (default is "pixel_gam").
    #(1) "avg_NDSI":     Calculates the average NDSI, NDVI and NDMI values within each point's bufferzone over time and extracts the moment 
    #                    when the NDSI value is equal to 'NDSI_threshold'.
    #(2) "snowfraction": Counts the fraction of pixels within the aoi with NDSI > 'NDSI_threshold' over time and extracts the moment 
    #                    when this fraction is equal to 'Snowfraction_threshold'.
    #(3) "pixel_gam":    Fits a GAM through the NDSI data per pixel and calculates when this function passes NDSI_threshold. Then use these
    #                    pixel-specific dates of snowmelt to calculate a fraction of snowcovered pixels for each day of year, and extract 
    #                    the moment when this fraction is equal to 'Snowfraction_threshold'.

    #should pixel-level image of date of snowmelt be exported to Google Drive (only applicable when method="pixel_gam")
    export_pixel_image_to_drive=TRUE 
    
  #(e): Cloud masking

    #Should clouds be masked from the analysis (default=TRUE).
    mask_clouds=TRUE
    
    #Define which MODIS cloud masking algorithm has been used ("PGE11", "MOD35", or "Combined)
    MODIS_cloud_masking_algorithm = 'MOD35' #default is MOD35
    
    #Maximum fraction of cloud cover allowed in each image
    max_cloud_fraction=0.75
    
    #Specify the resolution of the QABand   
    resolution_qaband=1000 #resolution of QA_band (MODIS=1000m)   

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
    #Larger values result in a more precise GAM-fit, but might result in overfitting.
    
    #Specify whether plots of the GAM NDSI fit per pixel should be generated as a pdf file (only applicable when method="pixel_gam")
    pixel_gam_plots=TRUE
    
  #(h): Counts of unmasked pixels
    
    #Should counts of the number of unmasked pixels per doy within the shapefile area be conducted (increases computation time)
    pixel_counts=TRUE
    

##################################################################################################################################
         
#III: Define some additional parameters (automated)
         
##################################################################################################################################
         
   #Create a unique data_ID
    if(nchar(area_name)>3){area_name <- substr(area_name, start = 1, stop = 3)}
    data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
    data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)
    
   #Create a timestamp variable
    timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")
    
   #First and last day of year in dataset  
    start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
    end_date_doy <- as.numeric(strftime(end_date, format = "%j"))  
        
   #Create a unique Asset folder for storing the generated datafiles (delete this folder if already present) 
    path_asset <- paste0(ee_get_assethome(), "/", data_ID)
    #tryCatch(ee_manage_assetlist(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
    tryCatch(ee_manage_delete(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
    ee_manage_create(path_asset=path_asset, asset_type="Folder")
    ee_manage_assetlist()

   #Create output folder
    if(dir.exists(paste0(here(), "/Output/MODIS/02_Points_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/MODIS/02_Points_Snowmelt"), recursive = TRUE)}
    if(dir.exists(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations"))==FALSE){dir.create(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations"), recursive = TRUE)}
    
   #Save all parameters and their values in the environment to a text file 
    file_conn <- file(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Parameters.txt"), "w")
    for (obj in setdiff(ls(), lsf.str())) {cat(paste(obj, "=", get(obj)), file = file_conn) ; cat("\n", file = file_conn)}
    close(file_conn)
    
    
##################################################################################################################################

#IV: Read dataframe with points of interest

##################################################################################################################################

  #Read a dataframe with points for which MODIS data needs to be extracted
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

#V: Iterate through all locations of interest and extract either (I) mean bandvalues, (II) the fraction of snowcover, or (III) pixel-level fraction of snowcover within each point's bufferzone over time 

################################################################################################################################################################################

 #Create an empty dataframe for storing the mean Bandvalue results per location (method="avg_NDSI")
  df_Locations_BandValues <- data.frame(LocationID=character(),
                                        Date=character(),
                                        doy=numeric(),
                                        NDSI=numeric(),
                                        NDVI=numeric(),
                                        NDMI=numeric(),
                                        LON_x=numeric(),
                                        LAT_y=numeric())

 #Create an empty dataframe for storing the Snowfraction results per location (method="SnowFraction")
  df_Locations_SnowFraction <- data.frame(LocationID=character(),
                                          Date=character(),
                                          doy=numeric(),
                                          SnowFraction=numeric(),
                                          NDSI_threshold=character(),
                                          LON_x=numeric(),
                                          LAT_y=numeric())
    
 #Create a dataframe for storing the Pixel-level Snowfraction results per location (method="pixel_gam")
  df_Locations_Pixel_SnowFraction <- data.frame(LocationID=character(),
                                                doy=numeric(),
                                                SnowFraction=numeric(),
                                                NDSI_threshold=character(),
                                                LON_x=numeric(),
                                                LAT_y=numeric())
  
 #Loop through all locations and extract the average band value and/or snowfraction within the location including buffer zone
  for(Location_i in 1:nrow(df_locations)){

 #(A): Print progress:
    cat("\n")
    print("#########################################################################################################")
    print(paste0("                  STARTING ANALYSIS FOR LOCATION ", Location_i, " OUT OF ", nrow(df_locations), " LOCATIONS"))
    print("#########################################################################################################")
    
 #(B): Load current point of interest and add required buffer zone:

    #Store current location in separate dataframe and store its locationID
     Locations <- df_locations[Location_i,]
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

    #Create aoi_Shapefile (corresponding to buffer zone of Location_i) as this is required input for several auxiliary functions.
     aoi_Shapefile <- ee$FeatureCollection(Locations)
    
 #(C): Extract MODIS data for the point of interest

    #Extract MODIS data for a given area and date range:
     MODIS_col <- ee$ImageCollection(MODIS_dataset)
     MODIS_col <- MODIS_col$
       filterBounds(aoi)$
       filterDate(start_date, end_date)

     #Note that MODIS_col is so far only filtered by 'aoi' and not yet clipped by aoi_Shapefile. Moreover, because
     #MODIS images are global composites, this single composite always intersects with aoi and this full image
     #is thus retained (i.e. filterBounds does not do anything).
     
    #Add a NDSI, NDVI, NDMI and NDWI band to the image collection
     MODIS_col <- MODIS_col$
       map(getNDSI)$
       map(getNDVI)$
       map(getNDMI)$
       map(getNDWI)

    # #Plot the first image in MODIS_col (for Debugging):
    #  image <- MODIS_col$filterDate(paste0(year_ID, "-06-10"), end_date)$first()
    # 
    # #Plot an RGB, NDSI and NDVI, NDMI and NDWI image for the first extracted satellite image
    #  Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
    #  Map$addLayer(image,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+ 
    #   Map$addLayer(image, list(bands=c("NDSI"), min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
    #   Map$addLayer(image, list(bands=c("NDVI"),min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
    #   Map$addLayer(image, list(bands=c("NDMI"),min=-1, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDMI')+
    #   Map$addLayer(image, list(bands=c("NDWI"),min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')

    #Create a timeseries gif of RGB images for the area of interest (For debugging)

      # #Check number of images in collection
      # MODIS_col$size()$getInfo()

      # #Create a timelapse video of RGB band
      # videoArgs <- list(dimensions=300, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
      # tryCatch({browseURL(MODIS_col$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})
    
 #(D): Filter and mask clouds within the image collection

    #There are many images that are almost fully covered in dense clouds. We thus first filter and mask clouds in the image collection.
    #We delete all images with a cloud cover fraction >= max_cloud_fraction. In the remaining images we mask all pixels that are
    #covered by clouds.
     if(mask_clouds==TRUE){

       #print message
       cat("\n")
       print("Cloud masking = TRUE")

       #The MODIS quality band 'state_1km' provides several cloud flag algorithms (See Wilson 2014, MODIS user guide). The cloud flags from
       #the 'MOD35'  algorithm are contained in bits 0-1 and have four categories: (0) clear, (1) cloudy, (2) mixed, (3) not set, assumed clear.
       #The cloud flags from the MODIS Land Team 'PGE11' internal cloud masking algorithm are contained in bit 10 and have two categories:
       #(0) no cloud, (1) cloud. Cirrus cloud flags are stored in bits 8-9 and have four categories: (0) none, (1) small, (2) average,
       #(3) high. The MOD35 cloud mask uses 22 MODIS bands, ecosystem/landcover type and other environmental data to identify clouds
       #(Wilson et al 2014), while the PGE11 internal cloud mask detects clouds based on two reflective testes and a thermal test. Usually
       #clouds are filtered if one or both of these two algorithms detects clouds (Wilson et al 2014).

       #Add the cloud fraction within the buffer zone surrounding Location to each separate image by mapping the cloud functions over the image collection
       MODIS_col <- MODIS_col$
         #Determine which pixels are either opague or cirrus clouds
         map(computeClouds)$
         #Add the fraction of cloud-covered pixels within the area of interest as image property
         map(AddCloudFraction)$
         #Add NULL to those images in which cloudfraction could not be calculated
         map(AddNULLCloudFraction)$
         #Add date characteristics to each image
         map(add_Date)

        # #Check if Cloud information has been added to the properties of each image (for debugging)
        # MODIS_col$first()$propertyNames()$getInfo()
        # MODIS_col$first()$bandNames()$getInfo()

        # # Visual check of cloud mask (for debugging) 
        #  image <- MODIS_col$filterDate(paste0(year_ID, "-08-15"), end_date)$first()$
        #    clipToCollection(aoi_Shapefile)$
        #    select("MOD35_clouds", "PGE11_clouds", "Combined_clouds", "MOD35_clouds_inv", "PGE11_clouds_inv", "Combined_clouds_inv",
        #           "sur_refl_b01", "sur_refl_b04", "sur_refl_b03")
        #  Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
        #  Map$addLayer(image,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
        #    Map$addLayer(image, list(bands='PGE11_clouds', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'PGE11_clouds')+
        #    Map$addLayer(image, list(bands='MOD35_clouds', min=0, max=1, opacity=0.5, palette=c('000000', 'orange')), 'MOD35_clouds')
        
        # #(for debugging) Compare the different MODIS cloud algorithms
        # 
        #   #Extract cloud fraction of all images in image collection for the area of interest (for debugging)
        #   Combined_clouds_Fraction <- unlist(MODIS_col$aggregate_array('Combined_clouds_Fraction')$getInfo())
        #   PGE11_clouds_Fraction <- unlist(MODIS_col$aggregate_array('PGE11_clouds_Fraction')$getInfo())
        #   MOD35_clouds_Fraction <- unlist(MODIS_col$aggregate_array('MOD35_clouds_Fraction')$getInfo())
        #   doy <- unlist(MODIS_col$aggregate_array('doy')$getInfo())
        # 
        #   #Replace -9999 values by NA (for debugging)
        #   Combined_clouds_Fraction[Combined_clouds_Fraction < -9000] <- NA
        #   PGE11_clouds_Fraction[PGE11_clouds_Fraction < -9000] <- NA
        #   MOD35_clouds_Fraction[MOD35_clouds_Fraction < -9000] <- NA
        #   doy[doy < -9000] <- NA
        # 
        #   #Combine all cloud fraction measures in a single dataframe (for debugging)
        #   clouds_MODIS <- data.frame(doy=doy,
        #                              Combined_clouds_Fraction=Combined_clouds_Fraction,
        #                              PGE11_clouds_Fraction=PGE11_clouds_Fraction,
        #                              MOD35_clouds_Fraction=MOD35_clouds_Fraction)
        # 
        #   #Assess relationship between both cloud algorithms (for debugging)
        #   ggplot()+
        #     geom_line(aes(x=doy, y=MOD35_clouds_Fraction), col="red", linewidth=1)+
        #     geom_line(aes(x=doy, y=PGE11_clouds_Fraction), col="black", linewidth=1)+
        #     #geom_line(aes(x=doy, y=Combined_clouds_Fraction), col="blue", linewidth=1)+
        #     theme_tom()
        # 
        #   #Correlation between MOD35 and PGE11 (for debugging)
        #   ggplot()+
        #     geom_point(aes(x=MOD35_clouds_Fraction, y=PGE11_clouds_Fraction), col="black")+
        #     geom_abline(slope=1, intercept=0, lty=2)+
        #     theme_tom()
        # 
        #   #Based on these plots we chose MOD35 as it seems to be a bit more conservative with detecting pixels as
        #   #being cloud-covered (i.e. more points above the line than below)

       #Automatically detect which cloud algorithm was used and what the name is of the clouds_fraction parameter
        cloud_algorithm <- paste0(MODIS_cloud_masking_algorithm, "_clouds")
        clouds_fraction <- paste0(MODIS_cloud_masking_algorithm, "_clouds_Fraction")
       
       #Exclude all images from the image collection that have a CloudFraction value > max_cloud_fraction, and mask all cloud pixels within the remaining images
        MODIS_clouds_filtered <- MODIS_col$
          #Only remove images that are fully cloud-covered (PGE11_clouds_Fraction=1.0)
          filter(ee$Filter$lt(clouds_fraction, max_cloud_fraction))$
          #Apply cloudmask for individual pixels
          map(AddCloudMask)

       # #Create timelapse video of the cloud filtered/masked RGB images (for debugging)
       #  videoArgs <- list(dimensions=380, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=10000, gamma=c(1.9, 1.7, 1.7))
       #  tryCatch({browseURL(MODIS_clouds_filtered$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$getVideoThumbURL(videoArgs))}, error = function(cond){return("Too many pixels. Reduce dimensions.")})

       }
     if(mask_clouds==FALSE){

       #print message
        print("Cloud masking = FALSE")

       #Add Date and Time characteristics to each image
       MODIS_clouds_filtered <- MODIS_col$
         #Add date characteristics to each image
         map(add_Date)

     }

    #Store default MODIS image projection
     modisProjection <- MODIS_clouds_filtered$first()$select("NDSI")$projection()
     
    #Note that MODIS_clouds_filtered is not yet clipped by aoi_Shapefile (only by aoi)!

 #(E): Mask permanent waterbodies (ponds, lakes, rivers, sea) within the image collection   
    
    #Mask permanent waterbodies if mask_water==TRUE
     if(mask_water==TRUE){

       #(A): print message
         print("Water masking = TRUE")

       #(B): Extract permanent waterbodies from the Terra Land Water Mask map (250m resolution)

         #Load auxiliary function
         compute_Water_MODIS=compute_Water_MODIS #sourced

         #Extract pixels corresponding to permanent waterbodies
         water_mask <- compute_Water_MODIS()

       #(C): Apply the water masking function to each image in the collection:

         # #Plot permanent waterbodies(for debugging)
         # Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
         # Map$addLayer(MODIS_clouds_filtered$filterDate(paste0(year_ID, "-08-04"), end_date)$first(), list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
         # Map$addLayer(water_mask,list(min=0, max = 1, palette = c('ffffff', '000000')), 'Water_mask')

         #Load auxiliary function to mask water-pixels
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
    
 #(F): Create a composite image for each day by mosaicking all images from the that day    
    
     #For Sentinel-2 data there might be days for which multiple satellite photos are available that can slightly overlap. In
     #that case we need to deal with these overlapping pixels and we can do that by making a composite image by selecting 
     #the overlaping pixel with e.g. the least cloudcover. In contrast, MODIS data by default already is composite image for 
     #each day, thus there is no need to conduct this step.

 #(G): Count the total number of unmasked pixels and the total number of pixels per doy within aoi_Shapefile
     
     if(pixel_counts==TRUE){
     
     #Add pixel counts within the area of interest to each separate image by mapping the pixel count functions over the image collection
      MODIS_clouds_filtered <- MODIS_clouds_filtered$
        #Count number of unmasked pixels within aoi_Shapefile
        map(AddPixelCount)$
        #Add NULL to those images in which pixel count could not be calculated
        map(AddNULLPixelCount)
     
     #Extract pixel_counts of all images in image collection for aoi_Shapefile
     
       #create a current timestamp to prevent identical names on Google Drive
        current_timestamp0 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
       
       #We use ee_table_to_drive() to prevent memory limits
        task_vector0 <- ee_table_to_drive(
          collection = MODIS_clouds_filtered,
          description = paste0(current_timestamp0, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Data_Pixel_Counts_polygon"),
          fileFormat = "CSV",
          selectors = c('doy', 'unmasked', 'total')
          )
       
       #Monitor the task
        task_vector0$start()
        print("Count the number of unmasked pixels within the shapefile per doy:")
        ee_monitoring(task_vector0, quiet=T, max_attempts=1000000)
       
       #Export results to local folder
        exported_stats <- ee_drive_to_local(task = task_vector0, dsn=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Data_Pixel_Counts_polygon"))
        df_pixelcount <- read.csv(exported_stats)
        unlink(exported_stats)
     
     #Replace -9999 values by NA (for debugging)
      df_pixelcount$unmasked[df_pixelcount$unmasked < -9000] <- NA
      df_pixelcount$total[df_pixelcount$total < -9000] <- NA
      df_pixelcount$doy[df_pixelcount$doy < -9000] <- NA
     
     #Add missing dates with 0 unmasked pixels to the dataframe
      df_doy_missing <- data.frame(doy=seq(start_date_doy, end_date_doy)[!(seq(start_date_doy, end_date_doy) %in% df_pixelcount$doy)],
                                   unmasked=0,
                                   total=max(df_pixelcount$total))
      df_pixelcount <- rbind(df_pixelcount, df_doy_missing)
      df_pixelcount <- df_pixelcount[order(df_pixelcount$doy),]
      write.csv(df_pixelcount, file=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Data_Pixel_Counts_polygon.csv"), quote=FALSE, row.names=FALSE)
     
     #Create barplot with the pixel counts per day of year within aoi_Shapefile
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
     
     #Save barplot
      pdf(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Plot_Pixel_Counts_polygon.pdf"), width=12, height=8)
      print(p_pixelcounts)
      dev.off()
      
     }
     
 #(H): Extract MODIS mean band values (NDSI, NDVI, NDMI) within the buffer zone of 'Location' for all images in the image collection

     if("avg_NDSI" %in% method){

       #Print message
       cat("\n")
       print("--------------------------------------------------------------------------------------------------------------------------")
       print(paste0("METHOD: 'avg_NDSI' - CALCULATING THE AVERAGE NDSI, NDVI AND NDMI FOR LOCATION ", Location_i))
       print("--------------------------------------------------------------------------------------------------------------------------")
       
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
       #has been defined. This function calculates mean NDSI, NDVI and NDMI within aoi_Shapefile (buffer zone of Location_i)
        Extract_BandValuesAtPoins = Extract_BandValuesAtPoins #sourced

       #Iterate over the ImageCollection and select appropriate band values
        FC_merged <- ee$FeatureCollection(MODIS_clouds_filtered$select("NDSI", "NDVI", "NDMI")$iterate(Extract_BandValuesAtPoins, FC_initial))

       #Transform the feature collection with average Band values per day per location to a dataframe by exporting it

         #(I): Google Drive (preferred option, slightly faster than using Asset folder)

           #create a current timestamp to prevent identical names on Google Drive
            current_timestamp1 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
        
           #Setup task
            task_vector1 <- ee_table_to_drive(
             collection= FC_merged,
             description = paste0(current_timestamp1, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Data_MeanBandvalues"),
             folder="RGEE_tmp",
             fileFormat="CSV",
             selectors=c('LocationID', 'Date', 'NDSI', 'NDVI', 'NDMI')
             )

           #Run and monitor task
            task_vector1$start()
            ee_monitoring(task_vector1, quiet=TRUE, max_attempts=1000000)

           #Import results
            exported_stats <- ee_drive_to_local(task=task_vector1, dsn=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Data_MeanBandvalues"))
            df_Locations_Bandvalues_new <- read.csv(exported_stats)
            unlink(exported_stats)
            
           #Add day of year
            df_Locations_Bandvalues_new$doy <- as.numeric(format(as.POSIXct(df_Locations_Bandvalues_new$Date, format = "%Y-%m-%d %H:%M:%S"), "%j"))
            
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
         #    ee_monitoring(task_vector, quiet=TRUE, max_attempts=1000000)
         #
         #   #Import results
         #    ee_fc <- ee$FeatureCollection(assetid)
         #    df_Locations_Bandvalues_new <- data.frame(LocationID=unlist(ee_fc$aggregate_array('LocationID')$getInfo()),
         #                                           Date=unlist(ee_fc$aggregate_array('Date')$getInfo()),
         #                                           doy=unlist(ee_fc$aggregate_array('doy')$getInfo()),
         #                                           NDSI=round(as.numeric(unlist(ee_fc$aggregate_array('NDSI')$getInfo())), 5),
         #                                           NDVI=round(as.numeric(unlist(ee_fc$aggregate_array('NDVI')$getInfo())), 5),
         #                                           NDMI=round(as.numeric(unlist(ee_fc$aggregate_array('NDMI')$getInfo())), 5)
         #                                           )

           #Change -9999 to NA
            df_Locations_Bandvalues_new$NDSI[df_Locations_Bandvalues_new$NDSI < -9000] <- NA
            df_Locations_Bandvalues_new$NDVI[df_Locations_Bandvalues_new$NDVI < -9000] <- NA
            df_Locations_Bandvalues_new$NDMI[df_Locations_Bandvalues_new$NDMI < -9000] <- NA

           #Sort dataframe by LocationID and doy
            index <- with(df_Locations_Bandvalues_new, order(LocationID, doy))
            df_Locations_Bandvalues_new <- df_Locations_Bandvalues_new[index,]

           #Add coordinates to each location
            df_Locations_Bandvalues_new <- left_join(df_Locations_Bandvalues_new, df_locations[,c("LON_x", "LAT_y", "LocationID")], by=c("LocationID"))

           #Add dataframe for current location to dataframe from previous iterations:
            df_Locations_BandValues <- rbind(df_Locations_BandValues, df_Locations_Bandvalues_new)

           # #Save dataframe for current location
           #  write.csv(df_Locations_Bandvalues_new, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Data_MeanBandvalues_polygon.csv"), row.names = FALSE)

       }
     
 #(I): Extract the fraction of snow covered pixels within the buffer zone of 'Location' for all images in the image collection
     
     if("snowfraction" %in% method){      
    
       #Print message
       cat("\n")
       cat("\n")
       print("--------------------------------------------------------------------------------------------------------------------------")
       print(paste0("METHOD: 'snowfraction' - CALCULATING THE FRACTION OF PIXELS WITH NDSI > NDSI_threshold FOR LOCATION ", Location_i))
       print("--------------------------------------------------------------------------------------------------------------------------")
       
       #Run the analysis for each level of NDSI_threshold_vector
       for(NDSI_threshold in NDSI_threshold_vector){  
         
           #Print message
           print(paste0("  -Start analysis for NDSI_threshold = ", NDSI_threshold))
         
           #Store NDSI_threshold as a character (used for naming of outputs)
           NDSI_threshold_char <- gsub("\\.", "_", as.character(NDSI_threshold))
           
         #(I.1): Add a binary snow cover band to each image and calculate the fraction of snow covered pixels within the buffer zone (aoi_Shapefile)
  
           #Map Snow computation functions over the mosaicked image collection
           MODIS_clouds_filtered_snow <- MODIS_clouds_filtered$
             #Determine which pixels are snow-covered (NDSI > NDSI threshold)
             map(computeSnow)$
             #add the fraction of snow covered pixels within aoi_Shapefile as an image property (excluding cloud masked pixels from the calculations)
             map(AddSnowFraction)$
             #Add a NULL value to images for which the snow fraction could not be calculated
             map(AddNULLSNOW)
  
             # #Mask and display the binary layer (for debugging).
             # img_snow <- MODIS_clouds_filtered_snow$filterDate(paste0(year_ID, "-05-25"), end_date)$first()
             # Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
             # Map$addLayer(img_snow,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
             # Map$addLayer(img_snow$select('SNOW'), list(min=0, max = 1, palette = c('ffffff', 'orange')))+
             # Map$addLayer(img_snow$select('NDSI'), list(min=-1, max = 1, palette = c('lightblue', 'darkblue')))
         
         #(I.2): Extract fraction of snowcover within the buffer zone for each image in the image collection
  
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
           FC_merged <- ee$FeatureCollection(MODIS_clouds_filtered_snow$iterate(extract_snowcover, FC_initial))
  
           #Export the feature collection as a .csv table
           #We export the data instead of using aggregate_array() as the latter might fail due to computation timeouts.
  
             #create a current timestamp to prevent identical names on Google Drive
              current_timestamp2 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
           
             #Setup task
              task_vector2 <- ee_table_to_drive(
                collection= FC_merged,
                description = paste0(current_timestamp2, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_NDSI", NDSI_threshold_char, "_Data_SnowFraction"),
                folder="RGEE_tmp",
                fileFormat="CSV",
                selectors=c('SnowFraction', 'Date')
                )
  
             #Run and monitor task
              task_vector2$start()
              ee_monitoring(task_vector2, quiet=TRUE, max_attempts=1000000)
  
             #Import results
              exported_stats <- ee_drive_to_local(task=task_vector2, dsn=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_NDSI", NDSI_threshold_char, "_Data_SnowFraction"))
              df_Locations_SnowFraction_new <- read.csv(exported_stats)
              unlink(exported_stats)
              
             #Add day of year
              df_Locations_SnowFraction_new$Date <- as.POSIXct(df_Locations_SnowFraction_new$Date, format="%Y-%m-%d %H:%M:%S")
              df_Locations_SnowFraction_new$doy <- as.numeric(strftime(df_Locations_SnowFraction_new$Date, format = "%j"))
  
             #Remove NAs in the SnowFraction variable
              df_Locations_SnowFraction_new$SnowFraction[df_Locations_SnowFraction_new$SnowFraction < -9000] <- NA #replace -9999 by NA
              df_Locations_SnowFraction_new <- df_Locations_SnowFraction_new[!is.na(df_Locations_SnowFraction_new$SnowFraction), ]
  
             #Add LocationID
              df_Locations_SnowFraction_new$LocationID <- Location_ID
   
             #Sort dataframe by LocationID and doy
              index <- with(df_Locations_SnowFraction_new, order(LocationID, doy))
              df_Locations_SnowFraction_new <- df_Locations_SnowFraction_new[index,]
              df_Locations_SnowFraction_new <- df_Locations_SnowFraction_new[,c("LocationID", "Date", "doy", "SnowFraction")]
  
             #Add NDSI_threshold as a new column
              df_Locations_SnowFraction_new$NDSI_threshold <- as.factor(NDSI_threshold)
              
             #Add coordinates to each location
              df_Locations_SnowFraction_new <- left_join(df_Locations_SnowFraction_new, df_locations[,c("LON_x", "LAT_y", "LocationID")], by=c("LocationID"))
   
             #Add dataframe for current location to dataframe from previous iterations:
              df_Locations_SnowFraction <- rbind(df_Locations_SnowFraction, df_Locations_SnowFraction_new)
  
             # #For debugging
             #  ggplot() + geom_point(data=df_Locations_SnowFraction, aes(x=doy, y=SnowFraction, col=NDSI_threshold)) + theme_classic()
               
             # #Save dataframe for current location
             #  write.csv(df_Locations_SnowFraction_new, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_", Location_ID, "_Data_SnowFraction_polygon.csv"), row.names = FALSE)
           }
       
       }
           
 #(J): Calculate the date of snowmelt for every pixel within aoi_Shapefile by fitting a GAM through the pixel-specific NDSI data
     
     if("pixel_gam" %in% method){
       
        #(J.1): Transform each MODIS image to a feature Collection of NDSI values for all pixels within aoi (i.e. bounding box!)

          #Print message
          cat("\n")
          print("--------------------------------------------------------------------------------------------------------------------------")
          print(paste0("METHOD: 'pixel_gam' - CALCULATING THE FRACTION OF SNOWCOVER ON A PIXEL LEVEL FOR LOCATION ", Location_i))
          print("--------------------------------------------------------------------------------------------------------------------------")
       
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
           Extract_BandValuesAtPixels = Extract_BandValuesAtPixels 
           
           #Note that within 'Extract_BandValuesAtPixels' the region is set to 'aoi'. We do this because we can then calculate the
           #date of snowmelt for all pixels including those that fall on the border of the buffer zone 'aoi_Shapefile'. Otherwise 
           #these 'border-pixels' would have been excluded which would create an incomplete plot at the edge of the image. Instead, 
           #we thus first calculate snowmelt for all pixels in 'aoi', then crop this to 'aoi_Shapefile' which retains these border-
           #pixels, then create a plot of this image and then store the date of snowmelt for the subset of pixels within this image
           #(i.e. within aoi_Shapefile). IMPORTANTLY, if this plot is not required, we can simply change the function 'Extract_Band
           #ValuesAtPixels' above to use region=aoi_Shapefile instead of region=aoi to directly calculate the date of snowmelt only
           #for those pixels within aoi_Shapefile.

          #Iterate over the ImageCollection (output is a large feature collection)
           FC_merged <- ee$FeatureCollection(MODIS_clouds_filtered$select("NDSI")$iterate(Extract_BandValuesAtPixels, FC_initial))
           #FC_merged$first()$getInfo() #for debugging

           #Note that at this point all pixels within 'aoi' are included (not only those in 'aoi_Shapefile'!)
           
        #(J.2): Transform feature collection to a dataframe:

          #create a current timestamp to prevent identical names on Google Drive
           current_timestamp3 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
           
          #We use ee_table_to_drive() to prevent memory limits
           a=Sys.time()
           task_vector3 <- ee_table_to_drive(
             collection = FC_merged,
             description = paste0(current_timestamp3, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_NDSI_bbox"),
             fileFormat = "CSV",
             selectors = c('NDSI', 'Date', 'lat', 'lon')
             )

           task_vector3$start()
           cat("\n")
           print("  -Step 1: Transform each MODIS image to a feature Collection of NDSI values for all pixels:")
           ee_monitoring(task_vector3, quiet=TRUE, max_attempts = 100000)

           exported_stats <- ee_drive_to_local(task = task_vector3, dsn=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_NDSI_bbox"))
           df_pixel_ndsi <- read.csv(exported_stats)
           b=Sys.time()
           #print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))

          # #Load dataframe (takes ca 2 minutes):
          #  df_pixel_ndsi <- read.csv(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_NDSI_bbox.csv"))

          #Add day of year
           df_pixel_ndsi$doy <- as.numeric(format(as.POSIXct(df_pixel_ndsi$Date, format = "%Y-%m-%d %H:%M:%S"), "%j"))
          
          #Make sure each latitude/longitude combination gets its own pixel_ID (takes ca 1 minute):
           df_pixel_ndsi$pixel_ID <- paste0(format(round(df_pixel_ndsi$lat, 5), nsmall = 5), "_", format(round(df_pixel_ndsi$lon, 5), nsmall = 5))
       
        #(J.3): Calculate the date of snowmelt for each pixel in the dataframe (i.e. within area 'aoi')
  
           #Loop through all pixel_ID's, select dataframe for that pixel containing NDSI values
           #as measured in all images in the image collection, fit gam through the NDSI ~ doy data,
           #determine date when NDSI<NDSI_threshold and store this date of snowmelt together with
           #the pixel_ID in a new dataframe.
  
           #Specify NDSI threshold (values larger than this threshold are considered snow covered):
            NDSI_threshold_vector=NDSI_threshold_vector
  
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
           #year the predicted NDSI value of this GAM changes from above outlier_threshold to below. The code employs parallel 
           #processing using foreach and %dopar% on four local computer cores.
            f_detect_threshold_date_parallel <- f_detect_threshold_date_parallel #sourced
  
           #Run the 'f_detect_threshold_date_parallel' function over all data subsets, combine the results and save the resulting dataframe and plots
            cat("\n")
            print("  -Step 2: Calculate the date of snowmelt for each pixel within each location's bounding box:")
            results <- lapply(1:length(pixelIDs_split), FUN=f_detect_threshold_date_parallel,
                              pixelIDs_split=pixelIDs_split, df_pixel_y=df_pixel_ndsi, pixel_ID_column="pixel_ID",
                              y="NDSI", x="doy", pixel_gam_plots=pixel_gam_plots, y_threshold=NDSI_threshold_vector)
  
            #Clean up the cluster after finishing the parallel runs
            stopCluster(cl)
  
            #Turn parallel processing off and run sequentially again after this point
            registerDoSEQ()
  
           #Store date of snowmelt per pixel within aoi (bounding box!) as a dataframe
            df_pixel_snowmelt <- lapply(results, "[[", 1)
            df_pixel_snowmelt <- as.data.frame(do.call(rbind, do.call(c, df_pixel_snowmelt)))
            colnames(df_pixel_snowmelt)[colnames(df_pixel_snowmelt)=="x_threshold"] <- "doy_snowmelt"
            colnames(df_pixel_snowmelt)[colnames(df_pixel_snowmelt)=="y_threshold"] <- "NDSI_threshold"
            #write.csv(df_pixel_snowmelt, file=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_Snowmelt_bbox.csv"), quote = FALSE, row.names=FALSE)
            #df_pixel_snowmelt <- read.csv(file=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_Snowmelt_bbox.csv"), header=TRUE)

           #Store GAM plots per pixel within aoi (bounding box!) as a list
            if(pixel_gam_plots==TRUE){
              plot_pixel_snowmelt <- lapply(results, "[[", 2)
              plot_pixel_snowmelt <- do.call(c, plot_pixel_snowmelt)
              }
  
           #The dataframe df_pixel_snowmelt now contains the date of snowmelt for each individual pixel_ID within the
           #aoi (i.e. the defined bounding box). To be able to plot these data we transform this dataframe to a feature
           #collection and then transform this feature collection to an image. We can then clip this image by aoi_Shapefile
           #to get only those pixels within the buffer zone of each point (and not only in the bounding box).
           
        #(J.4): Transform df_pixel_snowmelt to a feature collection with random geometry

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
                cat("\n")
                print("  -Step 3: Transform dataframe with date of snowmelt per pixel to a feature collection with random geometry:")
                df_sf_tmp <- st_as_sf(x = df_tmp,
                                      coords = c("lon", "lat"),
                                      crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

               #Change sf object to an earth engine feature collection by uploading it to the asset folder
                FC_tmp <- sf_as_ee(
                  x = df_sf_tmp,
                  assetId = paste0(path_asset, "/", timestamp, "_", data_ID, "_Location_", Location_i, "_df_pixel_snowmelt_", i),
                  overwrite = FALSE,
                  monitoring = TRUE,
                  quiet = TRUE,
                  via = 'getInfo_to_asset')
                #FC_tmp <- sf_as_ee(df_sf_tmp)

               #Add feature collection to an expanding feature collection:
                FC_initial <- FC_initial$merge(FC_tmp)
                FC_initial <- ee$FeatureCollection(FC_initial)

               #print progress
                #print(paste0("Progress: ", (100*i)/(length(rowIDs_split))))
                #print(paste0("Size of FC_initial: ", FC_initial$size()$getInfo()))

             }
             FC_pixels_snowmelt <- ee$FeatureCollection(FC_initial)
             #FC_pixels_snowmelt$first()$getInfo()
             #FC_pixels_snowmelt$size()$getInfo()

             #Inspect assets folder:
             #  ee_manage_quota()
             #  ee_manage_assetlist(path_asset)

        #(J.5): Add geometry (latitude and longitude) of each pixel_ID to FC_pixels_snowmelt

            #The feature collection that we want to construct should contain 'pixel_ID' and 'doy_snowmelt' as properties and should contain
            #the original geometry corresponding to each pixel_ID. So far, FC_pixels_snowmelt contains a separate feature for each pixel, where
            #each feature contains doy_snowmelt and pixel_ID as a property. However, the geometry (lan/lon) of each feature is randomly chosen.
            #We need to make sure that the actual geometry matching each pixel_ID is added instead of this random geometry.

            #To obtain the corresponding geometry (lat/lon) of every pixel_ID within aoi, we sample from a single image on a 'resolution'
            #resolution using img$sample(). This gives a single distinct feature for each pixel, including their geometry. We can then re-construct
            #the property pixel_ID for every feature in this feature collection. The resulting feature collection is called FC_pixels_distinct.
            #The final step is then to join FC_pixels_snowmelt to FC_pixels_distinct based on an inner join with pixel_ID.

             #Print progress:  
              cat("\n")
              print("  -Step 4: Add geometry (latitude and longitude) to the feature collection with dates of snowmelt per pixel:")
              
             #(A): Select a single image from the image collection
              img <- MODIS_clouds_filtered$first()$select("NDSI")

             #(B): Create a feature collection of points within aoi (where each point corresponds to the center of a pixel in the image)
              FC_pixels_distinct <- img$sample(
               region=aoi, #Sample all pixels within aoi (i.e. bounding box)
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

              #FC_pixels_distinct$first()$getInfo()
              #FC_pixels_distinct$size()$getInfo()

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
                current_timestamp4 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
                
               #Delete FC_pixels_snowmelt_optimized if it already occured in the asset folder:
                tryCatch({ee_manage_delete(paste0(path_asset, "/", current_timestamp4, "_", data_ID, "_Location_", Location_i, "_FC_pixels_snowmelt_optimized"))},
                         error = function(cond){return("Path did not yet exist - no folder deleted")})

               #Upload to asset folder:
                assetid2 <- paste0(path_asset, "/", current_timestamp4, "_", data_ID, "_Location_", Location_i, "_FC_pixels_snowmelt_optimized")
                task_vector4 <- ee_table_to_asset(
                    collection = FC_Combined,
                    overwrite = FALSE,
                    assetId = assetid2
                    )
                task_vector4$start()
                cat("\n")
                print("  -Step 5: Optimize further calculations with this feature collection by uploading it to the asset folder:")
                ee_monitoring(task_vector4, quiet=TRUE, max_attempts = 100000)

               #Check assets folder:
                #ee_manage_quota()
                #ee_manage_assetlist(path_asset)

               #Get feature collection from asset folder and create FC_pixels_snowmelt_optimized
                #assetid2=paste0(path_asset, "/", current_timestamp4, "_", data_ID, "_Location_", Location_i, "_FC_pixels_snowmelt_optimized")
                FC_pixels_snowmelt_optimized <- ee$FeatureCollection(assetid2)
                #FC_pixels_snowmelt_optimized$first()$getInfo()
                #FC_pixels_snowmelt_optimized$size()$getInfo()

         #Conduct the following steps separately for all levels of NDSI_threshold_vector    
         for(NDSI_threshold in NDSI_threshold_vector){  
                  
           #Print message
           cat("\n")
           print(paste0("START ANALYSIS FOR NDSI-THRESHOLD: ", NDSI_threshold))
           
           #Store NDSI_threshold as a character (used for naming of outputs)
           NDSI_threshold_char <- gsub("\\.", "_", as.character(NDSI_threshold))
           
           #(J.6): Transform the Feature collection FC_pixels_snowmelt_optimized to an image (with doy_snowmelt as an image band)

               #(A): Reduce feature collection FC_pixels_snowmelt_optimized to an Image with a 500m resolution:
                image_snowmelt <- FC_pixels_snowmelt_optimized$
                  filter(ee$Filter$notNull(list('doy_snowmelt')))$ #pre-filter data for nulls that cannot be turned into an image
                  filter(ee$Filter$neq(name='doy_snowmelt', value=-9999))$ #pre-filter data for -9999 values
                  filter(ee$Filter$eq(name='NDSI_threshold', value=NDSI_threshold))$ #filter for a specific NDSI threshold
                  reduceToImage(properties=list('doy_snowmelt'), reducer=ee$Reducer$first()$setOutputs(list('doy_snowmelt')))$
                  #reproject(crs=crs, crsTransform=NULL, scale=resolution) #ensures sure that reduceToImage above is done in the crs projection
                  reproject(crs=modisProjection, crsTransform=NULL, scale=resolution) #ensures sure that reduceToImage above is done in the MODIS projection

                #image_snowmelt$projection()$getInfo()

               #(B): Clip image along The specified shapefile.
                image_snowmelt <- image_snowmelt$clipToCollection(aoi_Shapefile)

               # #(C): Extract day of snowmelt in year of interest for a single point (for debugging)
               #  ee_extract(x=image_snowmelt, y=ee$Geometry$Point(coordinates_point[1], coordinates_point[2]), fun=ee$Reducer$first(), scale=resolution, sf=TRUE)

               # #(D): Plot snowmelt day of year as a coloured image (for debugging)
               #  MODIS_image <- MODIS_clouds_filtered$filterDate(paste0(year_ID, "-07-20"), end_date)$first()#$reproject(crs=crs, crsTransform=NULL, scale=resolution)
               #  Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
               #  Map$addLayer(MODIS_image,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
               #  Map$addLayer(image_snowmelt,list(bands="doy_snowmelt", min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red')), 'Snowmelt_doy')+
               #  Map$addLayer(water_mask$clipToCollection(aoi_Shapefile)$updateMask(water_mask$neq(0)),list(min=0, max = 1, palette = c('ffffff', 'darkblue')), 'Water_mask')
                 
               #(E): Export snowmelt image to Google Drive (takes c.a. 2 minutes):
                if(export_pixel_image_to_drive==TRUE){
                
                   #create a current timestamp to prevent identical names on Google Drive
                    current_timestamp5 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
                  
                   #Create task to export the original doy_snowmelt image to Google Drive
                    task_vector5 <- ee_image_to_drive(
                     fileFormat='GeoTIFF',
                     image=image_snowmelt,
                     description=paste0(current_timestamp5, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt'),
                     region=aoi,
                     #scale=ee$Number(resolution), #defaults to native resolution of image asset.
                     crs="EPSG:3857", #Coordinate reference system of projection of exported image
                     maxPixels=1e9, #maximum allowed number of pixels in exported image
                     dimensions=ee$Number(1024) #maximum dimension
                     )

                   #Start and monitor export task:
                    task_vector5$start()
                    cat("\n")
                    print("  -Step 6: Transform the featurecollection to a snowmelt image and export it to Google Drive:")
                    ee_monitoring(task_vector5, quiet=TRUE, max_attempts = 1000000)
                    ee_drive_to_local(task = task_vector5, dsn=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt'))
                    
                 #(F): Export RGB image to Google Drive (takes c.a. 2 minutes):

                   #Convert original image to an RGB image:
                    image_snowmelt_RGB <- image_snowmelt$visualize(bands=c('doy_snowmelt'), min=start_date_doy, max=end_date_doy, palette=c('green', 'yellow', 'red'))
                    #ee_print(image_snowmelt_RGB)
                    #image_snowmelt_RGB$projection()$getInfo()

                   #create a current timestamp to prevent identical names on Google Drive
                    current_timestamp6 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
                    
                   #Create task to export RGB image to Google Drive:
                    task_vector6 <- ee_image_to_drive(
                      fileFormat='GeoTIFF',
                      image=image_snowmelt_RGB,
                      description=paste0(current_timestamp6, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt_RGB'),
                      region=aoi,
                      #scale=ee$Number(resolution), #defaults to native resolution of image asset.
                      crs="EPSG:3857", #Coordinate reference system of projection of exported image
                      maxPixels=1e9, #maximum allowed number of pixels in exported image
                      dimensions=ee$Number(1024) #maximum dimension
                      )

                   #Start and monitor export task:
                    cat("\n")
                    print("  -Step 7: Transform snowmelt image to an RGB image and export it to Google Drive:")
                    task_vector6$start()
                    ee_monitoring(task_vector6, quiet=TRUE, max_attempts = 1000000)
                    ee_drive_to_local(task = task_vector6, dsn=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, '_Pixel_Image_DoySnowmelt_RGB'))
                    
                }
                      
           #(J.7): Extract the date of snowmelt for all pixels within image_snowmelt (i.e. clipped by aoi_Shapefile)

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
                  current_timestamp7 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
                
                 #We use ee_table_to_drive() to prevent memory limits
                  a=Sys.time()
                  task_vector7 <- ee_table_to_drive(
                    collection = FC_snowmelt,
                    description = paste0(current_timestamp7, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Snowmelt_shapefile"),
                    fileFormat = "CSV",
                    selectors = c('doy_snowmelt', 'lat', 'lon')
                    )

                 #Execute task
                  task_vector7$start()
                  cat("\n")
                  print("  -Step 8: Extract the date of snowmelt from the snowmelt image for all pixels within aoi_Shapefile:")
                  ee_monitoring(task_vector7, quiet=TRUE, max_attempts = 1000000)
                  exported_stats <- ee_drive_to_local(task = task_vector7, dsn=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Snowmelt_shapefile"))
                  df_pixel_snowmelt_shapefile <- read.csv(exported_stats)
                  unlink(exported_stats)
                  b=Sys.time()
                  #print(paste0("Computation finished in ",  round(as.numeric(difftime(b, a, units="mins")),2), " minutes"))

                 #Make sure each latitude/longitude combination gets its own pixel_ID
                  df_pixel_snowmelt_shapefile$pixel_ID <- paste0(format(round(df_pixel_snowmelt_shapefile$lat, 5), nsmall = 5), "_", format(round(df_pixel_snowmelt_shapefile$lon, 5), nsmall = 5))

                 #Add NDSI_threshold as a new column
                  df_pixel_snowmelt_shapefile$NDSI_threshold <- as.factor(NDSI_threshold)
                  
                 #Only select columns "pixel_ID",  "doy_snowmelt" and "NDSI_threshold"
                  df_pixel_snowmelt_shapefile <- df_pixel_snowmelt_shapefile[,c("pixel_ID", "doy_snowmelt", "NDSI_threshold")]

                 #Store the pixelID of all pixels within aoi_Shapefile        
                  pixelIDs_shapefile <- unique(df_pixel_snowmelt_shapefile$pixel_ID)     
                  
               #(E): Save dataframe (is done for all NDSI-thresholds simultaneously further below)
                 #write.csv(df_pixel_snowmelt_shapefile, file=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, "_Data_Pixel_Snowmelt_polygon.csv"), quote = FALSE, row.names=FALSE)

           #(J.8): Calculate the fraction of snowcovered pixels within aoi_Shapefile per doy
                 
                #Print progress:  
                 cat("\n")
                 print("  -Step 9: Calculate the fraction of snowcovered pixels within aoi_Shapefile per doy:")
                 
               #(A) Create dataframe df_Locations_Pixel_SnowFraction_new  
                 df_Locations_Pixel_SnowFraction_new <- data.frame(LocationID=character(length(start_date_doy:end_date_doy)),
                                                                   doy=start_date_doy:end_date_doy,
                                                                   SnowFraction=numeric(length(start_date_doy:end_date_doy)))
                 
                 
                 #Remove NAs from df_pixel_snowmelt_shapefile
                 df_pixel_snowmelt_shapefile <- df_pixel_snowmelt_shapefile[df_pixel_snowmelt_shapefile$doy_snowmelt>0,]
                 
                 #Calculate the fraction of snowcovered pixels for each doy
                 for(i in 1:nrow(df_Locations_Pixel_SnowFraction_new)){
                   
                   #Select doy
                   doy_i <- df_Locations_Pixel_SnowFraction_new[i,"doy"]
                   
                   #Count fraction of pixels in df_pixel_snowmelt_shapefile which are still snowcovered on this doy
                   doy_i_SnowcoverFraction <- length(which(df_pixel_snowmelt_shapefile$doy_snowmelt > doy_i)) / length(df_pixel_snowmelt_shapefile$doy_snowmelt)
                   
                   #Store the snowfraction at this doy
                   df_Locations_Pixel_SnowFraction_new[i,"SnowFraction"] <- doy_i_SnowcoverFraction
                   
                 }
                 
                 #Add LocationID
                 df_Locations_Pixel_SnowFraction_new$LocationID <- Location_ID
                 
                 #Add NDSI_threshold as a new column
                 df_Locations_Pixel_SnowFraction_new$NDSI_threshold <- as.factor(NDSI_threshold)
                 
                 #Add coordinates to each location
                 df_Locations_Pixel_SnowFraction_new <- left_join(df_Locations_Pixel_SnowFraction_new, df_locations[,c("LON_x", "LAT_y", "LocationID")], by=c("LocationID"))
                 
                 #Add dataframe for current location to dataframe from previous iterations:
                 df_Locations_Pixel_SnowFraction <- rbind(df_Locations_Pixel_SnowFraction, df_Locations_Pixel_SnowFraction_new)

                 # #Save dataframe for current location
                 # write.csv(df_Locations_Pixel_SnowFraction_new, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_NDSI", NDSI_threshold_char, "_Data_Pixel_SnowFraction_polygon.csv"), row.names = FALSE)

           }
     
         #Store the date of snowmelt for all pixels within aoi_Shapefile (i.e. buffer zone)
          df_pixel_snowmelt_shapefile <- df_pixel_snowmelt[df_pixel_snowmelt$pixel_ID %in% pixelIDs_shapefile,]
          df_pixel_snowmelt_shapefile <- df_pixel_snowmelt_shapefile[,-which(colnames(df_pixel_snowmelt_shapefile) %in% c("lon", "lat"))]
          write.csv(df_pixel_snowmelt_shapefile, file=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_Snowmelt_polygon.csv"), quote = FALSE, row.names=FALSE)
          
         #Store NDSI timeseries for all pixels within aoi shapefile (i.e. buffer zone)
          df_pixel_ndsi_shapefile <- df_pixel_ndsi[df_pixel_ndsi$pixel_ID %in% pixelIDs_shapefile,]
          write.csv(df_pixel_ndsi_shapefile, file=paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_NDSI_polygon.csv"), quote = FALSE, row.names=FALSE)
          unlink(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Data_Pixel_NDSI_bbox.csv"))
                
         #Store plots of NDSI timeseries for all pixels within aoi_Shapefile (i.e. buffer zone)
          if(pixel_gam_plots==TRUE){
            plot_pixel_snowmelt_shapefile <- plot_pixel_snowmelt[which(pixelIDs %in% pixelIDs_shapefile)]
            plots_per_page = 25
            plot_pixel_snowmelt_shapefile <- split(plot_pixel_snowmelt_shapefile, ceiling(seq_along(plot_pixel_snowmelt_shapefile)/plots_per_page))
            pdf(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/DataLocations/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Location_", Location_i, "_Plot_Pixel_NDSI_Snowmelt_polygon.pdf"), width=20, height=16, onefile = TRUE)
            for (i in seq(length(plot_pixel_snowmelt_shapefile))) { do.call("grid.arrange", plot_pixel_snowmelt_shapefile[[i]]) }
            dev.off()
            }
          
        }
     
    }
   
 #Save combined dataframe for all locations
  
   #Store dataframe when method="avg_NDSI"
   if("avg_NDSI" %in% method){ 
    write.csv(df_Locations_BandValues, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_MeanBandvalues_polygon.csv"), row.names = FALSE)
    #df_Locations_BandValues <- read.csv(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_MeanBandvalues_polygon.csv"), header = T)
    }
    
   #Store dataframe when method="snowfraction"
   if("snowfraction" %in% method){ 
    write.csv(df_Locations_SnowFraction, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_SnowFraction_polygon.csv"), row.names = FALSE)
    #df_Locations_SnowFraction <- read.csv(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_SnowFraction_polygon.csv"), header=T)
    }
    
   #Store dataframe when method="pixel_gam"
   if("pixel_gam" %in% method){ 
    write.csv(df_Locations_Pixel_SnowFraction, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_Pixel_SnowFraction_polygon.csv"), row.names = FALSE)
    #df_Locations_Pixel_SnowFraction <- read.csv(paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_Pixel_SnowFraction_polygon.csv"), header=T)
    }

  
##################################################################################################################################

#VI: Fit a Generalized Additive Model (GAM) through the Mean NDSI and Snowfraction data for each location

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
      
###############################################################################################################################################################

 #(I - avg_NDSI) - Fit a Generalized Additive Model (GAM) through the NDSI data for each Location
  
  if("avg_NDSI" %in% method){

   #(A) Specify at which NDSI value(s), a location is considered snow-free
    NDSI_threshold_vector = NDSI_threshold_vector

   #(B) Create an empty dataframe
    df_Locations_NDSI <- data.frame(NDSI=numeric(),
                                    Date=factor(),
                                    doy=numeric(),
                                    LocationID=factor(),
                                    outliers=logical())

    df_Locations_NDSI_predictions <- data.frame(LocationID=character(),
                                                doy=numeric(),
                                                NDSI_gam_predict=numeric(),
                                                stringsAsFactors=FALSE)

   #(C) Loop through all Locations and fit a separate gam with sequential outlier removal to the location specific mean NDSI data
    for(i in unique(df_Locations_BandValues$LocationID)){

      #For debugging
      #i=unique(df_Locations_BandValues$LocationID)[1]

      #Select Location-specific subset of data:
       df_Location_NDSI_new <- df_Locations_BandValues[df_Locations_BandValues$LocationID==i & !is.na(df_Locations_BandValues$NDSI),
                                           c("NDSI", "Date", "doy", "LocationID")]

      #Fit a gam through the Location-specific mean_NDSI ~ doy data and employ sequential outlier removal
       df_Location_NDSI_new <- f_gam_SeqRemOutliers(data=df_Location_NDSI_new, y="NDSI", x="doy", outlier_removal=outlier_removal,
                                           outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                           default_k=gam_k_outlier)

      #Sort df_Location_NDSI_new by doy:
       df_Location_NDSI_new <- df_Location_NDSI_new[order(df_Location_NDSI_new$doy),]

      #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
       df_Locations_NDSI <- rbind(df_Locations_NDSI, df_Location_NDSI_new)

      #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

       #Refit GAM through data
        index <- which(df_Location_NDSI_new$outliers==FALSE)
        mod_gam <- with(df_Location_NDSI_new[index,], mgcv::gam(NDSI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

       #Use gam to make predictions on a more detailed (1-day) day of year interval
        df_Locations_NDSI_predictions_new <- data.frame(LocationID=i, doy=seq(min(df_Location_NDSI_new$doy), max(df_Location_NDSI_new$doy), 1))
        df_Locations_NDSI_predictions_new$NDSI_gam_predict <- stats::predict(mod_gam, newdata=df_Locations_NDSI_predictions_new, type="response")
        df_Locations_NDSI_predictions_new <- df_Locations_NDSI_predictions_new[order(df_Locations_NDSI_predictions_new$doy),]

       #Add predictions to df_Locations_NDSI_predictions dataframe:
        df_Locations_NDSI_predictions <- rbind(df_Locations_NDSI_predictions, df_Locations_NDSI_predictions_new)

        }

     #Change column LocationID to a factor:
      df_Locations_NDSI$LocationID <- as.factor(as.character(df_Locations_NDSI$LocationID))
      df_Locations_NDSI_predictions$LocationID <- as.factor(as.character(df_Locations_NDSI_predictions$LocationID))

     #Save dataframe with GAM fits for NDSI
      #write.csv(df_Locations_NDSI, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_NDSI_GAM_polygon.csv"), row.names = FALSE)
      write.csv(df_Locations_NDSI_predictions, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_GAM_Predictions_NDSI_polygon.csv"), row.names = FALSE)
      
   #(D) Plot the raw NDSI datapoints and gam predictions for each Location:

     # #Plot NDSI and model predictions for all locations in a single plot
     #  p_Locations_NDSI = ggplot()+
     #    geom_point(data=df_Locations_NDSI, aes(x=doy, y=NDSI, fill=LocationID, col=LocationID))+
     #    geom_line(data=df_Locations_NDSI_predictions, aes(x=doy, y=NDSI_gam_predict, col=LocationID)) +
     #    xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
     #    ylab("NDSI") +
     #    theme_tom()
     # 
     #  ggsave(plot=p_Locations_NDSI, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_NDSI_polygon.pdf"), width=10, height = 8)

     # Plot NDSI and model predictions in a separate plot per location
      p_Locations_NDSI_grid = ggplot()+
        geom_point(data=df_Locations_NDSI[df_Locations_NDSI$outliers==FALSE,], aes(x=doy, y=NDSI))+
        geom_point(data=df_Locations_NDSI[df_Locations_NDSI$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
        geom_line(data=df_Locations_NDSI_predictions, aes(x=doy, y=NDSI_gam_predict), col="#1620de", linewidth=1.25)+
        geom_vline(xintercept = 150, colour="grey", lty=2)+
        geom_hline(yintercept = NDSI_threshold_vector, colour="grey", lty=2)+
        facet_wrap(~LocationID, ncol=ceiling(length(unique(df_Locations_NDSI$LocationID))^0.5))+
        xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
        ylab("NDSI") +
        theme_tom()

     #  ggsave(plot=p_Locations_NDSI_grid, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_NDSI_polygon.pdf"), width=12, height = 10)

   #(E) Calculate at which day of year the average NDSI is equal to NDSI_threshold_vector for each location using predictions from mod_gam

     #Create an empty dataframe
      df_Snowmelt_NDSI_Locations <- data.frame(LocationID=character(),
                                          NDSI_threshold=character(),
                                          doy=numeric(),
                                          stringsAsFactors=FALSE)

      #Setup parallel processing
      numCores <- detectCores()
      cl <- makePSOCKcluster(numCores)
      registerDoSNOW(cl)
      
     #Extract the day of snowmelt for each location for all levels of NDSI_threshold_vector
      for(LocationID in unique(df_Locations_NDSI_predictions$LocationID)){

          #For debugging
          #LocationID=unique(df_Locations_NDSI_predictions$LocationID)[1] #for debugging

          #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in NDSI_threshold_vector
          results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                      pixelIDs_split=LocationID, #Current location (input during loop)
                                                      df_pixel_y=df_Locations_NDSI_predictions, #dataframe containing GAM predictions
                                                      pixel_ID_column="LocationID", #Grouping column
                                                      y="NDSI_gam_predict", #response variable in GAM
                                                      x="doy", #predictor variable in GAM
                                                      pixel_gam_plots=FALSE, #Should GAM plots be created
                                                      y_threshold=NDSI_threshold_vector) #Which threshold values for 'y' should be calculated
        
          #Store dates of snowmelt per Location
          df_Snowmelt_NDSI_Locations_new <- lapply(results, "[[", 1)
          df_Snowmelt_NDSI_Locations_new <- as.data.frame(do.call(rbind, df_Snowmelt_NDSI_Locations_new))
          colnames(df_Snowmelt_NDSI_Locations_new)[colnames(df_Snowmelt_NDSI_Locations_new)=="pixel_ID"] <- "LocationID"
          colnames(df_Snowmelt_NDSI_Locations_new)[colnames(df_Snowmelt_NDSI_Locations_new)=="x_threshold"] <- "doy"
          colnames(df_Snowmelt_NDSI_Locations_new)[colnames(df_Snowmelt_NDSI_Locations_new)=="y_threshold"] <- "NDSI_threshold"
         
          #Add snowmelt data for LocationID to general dataframe
          df_Snowmelt_NDSI_Locations <- rbind(df_Snowmelt_NDSI_Locations, df_Snowmelt_NDSI_Locations_new)
         
      }
      
     #Turn parallel processing off and run sequentially again after this point
      stopCluster(cl)
      registerDoSEQ()
         
     #Add coordinates to each location
      df_Snowmelt_NDSI_Locations <- left_join(df_Snowmelt_NDSI_Locations, df_locations, by=c("LocationID"))  
      
     #Save dates of snowmelt per Location as a .csv file
      write.csv(df_Snowmelt_NDSI_Locations, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Snowmelt_NDSI_polygon.csv"), row.names = FALSE)
      
     #Add dates of snowmelt to the plot 'p_Locations_NDSI_grid'
      p_Locations_NDSI_Snowmelt_grid <- p_Locations_NDSI_grid +
       geom_point(data=df_Snowmelt_NDSI_Locations[!is.na(df_Snowmelt_NDSI_Locations$doy),], aes(x=doy, y=NDSI_threshold), col="red", size=3)
      
      ggsave(plot=p_Locations_NDSI_Snowmelt_grid, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_Snowmelt_NDSI_polygon.pdf"), width=12, height = 10)
    
      }
  
 #(II - snowfraction) - Fit a Generalized Additive Model (GAM) through the SnowFraction data
  
  if("snowfraction" %in% method){

   #(A) Specify for which snowfraction the corresponding date needs to be extracted
     Snowfraction_threshold_vector = Snowfraction_threshold_vector

   #(B) Create an empty dataframe
     df_Locations_SnowFraction_GAM <- data.frame(NDSI_threshold=character(),
                                                 SnowFraction=numeric(),
                                                 Date=factor(),
                                                 doy=numeric(),
                                                 LocationID=factor(),
                                                 outliers=logical())

     df_Locations_SnowFraction_GAM_predictions <- data.frame(LocationID=character(),
                                                             NDSI_threshold=character(),
                                                             doy=numeric(),
                                                             SnowFraction_gam_predict=numeric(),
                                                             stringsAsFactors=FALSE)

   #(C) Loop through all Locations and fit a separate GAM with sequential outlier removal to the location specific SnowFraction data
    for(i in unique(df_Locations_SnowFraction$LocationID)){

      #For debugging
      #i=unique(df_Locations_SnowFraction$LocationID)[1]

      #Loop through all NDSI_thresholds
      for(j in unique(df_Locations_SnowFraction$NDSI_threshold)){
        
        #For debugging
        #j=unique(df_Locations_SnowFraction$NDSI_threshold)[1]
        
        #Select Location-specific subset of data:
        df_Location_SnowFraction_GAM_new <- df_Locations_SnowFraction[df_Locations_SnowFraction$LocationID==i & 
                                                                df_Locations_SnowFraction$NDSI_threshold==j & 
                                                                !is.na(df_Locations_SnowFraction$SnowFraction),
                                                              c("NDSI_threshold", "SnowFraction", "Date", "doy", "LocationID")]
        
        #Fit a gam through the Location-specific mean_NDSI ~ doy data and employ sequential outlier removal
        df_Location_SnowFraction_GAM_new <- f_gam_SeqRemOutliers(data=df_Location_SnowFraction_GAM_new, y="SnowFraction", x="doy", outlier_removal=outlier_removal,
                                                         outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                         default_k=gam_k_outlier)
        
        #Sort df_Location_SnowFraction_GAM_new by doy:
        df_Location_SnowFraction_GAM_new <- df_Location_SnowFraction_GAM_new[order(df_Location_SnowFraction_GAM_new$doy),]
        
        #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
        df_Locations_SnowFraction_GAM <- rbind(df_Locations_SnowFraction_GAM, df_Location_SnowFraction_GAM_new)
        
        #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves
        
          #Refit GAM through data
          index <- which(df_Location_SnowFraction_GAM_new$outliers==FALSE)
          mod_gam <- with(df_Location_SnowFraction_GAM_new[index,], mgcv::gam(SnowFraction ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))
          
          #Use gam to make predictions on a more detailed (1-day) day of year interval
          df_Locations_SnowFraction_GAM_predictions_new <- data.frame(LocationID=i, NDSI_threshold=j, doy=seq(min(df_Location_SnowFraction_GAM_new$doy), max(df_Location_SnowFraction_GAM_new$doy), 1))
          df_Locations_SnowFraction_GAM_predictions_new$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=df_Locations_SnowFraction_GAM_predictions_new, type="response")
          df_Locations_SnowFraction_GAM_predictions_new <- df_Locations_SnowFraction_GAM_predictions_new[order(df_Locations_SnowFraction_GAM_predictions_new$doy),]
          df_Locations_SnowFraction_GAM_predictions_new$year <- year_ID
          
          #Add predictions to df_Locations_SnowFraction_GAM_predictions dataframe:
          df_Locations_SnowFraction_GAM_predictions <- rbind(df_Locations_SnowFraction_GAM_predictions, df_Locations_SnowFraction_GAM_predictions_new)
        
      }

    }

      #Change column LocationID to a factor:
       df_Locations_SnowFraction_GAM$LocationID <- as.factor(as.character(df_Locations_SnowFraction_GAM$LocationID))
       df_Locations_SnowFraction_GAM_predictions$LocationID <- as.factor(as.character(df_Locations_SnowFraction_GAM_predictions$LocationID))

      #Save dataframe with GAM fits for the SnowFraction data
       #write.csv(df_Locations_SnowFraction_GAM, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_SnowFraction_GAM_polygon.csv"), row.names = FALSE)
       write.csv(df_Locations_SnowFraction_GAM_predictions, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_GAM_Predictions_SnowFraction_polygon.csv"), row.names = FALSE)
       
   #(D) Plot the raw SnowFraction datapoints and gam predictions for each Location:

    #Plot SnowFraction and model predictions for all locations in a single plot
      p_Locations_SnowFraction = ggplot()+
       geom_point(data=df_Locations_SnowFraction_GAM, aes(x=doy, y=SnowFraction, fill=LocationID, col=LocationID))+
       geom_line(data=df_Locations_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict, col=LocationID)) +
       xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
       ylab(paste0("Snowcover fraction per location in ", year_ID)) +
       facet_wrap(~NDSI_threshold, ncol=3)+
       theme_tom()
     
    #Plot SnowFraction and model predictions in a separate plot per location
      p_Locations_SnowFraction_grid = ggplot()+
       geom_point(data=df_Locations_SnowFraction_GAM[df_Locations_SnowFraction_GAM$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
       geom_point(data=df_Locations_SnowFraction_GAM[df_Locations_SnowFraction_GAM$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
       geom_line(data=df_Locations_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict), col="#1620de", lwd=1)+
       geom_vline(xintercept = 150, colour="grey", lty=2)+
       geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
       facet_wrap(~LocationID + NDSI_threshold, ncol=ceiling(length(unique(df_Locations_SnowFraction_GAM$LocationID))^0.5))+
       xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
       ylab(paste0("Snowcover fraction per location in ", year_ID)) +
       theme_tom()

      # ggsave(plot=p_Locations_SnowFraction_grid, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_SnowFraction_polygon.pdf"), width=12, height = 10)

   #(E) Calculate at which day of year the SnowFraction is equal to Snowfraction_threshold_vector for each location using predictions from mod_gam

     #Create an empty dataframe
      df_SnowFraction_Locations <- data.frame(NDSI_threshold=character(),
                                              Snowfraction_threshold=character(),
                                              doy=numeric(),
                                              LocationID=character(),
                                              stringsAsFactors=FALSE)

      #Setup parallel processing
      numCores <- detectCores()
      cl <- makePSOCKcluster(numCores)
      registerDoSNOW(cl)
      
     #Loop through all locations in the df_Locations_SnowFraction_GAM_predictions dataframe and extract doy for all Snowfraction levels:
      for(i in unique(df_Locations_SnowFraction_GAM_predictions$LocationID)){

       #For debugging
       #i=unique(df_Locations_SnowFraction_GAM_predictions$LocationID)[1] #for debugging

       #Select dataset with GAM predictions for current Location
       df_SnowFraction_gam <- df_Locations_SnowFraction_GAM_predictions[df_Locations_SnowFraction_GAM_predictions$LocationID==i & !is.na(df_Locations_SnowFraction_GAM_predictions$SnowFraction),]
        
       #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector
       results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                   pixelIDs_split = list(NDSI_threshold_vector), #levels of NDSI_threshold (input needs to be a list)
                                                   df_pixel_y = df_SnowFraction_gam, #dataframe containing GAM predictions
                                                   pixel_ID_column="NDSI_threshold", #Grouping column
                                                   y="SnowFraction_gam_predict", #response variable in GAM
                                                   x="doy", #predictor variable in GAM
                                                   pixel_gam_plots = FALSE, #Should GAM plots be created
                                                   y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated
          
       #Store dates of snowmelt per Location
       df_SnowFraction_Locations_new <- results[[1]]
       df_SnowFraction_Locations_new <- as.data.frame(do.call(rbind, df_SnowFraction_Locations_new))
       colnames(df_SnowFraction_Locations_new)[colnames(df_SnowFraction_Locations_new)=="pixel_ID"] <- "NDSI_threshold"
       colnames(df_SnowFraction_Locations_new)[colnames(df_SnowFraction_Locations_new)=="x_threshold"] <- "doy"
       colnames(df_SnowFraction_Locations_new)[colnames(df_SnowFraction_Locations_new)=="y_threshold"] <- "Snowfraction_threshold"
       df_SnowFraction_Locations_new$LocationID <- i
       
       #Add snowmelt data for LocationID to general dataframe
       df_SnowFraction_Locations <- rbind(df_SnowFraction_Locations, df_SnowFraction_Locations_new)
       
      }
      
      #Turn parallel processing off and run sequentially again after this point
      stopCluster(cl)
      registerDoSEQ()
       
      #Add coordinates to each location
      df_SnowFraction_Locations <- left_join(df_SnowFraction_Locations, df_locations, by=c("LocationID"))  
      
      #Save dates of snowmelt per Location per SnowFraction threshold as a .csv file
      write.csv(df_SnowFraction_Locations, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Snowmelt_Snowfraction_polygon.csv"), row.names = FALSE)
      
      #Add dates of snowmelt to the plot 'p_Locations_NDSI_grid'
      p_Locations_SnowFraction_Snowmelt_grid <- p_Locations_SnowFraction_grid +
        geom_point(data=df_SnowFraction_Locations[!is.na(df_SnowFraction_Locations$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)
      
      ggsave(plot=p_Locations_SnowFraction_Snowmelt_grid, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_Snowmelt_Snowfraction_polygon.pdf"), width=12, height = 10)
      
   }
     
 #(III - pixel_gam) - Fit a Generalized Additive Model (GAM) through the pixel-level SnowFraction data
  
  if("pixel_gam" %in% method){

   #(A) Specify for which snowfraction the corresponding date needs to be extracted
     Snowfraction_threshold_vector = Snowfraction_threshold_vector

   #(B) Create an empty dataframe
     df_Locations_Pixel_SnowFraction_GAM <- data.frame(NDSI_threshold=character(),
                                                       SnowFraction=numeric(),
                                                       Date=factor(),
                                                       doy=numeric(),
                                                       LocationID=factor(),
                                                       outliers=logical())

     df_Locations_Pixel_SnowFraction_GAM_predictions <- data.frame(LocationID=character(),
                                                                   NDSI_threshold=character(),
                                                                   doy=numeric(),
                                                                   SnowFraction_gam_predict=numeric(),
                                                                   stringsAsFactors=FALSE)

   #(C) Loop through all Locations and fit a separate gam with sequential outlier removal to the location specific pixel-level SnowFraction data
    for(i in unique(df_Locations_Pixel_SnowFraction$LocationID)){

      #For debugging
      #i=unique(df_Locations_Pixel_SnowFraction$LocationID)[1]

      #Loop through all NDSI_thresholds
      for(j in unique(df_Locations_Pixel_SnowFraction$NDSI_threshold)){
    
        #For debugging
        #j=unique(df_Locations_Pixel_SnowFraction$NDSI_threshold)[1]
          
      #Select Location-specific subset of data:
       df_Location_Pixel_SnowFraction <- df_Locations_Pixel_SnowFraction[df_Locations_Pixel_SnowFraction$LocationID==i & 
                                                                           df_Locations_Pixel_SnowFraction$NDSI_threshold==j & 
                                                                           !is.na(df_Locations_Pixel_SnowFraction$SnowFraction),
                                                                         c("NDSI_threshold", "SnowFraction", "doy", "LocationID")]

       #Fit a gam through the Location-specific mean_NDSI ~ doy data and employ sequential outlier removal
       df_Location_Pixel_SnowFraction <- f_gam_SeqRemOutliers(data=df_Location_Pixel_SnowFraction, y="SnowFraction", x="doy", outlier_removal=outlier_removal,
                                                        outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                        default_k=gam_k_outlier)

       #Sort df_Location_Pixel_SnowFraction by doy:
       df_Location_Pixel_SnowFraction <- df_Location_Pixel_SnowFraction[order(df_Location_Pixel_SnowFraction$doy),]

       #Bind the Location-specific dataframe with GAM estimates to the dataframe containing all dataframes from previous iterations
       df_Locations_Pixel_SnowFraction_GAM <- rbind(df_Locations_Pixel_SnowFraction_GAM, df_Location_Pixel_SnowFraction)

      #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

        #Refit GAM through data
         index <- which(df_Location_Pixel_SnowFraction$outliers==FALSE)
         mod_gam <- with(df_Location_Pixel_SnowFraction[index,], mgcv::gam(SnowFraction ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

        #Use gam to make predictions on a more detailed (1-day) day of year interval
         aoi_Pixel_SnowFraction_predicted <- data.frame(LocationID=i, NDSI_threshold=j, doy=seq(min(df_Location_Pixel_SnowFraction$doy), max(df_Location_Pixel_SnowFraction$doy), 1))
         aoi_Pixel_SnowFraction_predicted$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=aoi_Pixel_SnowFraction_predicted, type="response")
         aoi_Pixel_SnowFraction_predicted <- aoi_Pixel_SnowFraction_predicted[order(aoi_Pixel_SnowFraction_predicted$doy),]
         aoi_Pixel_SnowFraction_predicted$year <- year_ID

      #Add predictions to df_Locations_Pixel_SnowFraction_GAM_predictions dataframe:
       df_Locations_Pixel_SnowFraction_GAM_predictions <- rbind(df_Locations_Pixel_SnowFraction_GAM_predictions, aoi_Pixel_SnowFraction_predicted)

      }
    
    }
    
    #Change column LocationID to a factor:
    df_Locations_Pixel_SnowFraction_GAM$LocationID <- as.factor(as.character(df_Locations_Pixel_SnowFraction_GAM$LocationID))
    df_Locations_Pixel_SnowFraction_GAM_predictions$LocationID <- as.factor(as.character(df_Locations_Pixel_SnowFraction_GAM_predictions$LocationID))

    #Save dataframe with GAM fits for pixel level Snowfraction data
    #write.csv(df_Locations_Pixel_SnowFraction_GAM, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Data_Pixel_SnowFraction_GAM_polygon.csv"), row.names = FALSE)
    write.csv(df_Locations_Pixel_SnowFraction_GAM_predictions, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_GAM_Predictions_Pixel_SnowFraction_polygon.csv"), row.names = FALSE)
    
   #(D) Plot the raw pixel-level SnowFraction datapoints and gam predictions for each Location:

    # #Plot SnowFraction and model predictions for all locations in a single plot
    #   p_Locations_Pixel_SnowFraction = ggplot()+
    #    geom_point(data=df_Locations_Pixel_SnowFraction_GAM, aes(x=doy, y=SnowFraction, fill=LocationID, col=LocationID))+
    #    geom_line(data=df_Locations_Pixel_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict, col=LocationID)) +
    #    xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
    #    ylab(paste0("Snowcover fraction per location in ", year_ID)) +
    #    theme_tom()
    # 
    #   ggsave(plot=p_Locations_Pixel_SnowFraction, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_Pixel_SnowFraction_polygon.pdf"), width=10, height = 8)

    #Plot SnowFraction and model predictions in a separate plot per location
      p_Locations_Pixel_SnowFraction_grid = ggplot()+
       geom_point(data=df_Locations_Pixel_SnowFraction_GAM[df_Locations_Pixel_SnowFraction_GAM$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
       geom_point(data=df_Locations_Pixel_SnowFraction_GAM[df_Locations_Pixel_SnowFraction_GAM$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
       geom_line(data=df_Locations_Pixel_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict), col="#1620de", lwd=1.25)+
       geom_vline(xintercept = 150, colour="grey", lty=2)+
       geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
       facet_wrap(~LocationID + NDSI_threshold, ncol=ceiling(length(unique(df_Locations_Pixel_SnowFraction_GAM$LocationID))^0.5))+
       xlab(paste0("Day of year (starting at 01-01-", year_ID,")")) +
       ylab(paste0("Snowcover fraction per location in ", year_ID)) +
       theme_tom()

      # ggsave(plot=p_Locations_Pixel_SnowFraction_grid, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_Pixel_SnowFraction_polygon.pdf"), width=12, height = 10)

   #(E) Calculate at which day of year the pixel-level SnowFraction is equal to Snowfraction_threshold_vector for each location using predictions from mod_gam

     #Create an empty dataframe
      df_Pixel_SnowFraction_Locations <- data.frame(NDSI_threshold=character(),
                                                    Snowfraction_threshold=character(),
                                                    doy=numeric(),
                                                    LocationID=character(),
                                                    stringsAsFactors=FALSE)
      
      #Setup parallel processing
      numCores <- detectCores()
      cl <- makePSOCKcluster(numCores)
      registerDoSNOW(cl)
      
     #Loop through all locations in the df_Locations_Pixel_SnowFraction_GAM_predictions dataframe:
      for(i in unique(df_Locations_Pixel_SnowFraction_GAM_predictions$LocationID)){

       #For debugging
       #i=unique(df_Locations_Pixel_SnowFraction_GAM_predictions$LocationID)[1] #for debugging

       #Select dataset with GAM predictions for current Location
       df_Pixel_SnowFraction_gam <- df_Locations_Pixel_SnowFraction_GAM_predictions[df_Locations_Pixel_SnowFraction_GAM_predictions$LocationID==i & !is.na(df_Locations_Pixel_SnowFraction_GAM_predictions$SnowFraction),]
        
       #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector
       results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                   pixelIDs_split = list(NDSI_threshold_vector), #levels of NDSI_threshold (input needs to be a list)
                                                   df_pixel_y = df_Pixel_SnowFraction_gam, #dataframe containing GAM predictions
                                                   pixel_ID_column="NDSI_threshold", #Grouping column
                                                   y="SnowFraction_gam_predict", #response variable in GAM
                                                   x="doy", #predictor variable in GAM
                                                   pixel_gam_plots = FALSE, #Should GAM plots be created
                                                   y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated
        
       #Store dates of snowmelt per Location
       df_Pixel_SnowFraction_Locations_new <- results[[1]]
       df_Pixel_SnowFraction_Locations_new <- as.data.frame(do.call(rbind, df_Pixel_SnowFraction_Locations_new))
       colnames(df_Pixel_SnowFraction_Locations_new)[colnames(df_Pixel_SnowFraction_Locations_new)=="pixel_ID"] <- "NDSI_threshold"
       colnames(df_Pixel_SnowFraction_Locations_new)[colnames(df_Pixel_SnowFraction_Locations_new)=="x_threshold"] <- "doy"
       colnames(df_Pixel_SnowFraction_Locations_new)[colnames(df_Pixel_SnowFraction_Locations_new)=="y_threshold"] <- "Snowfraction_threshold"
       df_Pixel_SnowFraction_Locations_new$LocationID <- i 
        
       #Add snowmelt data for LocationID to general dataframe
       df_Pixel_SnowFraction_Locations <- rbind(df_Pixel_SnowFraction_Locations, df_Pixel_SnowFraction_Locations_new)
        
      }
      
     #Turn parallel processing off and run sequentially again after this point
      stopCluster(cl)
      registerDoSEQ()
        
     #Add coordinates to each location
      df_Pixel_SnowFraction_Locations <- left_join(df_Pixel_SnowFraction_Locations, df_locations, by=c("LocationID"))     
        
     #Save dates of snowmelt per Location per SnowFraction threshold as a .csv file
      write.csv(df_Pixel_SnowFraction_Locations, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Snowmelt_Pixel_Snowfraction_polygon.csv"), row.names = FALSE)
      
     #Add dates of snowmelt to the plot 'p_Locations_NDSI_grid'
      p_Locations_Pixel_SnowFraction_Snowmelt_grid <- p_Locations_Pixel_SnowFraction_grid +
        geom_point(data=df_Pixel_SnowFraction_Locations[!is.na(df_Pixel_SnowFraction_Locations$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)
      
      ggsave(plot=p_Locations_Pixel_SnowFraction_Snowmelt_grid, paste0(here(), "/Output/MODIS/02_Points_Snowmelt/", timestamp, "_", data_ID, "_Buffer", Buffer_radius_m, "_Res", resolution, "_Locations_Plot_Snowmelt_Pixel_Snowfraction_polygon.pdf"), width=12, height = 10)
      
  }

   
##########################################################################################################################################################################
  
#The End   
  
##########################################################################################################################################################################