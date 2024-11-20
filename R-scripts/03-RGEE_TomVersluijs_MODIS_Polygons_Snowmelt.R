##################################################################################################################################

#The date of snow melt is calculated based on MODIS data for a single **`polygon`**, or **`multi-polygon** specified using a
#shapefile. This shapefile can be created using e.g. QGIS (see manual). The user can specify whether clouds and permanent
#water bodies need to be masked. Snow melt is analysed within each polygon based on one of the following methods (specified
#by the user by setting the parameter 'method'):

  #(1) 'avg_NDSI': Calculate the average NDSI value over time within each polygon, fit a GAM through these data and calculate
  #     when this model passes the specified NDSI threshold representing the moment of snow melt. In addition, time series of
  #     the average NDVI and NDMI are extracted within each polygon. Also, time series of the average Fractional Snow Cover
  #     (FSC, a within-pixel estimate of the fraction of snow cover) are extracted within each polygon based on the formulas
  #     specified in Gascoin et al 2020 and Aalstad et al 2020.
  #(2) 'snowfraction': Calculate the fraction of pixels within each polygon over time where NDSI > 'NDSI_threshold', fit
  #     a GAM through these data and extract the moment when this model passes a user-specified 'Snowfraction_threshold'.

#Note that snow melt is not calculated based on pixel-level GAM fits (i.e. 'pixel_gam' method is not implemented). This approach
#is instead implemented for a single polygon in script "01-RGEE_TomVersluijs_MODIS_Pixels_Snowmelt.R" and involves fitting of
#GAMS through NDSI data per pixel and extracting the moment this GAM passes a user-defined NDSI-threshold. This results in a
#pixel-level map of the date of snow melt. Script "10-RGEE_TomVersluijs_ExtractSnowFraction.R" can then be used to extract
#time series of the fraction of snow cover for points/polygons of interest from this map.

#Copyright Tom Versluijs 2024-11-19. Do not use this code without permission. Contact information: tom.versluijs@gmail.com

#Before running this script make sure to install RGEE according to the instructions in script "00-RGEE_TomVersluijs_Installation.R".
#Note that a GoogleDrive is required. Important: make sure to run this script from within the "RGEE_Snowmelt.Rproj" project file.

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
     p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr, tidyr, foreach, parallel, doSNOW, gridExtra, rgeeExtra, magick)

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

    #(4): Initialize earth engine
     rgee::ee_Initialize(user = "tom.versluijs@gmail.com", drive = TRUE)

##################################################################################################################################

#II: Specify parameters of interest

##################################################################################################################################

 #(5): Manually specify parameters of interest

   #(a): MODIS satellite

     #MODIS dataset
     MODIS_dataset <- "MODIS/061/MOD09GA"

     #Specify resolution of images in meters
     resolution=500 #default maximum resolution for MODIS = 500m

   #(b) Areas of interest

     #Name of Shapefile(s) containing either a single polygon, or a multi-polygon (located in '/Input/Shapefiles' folder)
	   #Follow the guide "Manual_CreateShapefilePolygons.pdf" when creating these shapefiles. Specify multiple using c().
     shapefile_vector <- "ZAC_TenEqualSizeVoronoiPolygons_EPSG4326.shp"

     #Coordinate reference system used for calculations (EPSG:4326 is recommended)
     crs <- "EPSG:4326"

   #(c) Dates

     #Specify the year(s) of interest (specify multiple using c())
     year_ID_vector <- as.character(2000:2023)

     #Date range of all images considered for analysis (format "-mm-dd")
     start_month_day <- "-03-15" #choose date (well) before the first snow melt occurs within the study area
     end_month_day <- "-09-15" #choose date (well) after last snow melt occurs within the study area

   #(d) Snow detection

     #NDSI threshold for snow detection (specify multiple using c())
     NDSI_threshold_vector = c(0.3, 0.4, 0.5)

     #Define the snowcover fraction for which the date of its occurrence will be calculated (specify multiple using c())
     Snowfraction_threshold_vector = c(0.25, 0.5, 0.75)

     #Define the preferred method(s) for the snowmelt analysis
     method=c("avg_NDSI", "snowfraction") #either "avg_NDSI", "snowfraction", or a combination using c()
     #(1) "avg_NDSI":     Calculates the average FSC, NDSI, NDVI and NDMI values within each polygon over time and extracts the moment
     #                    when the NDSI value is equal to 'NDSI_threshold', and FSC is equal to 'Snowfraction_threshold'
     #(2) "snowfraction": Counts the fraction of pixels within each polygon with NDSI > 'NDSI_threshold' over time and extracts the moment
     #                    when this fraction is equal to 'Snowfraction_threshold'.

     #Define model output for method 'avg_NDSI' (either "FSC_Gascoin2020", "FSC_Aalstad2020", "NDSI", "NDVI", and/or "NDMI")
     method_output <- c("FSC_Gascoin2020", "FSC_Aalstad2020", "NDSI", "NDVI", "NDMI")

   #(e): Cloud masking

     #Should clouds be masked from the analysis (default=TRUE).
     mask_clouds=TRUE

     #Define which MODIS cloud masking algorithm has been used ("PGE11", "MOD35", or "Combined)
     MODIS_cloud_masking_algorithm = 'MOD35' #default is MOD35

     #Maximum fraction of cloud cover allowed in each image
     max_cloud_fraction=0.75 #0.75 equals 75% cloud cover

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
     #Larger values result in a more precise GAM-fit, but might result in overfitting.

   #(h): Counts of unmasked pixels

     #Should counts of the number of unmasked pixels per day of year within the shapefile area be conducted (increases computation time)
     pixel_counts=TRUE
     
   #(i): GIF animations
     
     #Should GIF animations be constructed (increases computation time)
     gif_output=FALSE
     
     #Maximum dimension (in pixels) of GIF
     gif_max_pixels=400

##################################################################################################################################

 #Store the names of all parameters        
 params_output <- setdiff(ls(), lsf.str())             
     
 #Run the analysis for all shapefiles
 for(shapefile in shapefile_vector){

  #For debugging
  #shapefile <- shapefile_vector[1]
  
  #Print progress
  cat("\n") ; cat("\n")
  print("##########################################################################################################################")
  print(paste0("STARTING ANALYSIS FOR SHAPEFILE: ", shapefile))
  print("##########################################################################################################################")
  cat("\n")

  #Create empty lists to store dataframes and plots per year

    #Pixelcounts
    list_df_pixelcount <- vector('list', length(year_ID_vector))
    list_p_pixelcounts <- vector('list', length(year_ID_vector))

    #Raw data
    list_df_Polygons_BandValues <- vector('list', length(year_ID_vector))
    list_df_Polygons_SnowFraction <- vector('list', length(year_ID_vector))

    #Snowfraction
    list_df_Polygons_SnowFraction_GAM_predictions <- vector('list', length(year_ID_vector))
    list_df_Polygon_Snowfraction <- vector('list', length(year_ID_vector))
    list_p_snowfraction <- vector('list', length(year_ID_vector))

    #FSC Gascoin2020
    list_df_Polygons_FSC_Gascoin2020_GAM_predictions <- vector('list', length(year_ID_vector))
    list_df_Polygon_FSC_Gascoin2020 <- vector('list', length(year_ID_vector))
    list_p_fsc_gascoin2020 <- vector('list', length(year_ID_vector))

    #FSC Aalstad2020
    list_df_Polygons_FSC_Aalstad2020_GAM_predictions <- vector('list', length(year_ID_vector))
    list_df_Polygon_FSC_Aalstad2020 <- vector('list', length(year_ID_vector))
    list_p_fsc_aalstad2020 <- vector('list', length(year_ID_vector))

    #NDSI
    list_df_Polygons_NDSI_GAM_predictions <- vector('list', length(year_ID_vector))
    list_df_Polygon_NDSI <- vector('list', length(year_ID_vector))
    list_p_ndsi <- vector('list', length(year_ID_vector))

    #NDVI
    list_df_Polygons_NDVI_GAM_predictions <- vector('list', length(year_ID_vector))
    list_p_ndvi <- vector('list', length(year_ID_vector))

    #NDMI
    list_df_Polygons_NDMI_GAM_predictions <- vector('list', length(year_ID_vector))
    list_p_ndmi <- vector('list', length(year_ID_vector))

  #Run the analysis for all years
  for(year_ID in year_ID_vector){

    #year_ID <- year_ID_vector[1]

    #Print progress
    cat("\n") ; cat("\n")
    print("##########################################################################################################################")
    print(paste0("STARTING ANALYSIS FOR YEAR: ", year_ID))
    print("##########################################################################################################################")
    cat("\n")

    #Define date range of all images considered for analysis
    start_date <- paste0(year_ID, start_month_day) #choose date (well) before the first snowmelt occurs within the study site
    end_date <- paste0(year_ID, end_month_day) #choose date (well) after last snowmelt occurs within the study site

    #Define area name (used as prefix in output files)
    area_name <- toupper(substr(shapefile, 1, 3))

##################################################################################################################################

#III: Define some additional parameters (automated)

##################################################################################################################################

   #(6): Automatically define some additional parameters

     #Create a unique data_ID
     if(nchar(area_name)>3){area_name <- substr(area_name, start = 1, stop = 3)}
     data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
     data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)

     #Create a timestamp variable
     timestamp <- format(Sys.time(), "%Y%m%d%H%M%S")

     #First and last day of year in dataset
     start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
     end_date_doy <- as.numeric(strftime(end_date, format = "%j"))

     #Create output folder
     if(dir.exists(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt"))==FALSE){dir.create(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt"), recursive = TRUE)}

     #Save all parameters and their values in the environment to a text file (only once, because parameters are identical between runs for different years)
     if(year_ID == year_ID_vector[1]){
      file_conn <- file(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", paste0(area_name, "_MODIS_", MODIS_cloud_masking_algorithm), "_Parameters.txt"), "w")
      for (obj in sort(c(params_output, "shapefile", "year_ID", "start_date", "end_date", "start_date_doy", "end_date_doy", "area_name", "data_ID", "timestamp"))) {cat(paste(obj, "=", get(obj)), file = file_conn) ; cat("\n", file = file_conn)}
      close(file_conn)
      }


##################################################################################################################################

#IV: Read and display the unfiltered data

##################################################################################################################################

    #(7): Load shapefile (i.e. a single polygon or a multipolygon)

       #Read and plot the shapefile
       root_fldr <- here()
       aoi_Polygons <- st_read(paste0(root_fldr, "/Input/Shapefiles/", shapefile), quiet=T)
       p1 <- ggplot() +
         geom_sf(data = aoi_Polygons, fill=sf.colors(nrow(aoi_Polygons)), col = "black")+
         theme_tom()+
         geom_sf_label(data = aoi_Polygons, aes(label=LocationID), colour="black")

       #Save the plot only once, because the shapefile is identical for all years
       if(year_ID == year_ID_vector[1]){
         tryCatch({print(p1)},
            #Catch the 'st_point_on_surface' warning
             warning = function(w){
               if(grepl("st_point_on_surface", w$message)){
                 suppressWarnings(ggsave(plot=p1, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_",paste0(area_name, "_MODIS_", MODIS_cloud_masking_algorithm), "_Polygons.pdf"), width=10, height=8))
                 suppressWarnings(print(p1))
                 }
               if(!grepl("st_point_on_surface", w$message)){
                 ggsave(plot=p1, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_",paste0(area_name, "_MODIS_", MODIS_cloud_masking_algorithm), "_Polygons.pdf"), width=10, height=8)
                 print(p1)
                 }
               })}
       
       #Convert the shapefile to an earthengine feature collection:
       aoi_Polygons <- st_transform(aoi_Polygons, crs="EPSG:4326")
       aoi_Polygons <- sf_as_ee(aoi_Polygons)
       aoi_Polygons <- ee$FeatureCollection(aoi_Polygons)
       #Note that the feature$property "LocationID" has a unique number for each polygon

       #Collapse all polygons (if any) into a single polygon and calculate the (convex hull) outline of the study area
       aoi_Polygons_merged <- ee$FeatureCollection(aoi_Polygons$geometry()$dissolve())
       aoi_Shapefile <- ee$FeatureCollection(aoi_Polygons_merged$geometry()$convexHull())

         # #Inspect merged polygon (for debugging)
         # Map$centerObject(aoi_Polygons)
         # Map$addLayer(aoi_Polygons,  list(color="blue"), 'Original Polygons')+
         #   Map$addLayer(aoi_Polygons_merged,  list(color="red"), 'Merged Polygon')+
         #   Map$addLayer(aoi_Shapefile,  list(color="green"), 'Shapefile outline')

       #Calculate central point
       coordinates_point <- aoi_Shapefile$geometry()$centroid()

       #Plot the shapefiles and the central point of the study area
       Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 6)
       Map$addLayer(aoi_Polygons, list(color="grey"), name='Polygons')+
         Map$addLayer(aoi_Shapefile, list(color="black"), name='Study area')+
         Map$addLayer(coordinates_point, list(color="blue"), name='Central point')

       #Calculate the size of the study area in km2:
       img <- ee$Image$pixelArea()$divide(1000000)
       area2 <- img$reduceRegion(
         reducer= ee$Reducer$sum(),
         geometry= aoi_Shapefile,
         crs=crs,
         scale= resolution,
         maxPixels= 1E13)
       paste0('Size of study area calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')

    #(8) Extract MODIS satellite Surface Reflectance images:

      #Extract MODIS satellite data
       MODIS_col <- ee$ImageCollection(MODIS_dataset)
       MODIS_col <- MODIS_col$
         filterBounds(aoi_Shapefile)$
         filterDate(start_date, end_date)

       #Note that MODIS images are global composites which will always intersect with aoi_Shapefile,
       #and this full image is thus retained (i.e. filterBounds does not do anything).

      #Clip all images to the area depicted by 'aoi_Shapefile'
       MODIS_col <- MODIS_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})

       #Note that we clip the imagecollection to aoi_Shapefile and not to aoi_Polygons, to prevent issues
       #when the shapefile consists of multiple non adjacent polygons. The final calculations below are
       #however conducted for the area depicted by 'aoi_Polygons'.

    #(9) Add new bands to the clipped image collection

      #Map Normalized difference functions and FSC functions over the image collection
      MODIS_col <- MODIS_col$
        map(getNDSI)$
        map(getNDVI)$
        map(getNDMI)$
        map(getNDWI)$
        map(get_FSC)$
        map(add_Date)

      # # #Plot all layers (for debugging)
      #  image <- MODIS_col$filterDate(paste0(year_ID, "-06-18"), end_date)$first()
      #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
      #  Map$addLayer(image,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
      #  Map$addLayer(image,list(bands=c("NDSI"), min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
      #  #Map$addLayer(image,list(bands=c("NDVI"), min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
      #  #Map$addLayer(image,list(bands=c("NDWI"), min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')+
      #  #Map$addLayer(image,list(bands=c("FSC_Gascoin2020"), min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'FSC_Gascoin2020')+
      #  #Map$addLayer(image,list(bands=c("FSC_Aalstad2020"), min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'FSC_Aalstad2020')+
      #  Map$addLayer(aoi_Polygons, list(color="grey"), name='Polygons')+
      #  Map$addLayer(aoi_Shapefile, list(color="red"), name='Study area')+
      #  Map$addLayer(coordinates_point, list(color="black"), name='Central point')
      
      # #Extract bandvalues at point location (for debugging)
      #  ee_extract(x=image$select("NDSI", "FSC_Gascoin2020", "FSC_Aalstad2020"), y=coordinates_point, fun=ee$Reducer$first(), scale=resolution, sf=TRUE)

    #(10) Create a timeseries GIF of unmasked RGB and NDSI images for aoi_Shapefile
      if(gif_output==TRUE){
        
        #(A): GIF of RGB-images
        f_img_col_to_gif(img_col=MODIS_col,  #image collection
                         RGB_bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), #names of RGB-bands in image collection
                         shapefile=aoi_Shapefile, #shapefile of area of interest
                         centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
                         gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
                         gif_fps=5, #Frames per second of GIF
                         gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
                         gif_max=12000, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
                         gif_gamma=c(1.9, 1.7, 1.7), #Gamma correction factors (one for each band)
                         gif_crs='EPSG:3857', #CRS project of the output
                         gif_text_position="northwest", #Location of text (datetime string)
                         gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
                         gif_text_size=14, #Size of text
                         gif_text_col="#FFFFFF",
                         output_fldr=paste0(root_fldr, "/Output/MODIS/03_Polygons_Snowmelt/"),
                         file_name="GIF_RGB")
        
        #(b): GIF of NDSI-images
        
          #Transform NDSI-band to RGB bands
          MODIS_col_ndsi <- MODIS_col$map(ee_utils_pyfunc(function(image){
            f_band_to_RGB(img=image, band='NDSI', min_value=-1, max_value=1.5, palette=c('#000000', '#0dffff', '#0524ff', '#ffffff'))}))
          
          #Create gif of the transformed NDSI-band
          f_img_col_to_gif(img_col=MODIS_col_ndsi,  #image collection
                           RGB_bands=c("vis-red", "vis-green", "vis-blue"), #names of RGB-bands in image collection
                           shapefile=aoi_Shapefile, #shapefile of area of interest
                           centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
                           gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
                           gif_fps=5, #Frames per second of GIF
                           gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
                           gif_max=255, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
                           gif_gamma=c(1, 1, 1), #Gamma correction factors (one for each band)
                           gif_crs='EPSG:3857', #CRS project of the output
                           gif_text_position="northwest", #Location of text (datetime string)
                           gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
                           gif_text_size=14, #Size of text
                           gif_text_col="#FFFFFF",
                           output_fldr=paste0(root_fldr, "/Output/MODIS/03_Polygons_Snowmelt/"),
                           file_name="GIF_NDSI")
        
      }

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
          map(AddNULLCloudFraction)

       # #Check if Cloud information has been added to the properties of each image (for debugging)
       # MODIS_col$first()$propertyNames()$getInfo()

       # #Extract cloud fraction of all images in image collection for the area of interest (for debugging)
       # Combined_clouds_Fraction <- unlist(MODIS_col$aggregate_array('Combined_clouds_Fraction')$getInfo())
       # PGE11_clouds_Fraction <- unlist(MODIS_col$aggregate_array('PGE11_clouds_Fraction')$getInfo())
       # MOD35_clouds_Fraction <- unlist(MODIS_col$aggregate_array('MOD35_clouds_Fraction')$getInfo())
       # doy <- unlist(MODIS_col$aggregate_array('doy')$getInfo())
       #
       # #Replace -9999 values by NA (for debugging)
       # Combined_clouds_Fraction[Combined_clouds_Fraction < -9000] <- NA
       # PGE11_clouds_Fraction[PGE11_clouds_Fraction < -9000] <- NA
       # MOD35_clouds_Fraction[MOD35_clouds_Fraction < -9000] <- NA
       # doy[doy < -9000] <- NA
       #
       # #Combine all cloud fraction measures in a single dataframe (for debugging)
       # clouds_MODIS <- data.frame(doy=doy,
       #                            Combined_clouds_Fraction=Combined_clouds_Fraction,
       #                            PGE11_clouds_Fraction=PGE11_clouds_Fraction,
       #                            MOD35_clouds_Fraction=MOD35_clouds_Fraction)
       #
       # #Assess relationship between both cloud algorithms (for debugging)
       # ggplot()+
       #   geom_line(aes(x=doy, y=MOD35_clouds_Fraction), col="red", lwd=1)+
       #   geom_line(aes(x=doy, y=PGE11_clouds_Fraction), col="black", lwd=1)+
       #   #geom_line(aes(x=doy, y=Combined_clouds_Fraction), col="blue", lwd=1)+
       #   theme_tom()
       #
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
       # image <- MODIS_col$
       #    filterDate(paste0(year_ID, "-06-17"), end_date)$
       #    first()$
       #    clipToCollection(aoi_Shapefile)$
       #    select('MOD35_clouds', 'PGE11_clouds', 'Combined_clouds', "sur_refl_b01", "sur_refl_b04", "sur_refl_b03")
       # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
       # Map$addLayer(image, list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
       # Map$addLayer(image, list(bands='MOD35_clouds', min=0, max=1, opacity=1), 'MOD35_clouds')+
       # Map$addLayer(image, list(bands='PGE11_clouds', min=0, max=1, opacity=1), 'PGE11_clouds')+
       # Map$addLayer(image, list(bands='Combined_clouds', min=0, max=1, opacity=1), 'MOD35 and PGE11 combined')

       #These images also indicate that the PGE11 algorithm does a more conservative job that better matches clouds in the RGB images than the MOD35 algorithm.

     #(12): Exclude all images from the image collection that have a CloudFraction value >= max_cloud_fraction, and mask all cloud pixels within the remaining images

       #Cloud function III: Define a cloudmask function for the selected cloud masking algorithm
        cloud_algorithm <- paste0(MODIS_cloud_masking_algorithm, "_clouds")

       #Apply a manually constructed cloudmask to each image in the collection:
        clouds_fraction <- paste0(MODIS_cloud_masking_algorithm, "_clouds_Fraction")
        MODIS_clouds_filtered <- MODIS_col$
          #Filter all images with 100% cloudcover (PGE11_clouds_Fraction=1.0)
          filter(ee$Filter$lt(clouds_fraction, max_cloud_fraction))$
          #Apply cloudmask for individual pixels
          map(AddCloudMask)
      }
  if(mask_clouds==FALSE){

    #print message
    print("Cloud masking = FALSE")

    #Set MODIS_clouds_filtered equal to MODIS_col
    MODIS_clouds_filtered <- MODIS_col

  }

  #Store default MODIS image projection
  modisProjection <- MODIS_clouds_filtered$first()$select("NDSI")$projection()
  #modisProjection$getInfo()

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

     }
  if(mask_water==FALSE){

     #(A): print message
     print("Water masking = FALSE")

   }

##################################################################################################################################

#VII: Count the total number of unmasked pixels per doy within aoi_Shapefile

##################################################################################################################################

  #(14): Count the total number of unmasked pixels and the total number of pixels per doy within aoi_Shapefile
  if(pixel_counts==TRUE){

    #(A): Add pixel counts within aoi_Shapefile to each separate image by mapping the pixel count functions over the image collection
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
         description = paste0(current_timestamp0, "_", data_ID, "_Res", resolution, "_Polygons_Data_Pixel_Counts_polygon"),
         fileFormat = "CSV",
         selectors = c('doy', 'unmasked', 'total')
         )

      #Monitor the task
       task_vector0$start()
       print("Count the number of unmasked pixels within the shapefile per doy:")
       ee_monitoring(task_vector0, quiet=T, max_attempts=1000000)

      #Export results to local folder
       exported_stats <- ee_drive_to_local(task = task_vector0, dsn=paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_Pixel_Counts_polygon"))
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
     df_pixelcount$location <- area_name
     df_pixelcount$year <- year_ID
     write.csv(df_pixelcount, file=paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_Pixel_Counts_polygon.csv"), quote=FALSE, row.names=FALSE)

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
      ggtitle(paste0("Year: ", year_ID))+
      theme_classic()

    #(F): Save barplot
     pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Plot_Pixel_Counts_polygon.pdf"), width=12, height=8)
     print(p_pixelcounts)
     dev.off()

  }

##################################################################################################################################

#VIII: Create timelapse videos of the RGB bands, and the NDSI-, NDVI-, NDMI- and NDMI-bands for aoi_Polygons (Cloud and Water masked)

##################################################################################################################################

  #(15): Create timeseries GIF animations
  if(gif_output==TRUE){
    
    #(a): GIF of RGB-images
    f_img_col_to_gif(img_col=MODIS_clouds_filtered,  #image collection
                     RGB_bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), #names of RGB-bands in image collection
                     shapefile=aoi_Shapefile, #shapefile of area of interest
                     centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
                     gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
                     gif_fps=5, #Frames per second of GIF
                     gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
                     gif_max=12000, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
                     gif_gamma=c(1.9, 1.7, 1.7), #Gamma correction factors (one for each band)
                     gif_crs='EPSG:3857', #CRS project of the output
                     gif_text_position="northwest", #Location of text (datetime string)
                     gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
                     gif_text_size=14, #Size of text
                     gif_text_col="#FFFFFF",
                     output_fldr=paste0(root_fldr, "/Output/MODIS/03_Polygons_Snowmelt/"),
                     file_name="GIF_RGB_masked")
    
    #(b): GIF of NDSI-images
    
      #Transform NDSI-band to RGB bands
      MODIS_clouds_filtered_ndsi <- MODIS_clouds_filtered$map(ee_utils_pyfunc(function(image){
        f_band_to_RGB(img=image, band='NDSI', min_value=-1, max_value=1.5, palette=c('#000000', '#0dffff', '#0524ff', '#ffffff'))}))
      
      #Create gif of the transformed NDSI-band
      f_img_col_to_gif(img_col=MODIS_clouds_filtered_ndsi,  #image collection
                       RGB_bands=c("vis-red", "vis-green", "vis-blue"), #names of RGB-bands in image collection
                       shapefile=aoi_Shapefile, #shapefile of area of interest
                       centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
                       gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
                       gif_fps=5, #Frames per second of GIF
                       gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
                       gif_max=255, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
                       gif_gamma=c(1, 1, 1), #Gamma correction factors (one for each band)
                       gif_crs='EPSG:3857', #CRS project of the output
                       gif_text_position="northwest", #Location of text (datetime string)
                       gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
                       gif_text_size=14, #Size of text
                       gif_text_col="#FFFFFF",
                       output_fldr=paste0(root_fldr, "/Output/MODIS/03_Polygons_Snowmelt/"),
                       file_name="GIF_NDSI_masked")
    
    # #(c): GIF of NDVI-images
    #    
    #    #Transform NDVI-band to RGB bands
    #    MODIS_clouds_filtered_ndvi <- MODIS_clouds_filtered$map(ee_utils_pyfunc(function(image){
    #      f_band_to_RGB(img=image, band='NDVI', min_value=-0.25, max_value=1, palette=c("#cccccc", "#f46d43", "#fdae61", "#fee08b", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850"))}))
    #    
    #    #Create gif of the transformed NDVI-band
    #    f_img_col_to_gif(img_col=MODIS_clouds_filtered_ndvi,  #image collection
    #                     RGB_bands=c("vis-red", "vis-green", "vis-blue"), #names of RGB-bands in image collection
    #                     shapefile=aoi_Shapefile, #shapefile of area of interest
    #                     centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
    #                     gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
    #                     gif_fps=5, #Frames per second of GIF
    #                     gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
    #                     gif_max=255, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
    #                     gif_gamma=c(1, 1, 1), #Gamma correction factors (one for each band)
    #                     gif_crs='EPSG:3857', #CRS project of the output
    #                     gif_text_position="northwest", #Location of text (datetime string)
    #                     gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
    #                     gif_text_size=14, #Size of text
    #                     gif_text_col="#FFFFFF",
    #                     output_fldr=paste0(root_fldr, "/Output/MODIS/03_Polygons_Snowmelt/"),
    #                     file_name="GIF_NDVI_masked")
    # 
    # #(d): GIF of NDMI-images
    #    
    #    #Transform NDMI-band to RGB bands
    #    MODIS_clouds_filtered_ndmi <- MODIS_clouds_filtered$map(ee_utils_pyfunc(function(image){
    #      f_band_to_RGB(img=image, band='NDMI', min_value=-0.75, max_value=1, palette=c("#d73027", "#f46d43", "#fdae61", "#fee08b", "#d9ef8b", "#a6d96a", "#66bd63", "#1a9850", "#6ad99e", "#387ad9", "#003dd6"))}))
    #    
    #    #Create gif of the transformed NDMI-band
    #    f_img_col_to_gif(img_col=MODIS_clouds_filtered_ndmi,  #image collection
    #                     RGB_bands=c("vis-red", "vis-green", "vis-blue"), #names of RGB-bands in image collection
    #                     shapefile=aoi_Shapefile, #shapefile of area of interest
    #                     centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
    #                     gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
    #                     gif_fps=5, #Frames per second of GIF
    #                     gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
    #                     gif_max=255, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
    #                     gif_gamma=c(1, 1, 1), #Gamma correction factors (one for each band)
    #                     gif_crs='EPSG:3857', #CRS project of the output
    #                     gif_text_position="northwest", #Location of text (datetime string)
    #                     gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
    #                     gif_text_size=14, #Size of text
    #                     gif_text_col="#FFFFFF",
    #                     output_fldr=paste0(root_fldr, "/Output/MODIS/03_Polygons_Snowmelt/"),
    #                     file_name="GIF_NDMI_masked")
    # 
    # #(e): GIF of NDWI-images
    #    
    #    #Transform NDWI-band to RGB bands
    #    MODIS_clouds_filtered_ndwi <- MODIS_clouds_filtered$map(ee_utils_pyfunc(function(image){
    #      f_band_to_RGB(img=image, band='NDWI', min_value=-0.5, max_value=1, palette=c('#000000', '#0dffff', '#0524ff', '#ffffff'))}))
    #    
    #    #Create gif of the transformed NDWI-band
    #    f_img_col_to_gif(img_col=MODIS_clouds_filtered_ndwi,  #image collection
    #                     RGB_bands=c("vis-red", "vis-green", "vis-blue"), #names of RGB-bands in image collection
    #                     shapefile=aoi_Shapefile, #shapefile of area of interest
    #                     centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
    #                     gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
    #                     gif_fps=5, #Frames per second of GIF
    #                     gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
    #                     gif_max=255, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
    #                     gif_gamma=c(1, 1, 1), #Gamma correction factors (one for each band)
    #                     gif_crs='EPSG:3857', #CRS project of the output
    #                     gif_text_position="northwest", #Location of text (datetime string)
    #                     gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
    #                     gif_text_size=14, #Size of text
    #                     gif_text_col="#FFFFFF",
    #                     output_fldr=paste0(root_fldr, "/Output/MODIS/03_Polygons_Snowmelt/"),
    #                     file_name="GIF_NDWI_masked")
    
  }

##################################################################################################################################

#IX: Extract Average FSC, NDSI, NDVI, and NDMI per Polygon

##################################################################################################################################

 #Extract MODIS band values (FSC, NDSI, NDVI, and NDMI) for each polygon over time
  if("avg_NDSI" %in% method){

    #Print message
    cat("\n")
    print("--------------------------------------------------------------------------------------------------------------------------")
    print(paste0("METHOD: 'avg_NDSI' - CALCULATING THE AVERAGE FSC, NDSI, NDVI AND NDMI per Polygon"))
    print("--------------------------------------------------------------------------------------------------------------------------")

    #Create an iteration function that we will use to iterate through all images of the image collection. For each image,
    #the value of certain bands is extracted for each feature (i.e. Polygon) within the feature collection aoi_Polygons.
    #The resulting band values are added as a property to each feature (i.e. Polygon). This results in an updated feature
    #collection specific for the current image. This feature collection is then appended to a list of feature collections
    #from previous image iterations. After iterating through n-images, the result is thus a feature collection list of
    #length n * the number of features within the feature collection.

    #Thus, at each iteration we extract band values of interest for all Polygons within the current image, store this as a
    #feature collection and add this to an expanding list of feature collections from previous iterations.

    #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
     FC_initial <- ee$FeatureCollection(ee$List(list()))

    #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
    #(in this case the current iteration image) and the second element takes the output value from the iteration that preceded it. The latter
    #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with
    #has been defined.
     Extract_BandValuesAtPolygons = Extract_BandValuesAtPolygons #sourced

    #Iterate this function over all images in the ImageCollection
     FC_merged <- ee$FeatureCollection(MODIS_clouds_filtered$select("FSC_Gascoin2020", "FSC_Aalstad2020", "NDSI", "NDVI", "NDMI")$iterate(Extract_BandValuesAtPolygons, FC_initial))
     #FC_merged$getInfo() #for debugging
     #FC_merged$first()$getInfo() #for debugging

    #Transform feature collection with MODIS Band values for each polygon to a dataframe

      #create a current timestamp to prevent identical names on Google Drive
       current_timestamp1 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))

      #Setup task
       task_vector1 <- ee_table_to_drive(
         collection= FC_merged,
         description = paste0(current_timestamp1, "_", data_ID, "_Res", resolution, "_Data_MeanBandValues_Polygons"),
         folder="RGEE_tmp",
         fileFormat="CSV",
         selectors=c('FSC_Gascoin2020', 'FSC_Aalstad2020','NDSI', 'NDVI', 'NDMI', 'Date', 'LocationID')
         )

      #Run and monitor task
       task_vector1$start()
       ee_monitoring(task_vector1, task_time=60, max_attempts=1000000)

      #Import results
       exported_stats <- ee_drive_to_local(task=task_vector1, dsn=paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Data_MeanBandValues_Polygons"))
       df_Polygons_BandValues <- read.csv(exported_stats)
       unlink(exported_stats)

      #Restructure dataframe
       df_Polygons_BandValues$doy <- as.numeric(strftime(df_Polygons_BandValues$Date, format = "%j"))
       colnames(df_Polygons_BandValues) <- c('FSC_Gascoin2020', 'FSC_Aalstad2020', 'NDSI', 'NDVI', 'NDMI', 'Date', 'Polygon', 'doy')
       df_Polygons_BandValues <- df_Polygons_BandValues[ ,c('FSC_Gascoin2020', 'FSC_Aalstad2020', 'NDSI', 'NDVI', 'NDMI', 'Date', 'doy', 'Polygon')]

      #Change -9999 to NA
       df_Polygons_BandValues$FSC_Gascoin2020[df_Polygons_BandValues$FSC_Gascoin2020 < -9000] <- NA
       df_Polygons_BandValues$FSC_Aalstad2020[df_Polygons_BandValues$FSC_Aalstad2020 < -9000] <- NA
       df_Polygons_BandValues$NDSI[df_Polygons_BandValues$NDSI < -9000] <- NA
       df_Polygons_BandValues$NDVI[df_Polygons_BandValues$NDVI < -9000] <- NA
       df_Polygons_BandValues$NDMI[df_Polygons_BandValues$NDMI < -9000] <- NA
       df_Polygons_BandValues$doy[df_Polygons_BandValues$doy < -9000] <- NA
       df_Polygons_BandValues$Polygon[df_Polygons_BandValues$Polygon < -9000] <- NA

      #Sort dataframe by Polygon and doy
       index <- with(df_Polygons_BandValues, order(Polygon, doy))
       df_Polygons_BandValues <- df_Polygons_BandValues[index,]

      #Add LocationID and Year
       df_Polygons_BandValues$location <- area_name
       df_Polygons_BandValues$year <- year_ID

      #Store dataframe with average bandValues and Snowfraction for all polygons
       #write.csv(df_Polygons_BandValues, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_MeanBandValues.csv"), row.names = FALSE)
       #df_Polygons_BandValues <- read.csv(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_MeanBandValues.csv"), header=T)

  }


##################################################################################################################################

#X: Extract the fraction of pixels with NDSI > NDSI_threshold over time per Polygon

##################################################################################################################################

 #Extract the fraction of pixels with NDSI > NDSI_threshold over time per Polygon for all values of NDSI_threshold_vector
  if("snowfraction" %in% method){

    #Print message
    cat("\n")
    print("--------------------------------------------------------------------------------------------------------------------------")
    print(paste0("METHOD: 'snowfraction' - CALCULATING THE FRACTION OF PIXELS WITH NDSI > NDSI_threshold PER POLYGON"))
    print("--------------------------------------------------------------------------------------------------------------------------")

    #Specify the NDSI threshold(s) at which an area is perceived as snow-free
    NDSI_threshold_vector = NDSI_threshold_vector

    #Create an empty dataframe for storing output per Polygon
    df_Polygons_SnowFraction <- data.frame(SnowFraction=numeric(),
                                           Date=character(),
                                           doy=numeric(),
                                           Polygon=character(),
                                           NDSI_threshold=character())

    #Run the analysis for each level of NDSI_threshold_vector
    for(NDSI_threshold in NDSI_threshold_vector){

       #Print message
       print(paste0("  -START ANALYSIS FOR NDSI_threshold = ", NDSI_threshold))

       #Store NDSI_threshold as a character (used for naming of outputs)
       NDSI_threshold_char <- gsub("\\.", "_", as.character(NDSI_threshold))

       #Add binary SNOW-band to the image collection that identifies for each pixel whether it is snowcovered (NDSI > NDSI_threshold)

           #Map computeSnow function over the image collection
           modis_snow_masked <- MODIS_clouds_filtered$
               map(computeSnow)

           # #Visualize the binary layer (for debugging).
           #  image <- modis_snow_masked$filterDate(paste0(year_ID, "-06-15"), end_date)$first()
           #  Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
           #  Map$addLayer(image, list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
           #  Map$addLayer(image$select('SNOW'), list(min=0, max = 1, palette = c('ffffff', 'orange')))

       #Extract fraction of snowcover within each polygon for each MODIS image

           #Create a function that iterates over the image collection, extracts the fraction of snowcover within all Polygons
           #of each image image, stores this as a feature collection and adds this to an expanding list of feature collections
           #from previous iterations.

           #Create an empty FeatureCollection list. This list is used as input for the first iteration of the iteration function below.
           FC_initial_SnowFraction <- ee$FeatureCollection(ee$List(list()))

           #Specify the iteration function. This function takes two arguments. The first argument is the current element of the image collection
           #(in this case the current iteration image) and the second element takes the output value from the iteration that preceeded it. The latter
           #is not possible for the first iteration, that's why an initial object (empty feature collection) to start the iteration with
           #has been defined.
           Extract_SnowFractionAtPolygons = function(img, FC_initial_SnowFraction){

                #Take the mean of the binary SNOW band to get the fraction of snow covered pixels for each image. ReduceRegions does not include masked pixels (i.e. pixels
                #defined as cloud or water) when calculating these means.
                FC_image <- img$reduceRegions(collection=aoi_Polygons,
                                              reducer=ee$Reducer$mean(),
                                              scale=resolution,
                                              crs=crs,
                                              crsTransform=NULL)

                #Add datetime of current image as a property to each feature within the feature collection FC_image
                date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
                FC_image <- FC_image$map(function(feature){return(feature$set("Date", date))})

                #Make sure there is a value at each feature within the feature collection:
                FC_image <- FC_image$map(function(feature){
                  snow <- ee$List(list(feature$get('SNOW'), -9999))$reduce(ee$Reducer$firstNonNull())
                  date <- ee$List(list(feature$get('Date'), -9999))$reduce(ee$Reducer$firstNonNull())
                  return(feature$
                           set("SNOW", snow)$
                           set("Date", date)
                           )})

                #Merge the feature collection of the current image (FC_image) onto the feature collection FC_initial_SnowFraction.
                return (ee$FeatureCollection(FC_initial_SnowFraction)$merge(FC_image))

              }

           #Iterate this function over all images in the ImageCollection
           FC_merged_SnowFraction <- ee$FeatureCollection(modis_snow_masked$select("SNOW", "NDSI")$iterate(Extract_SnowFractionAtPolygons, FC_initial_SnowFraction))
           #Note that for some reason, including a second band next to SNOW is required for the code to function.
           #FC_merged_SnowFraction$getInfo() #for debugging
           #FC_merged_SnowFraction$first()$getInfo() #for debugging

       #Transform feature collection with MODIS Band values for each polygon to a dataframe

           #We export the data instead of using aggregate_array() as the latter might fail due to computation timeouts.

           #create a current timestamp to prevent identical names on Google Drive
           current_timestamp2 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))

           #Setup task
           task_vector2 <- ee_table_to_drive(
              collection= FC_merged_SnowFraction,
              description = paste0(current_timestamp2, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_SnowFraction_Polygons"),
              folder="RGEE_tmp",
              fileFormat="CSV",
              selectors=c('SNOW', 'Date', 'LocationID')
              )

           #Run and monitor task
           task_vector2$start()
           ee_monitoring(task_vector2, task_time=60, max_attempts=1000000)

           #Import results
           exported_stats <- ee_drive_to_local(task=task_vector2, dsn=paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI", NDSI_threshold_char, "_Data_SnowFraction_Polygons"))
           df_Polygons_SnowFraction_new <- read.csv(exported_stats)
           unlink(exported_stats)

           #Restructure dataframe
           df_Polygons_SnowFraction_new$doy <- as.numeric(strftime(df_Polygons_SnowFraction_new$Date, format = "%j"))
           colnames(df_Polygons_SnowFraction_new) <- c('SnowFraction', 'Date', 'Polygon', 'doy')
           df_Polygons_SnowFraction_new <- df_Polygons_SnowFraction_new[ ,c('SnowFraction', 'Date', 'doy', 'Polygon')]

           #Change -9999 to NA
           df_Polygons_SnowFraction_new$SnowFraction[df_Polygons_SnowFraction_new$SnowFraction < -9000] <- NA
           df_Polygons_SnowFraction_new$doy[df_Polygons_SnowFraction_new$doy < -9000] <- NA
           df_Polygons_SnowFraction_new$Polygon[df_Polygons_SnowFraction_new$Polygon < -9000] <- NA

           #Sort dataframe by Polygon and doy
           index <- with(df_Polygons_SnowFraction_new, order(Polygon, doy))
           df_Polygons_SnowFraction_new <- df_Polygons_SnowFraction_new[index,]

           #Add NDSI_threshold as a new column
           df_Polygons_SnowFraction_new$NDSI_threshold <- as.factor(NDSI_threshold)

           #Add dataframe for current NDSI_threshold to dataframe from previous iterations:
           df_Polygons_SnowFraction <- rbind(df_Polygons_SnowFraction, df_Polygons_SnowFraction_new)

           # #For debugging
           #  ggplot() + geom_point(data=df_Polygons_SnowFraction, aes(x=doy, y=SnowFraction, col=NDSI_threshold)) + theme_classic()

           # #Save dataframe for current NDSI_threshold
           #  write.csv(df_Polygons_SnowFraction_new, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Resolution", resolution, "_NDSI", NDSI_threshold_char, "_Data_SnowFraction_Polygons.csv"), row.names = FALSE)

         }

    #Add LocationID and year
    df_Polygons_SnowFraction$location <- area_name
    df_Polygons_SnowFraction$year <- year_ID

    #Store dataframe with Snowfraction data for all polygons for all levels of NDSI_threshold
    #write.csv(df_Polygons_SnowFraction, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_SnowFraction.csv"), row.names = FALSE)
    #df_Polygons_SnowFraction <- read.csv(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_SnowFraction.csv"), header=T)

  }


##################################################################################################################################

#XI: Fit GAMS through the Average NDSI, NDVI, NDMI and SnowFraction data per Polygon

##################################################################################################################################

   #(23): Fit GAMS

      #We fit a Generalized Additive Model (GAM) through the data within each polygon. We do this using a sequential
      #process. We first fit a GAM through the data, calculating model predictions and residuals. We then exclude all
      #rows from the dataframe where the residual >= (0.4 * the range of the data). We then re-fit a GAM to this reduced
      #dataset, make predictions and calculate residuals. We then exclude all rows from the reduced dataframe where the
      #residual >= >= (0.2 * the range of the data). This gives us a final dataframe through which we fit a third
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

        #Specify at which NDSI value, the Polygon is considered snow-free
         NDSI_threshold_vector = NDSI_threshold_vector

########################################################################################################################################################################################

  #(I): SNOWFRACTION - Fit a Generalized Additive Model (GAM) through the snow fraction data within each polygon
    if("snowfraction" %in% method){

      #Print message
      cat("\n")
      print("--------------------------------------------------------------------------------------------------------------------------")
      print(paste0("FIT A GAM THROUGH THE FRACTION OF OF PIXELS WITH NDSI > NDSI_threshold PER POLYGON"))
      print("--------------------------------------------------------------------------------------------------------------------------")

      #(1) Create an empty dataframe
        df_Polygons_SnowFraction_GAM <- data.frame(NDSI_threshold=character(),
                                                   SnowFraction=numeric(),
                                                   Date=factor(),
                                                   doy=numeric(),
                                                   Polygon=factor(),
                                                   outliers=logical())

        df_Polygons_SnowFraction_GAM_predictions <- data.frame(Polygon=character(),
                                                               NDSI_threshold=character(),
                                                               doy=numeric(),
                                                               SnowFraction_gam_predict=numeric(),
                                                               stringsAsFactors=FALSE)

      #(2) Loop through all Polygons and fit a separate GAM (with sequential outlier removal) to the Polygon-specific SnowFraction data

          #Loop through all Polygons
          for(i in unique(df_Polygons_SnowFraction$Polygon)){

            #For debugging
            #i=unique(df_Polygons_SnowFraction$Polygon)[2]

            #Loop through all NDSI_thresholds
            for(j in unique(df_Polygons_SnowFraction$NDSI_threshold)){

              #For debugging
              #j=unique(df_Polygons_SnowFraction$NDSI_threshold)[1]

              #Select Polygon-specific subset of data:
              df_Polygon_SnowFraction_GAM_new <- df_Polygons_SnowFraction[df_Polygons_SnowFraction$Polygon==i &
                                                                          df_Polygons_SnowFraction$NDSI_threshold==j &
                                                                          !is.na(df_Polygons_SnowFraction$SnowFraction),
                                                                          c("NDSI_threshold", "SnowFraction", "Date", "doy", "Polygon")]

              #If there are at least two datapoints, then continue to fit the GAM
              if(nrow(df_Polygon_SnowFraction_GAM_new) > 1){

                #Fit a GAM through the Polygon-specific SnowFraction ~ doy data and employ sequential outlier removal
                df_Polygon_SnowFraction_GAM_new <- f_gam_SeqRemOutliers(data=df_Polygon_SnowFraction_GAM_new, y="SnowFraction", x="doy", outlier_removal=outlier_removal,
                                                                        outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                                        default_k=gam_k_outlier)

                #Sort df_Polygon_SnowFraction_GAM_new by doy:
                df_Polygon_SnowFraction_GAM_new <- df_Polygon_SnowFraction_GAM_new[order(df_Polygon_SnowFraction_GAM_new$doy),]

                #Bind the Polygon-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
                 df_Polygons_SnowFraction_GAM <- rbind(df_Polygons_SnowFraction_GAM, df_Polygon_SnowFraction_GAM_new)

                #Create more detailed predictions (not only at the doy present in the dataframe) to plot more smooth curves

                  #Refit GAM through data
                   index <- which(df_Polygon_SnowFraction_GAM_new$outliers==FALSE)
                   mod_gam <- with(df_Polygon_SnowFraction_GAM_new[index,], mgcv::gam(SnowFraction ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

                  #Use GAM to make predictions on a more detailed (1-day) day of year interval
                   df_Polygon_SnowFraction_GAM_predictions_new <- data.frame(Polygon=i, NDSI_threshold=j, doy=seq(min(df_Polygon_SnowFraction_GAM_new$doy), max(df_Polygon_SnowFraction_GAM_new$doy), 0.01))
                   df_Polygon_SnowFraction_GAM_predictions_new$SnowFraction_gam_predict <- stats::predict(mod_gam, newdata=df_Polygon_SnowFraction_GAM_predictions_new, type="response")
                   df_Polygon_SnowFraction_GAM_predictions_new <- df_Polygon_SnowFraction_GAM_predictions_new[order(df_Polygon_SnowFraction_GAM_predictions_new$doy),]

                  #Add predictions to df_Polygons_SnowFraction_GAM_predictions dataframe:
                   df_Polygons_SnowFraction_GAM_predictions <- rbind(df_Polygons_SnowFraction_GAM_predictions, df_Polygon_SnowFraction_GAM_predictions_new)

                  }

              #If there are less than two datapoints, then append empty dataframes
              if(nrow(df_Polygon_SnowFraction_GAM_new) < 2){

                #Bind the empty dataframe to the dataframe containing all dataframes from previous iterations
                if(nrow(df_Polygon_SnowFraction_GAM_new) < 1){df_Polygon_SnowFraction_GAM_new <- data.frame(NDSI_threshold=j, SnowFraction=NA, Date=NA, doy=NA, Polygon=i, outliers=NA)}
                if(nrow(df_Polygon_SnowFraction_GAM_new) > 0){df_Polygon_SnowFraction_GAM_new$outliers <- NA}
                df_Polygons_SnowFraction_GAM <- rbind(df_Polygons_SnowFraction_GAM, df_Polygon_SnowFraction_GAM_new)

                #Add empty GAM predictions to df_Polygons_SnowFraction_GAM_predictions dataframe:
                df_Polygon_SnowFraction_GAM_predictions_new <- df_Polygon_SnowFraction_GAM_new[,c("Polygon", "NDSI_threshold", "doy")]
                df_Polygon_SnowFraction_GAM_predictions_new$SnowFraction_gam_predict <- NA
                df_Polygons_SnowFraction_GAM_predictions <- rbind(df_Polygons_SnowFraction_GAM_predictions, df_Polygon_SnowFraction_GAM_predictions_new)

                }

              }

            }

          #Change Polygon column to a factor
          df_Polygons_SnowFraction_GAM$Polygon <- as.factor(as.character(df_Polygons_SnowFraction_GAM$Polygon))
          df_Polygons_SnowFraction_GAM_predictions$Polygon <- as.factor(as.character(df_Polygons_SnowFraction_GAM_predictions$Polygon))

          #Add LocationID and year
          df_Polygons_SnowFraction_GAM_predictions$location <- area_name
          df_Polygons_SnowFraction_GAM_predictions$year <- year_ID

          #Save dataframe with GAM fits for SnowFraction
          #write.csv(df_Polygons_SnowFraction_GAM, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_SnowFraction_GAM.csv"), row.names = FALSE)
          #write.csv(df_Polygons_SnowFraction_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_GAM_Predictions_SnowFraction.csv"), row.names = FALSE)

      #(3) Plot the raw Snowfraction datapoints and gam predictions for each Polygon:

          #Plot Snowfraction and model predictions for all Polygons in a single plot
          p_Polygon_SnowFraction = ggplot()+
            geom_point(data=df_Polygons_SnowFraction_GAM, aes(x=doy, y=SnowFraction, fill=Polygon, col=Polygon))+
            geom_line(data=df_Polygons_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict, col=Polygon)) +
            xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
            ylab("Fraction of snow-covered pixels") +
            facet_wrap(~NDSI_threshold, ncol=3)+
            theme_tom()

      #(4): Calculate at which day of year the SnowFraction value is equal to snowFraction_threshold for each Polygon using predictions from mod_gam

          #Create an empty dataframe to store snowfraction dates
          df_Polygon_Snowfraction <- data.frame(NDSI_threshold=character(),
                                                Snowfraction_threshold=character(),
                                                doy=numeric(),
                                                Polygon=character(),
                                                stringsAsFactors=FALSE)

          #Setup parallel processing
          numCores <- detectCores()
          cl <- makePSOCKcluster(numCores)
          registerDoSNOW(cl)

          #Loop through all Polygons in the df_Polygons_SnowFraction_GAM_predictions dataframe and extract snowfraction dates for each Polygon:
          for(i in unique(df_Polygons_SnowFraction_GAM_predictions$Polygon)){

            #For debugging
            #i=unique(df_Polygons_SnowFraction_GAM_predictions$Polygon)[2] #for debugging

            #Select dataset with GAM predictions for current Polygon
            df_Polygon_gam <- df_Polygons_SnowFraction_GAM_predictions[df_Polygons_SnowFraction_GAM_predictions$Polygon==i & !is.na(df_Polygons_SnowFraction_GAM_predictions$SnowFraction_gam_predict),]

            #If this results in an empty dataframe due to filtering of NAs then we return the original dataframe (required for threshold detection function below)
            if(nrow(df_Polygon_gam) < 1){df_Polygon_gam <- df_Polygons_SnowFraction_GAM_predictions[df_Polygons_SnowFraction_GAM_predictions$Polygon==i,]}

            #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector (per level of NDSI_threshold)
            results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                        pixelIDs_split = list(NDSI_threshold_vector), #levels of NDSI_threshold (input needs to be a list)
                                                        df_pixel_y = df_Polygon_gam, #dataframe containing GAM predictions
                                                        pixel_ID_column="NDSI_threshold", #Grouping column
                                                        y="SnowFraction_gam_predict", #response variable in GAM
                                                        x="doy", #predictor variable in GAM
                                                        pixel_gam_plots = FALSE, #Should GAM plots be created
                                                        y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated

            #Store dates of snowmelt per Polygon
            df_Polygon_SnowFraction_new <- results[[1]]
            df_Polygon_SnowFraction_new <- as.data.frame(do.call(rbind, df_Polygon_SnowFraction_new))
            colnames(df_Polygon_SnowFraction_new)[colnames(df_Polygon_SnowFraction_new)=="pixel_ID"] <- "NDSI_threshold"
            colnames(df_Polygon_SnowFraction_new)[colnames(df_Polygon_SnowFraction_new)=="x_threshold"] <- "doy"
            colnames(df_Polygon_SnowFraction_new)[colnames(df_Polygon_SnowFraction_new)=="y_threshold"] <- "Snowfraction_threshold"
            df_Polygon_SnowFraction_new$Polygon <- i

            #Add snowmelt data for current Polygon to general dataframe
            df_Polygon_Snowfraction <- rbind(df_Polygon_Snowfraction, df_Polygon_SnowFraction_new)

          }

          #Turn parallel processing off and run sequentially again after this point
          stopCluster(cl)
          registerDoSEQ()

          #Change column Polygon to a factor:
          df_Polygon_Snowfraction$Polygon <- as.factor(as.character(df_Polygon_Snowfraction$Polygon))

          #Add LocationID and year
          df_Polygon_Snowfraction$location <- area_name
          df_Polygon_Snowfraction$year <- year_ID

          #Save dates of snowmelt per polygon per SnowFraction threshold as a .csv file
          #write.csv(df_Polygon_Snowfraction, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Snowmelt_Snowfraction.csv"), row.names = FALSE)

      #(5): Create a separate plot with GAM predictions per polygon and NDSI_threshold:

          #Create an empty list to store plots
          list_plots_snowfraction <- vector('list', length(unique(df_Polygons_SnowFraction_GAM_predictions$NDSI_threshold)))

          #Loop through all levels of NDSI_threshold
          for(i in unique(df_Polygons_SnowFraction_GAM_predictions$NDSI_threshold)){

            #i=unique(df_Polygons_SnowFraction_GAM_predictions$NDSI_threshold)[1]

            #Create an index variable for parameter i
            i_index <- which( unique(df_Polygons_SnowFraction_GAM_predictions$NDSI_threshold) == i)

            #Store i as a character (for naming in output files)
            NDSI_threshold_char <- gsub("\\.", "_", as.character(i))

            #Loop through all Polygons:
            for(j in unique(df_Polygons_SnowFraction_GAM_predictions$Polygon)){

              #j=unique(df_Polygons_SnowFraction_GAM_predictions$Polygon)[1]

              #Create an index variable for parameter j
              j_index <- which( unique(df_Polygons_SnowFraction_GAM_predictions$Polygon) == j)

              #Select datasets for current NDSI_threshold and Polygon:
              df_Polygon_SnowFraction_GAM <- df_Polygons_SnowFraction_GAM[df_Polygons_SnowFraction_GAM$NDSI_threshold==i & df_Polygons_SnowFraction_GAM$Polygon==j,]
              df_Polygon_SnowFraction_GAM_predictions <- df_Polygons_SnowFraction_GAM_predictions[df_Polygons_SnowFraction_GAM_predictions$NDSI_threshold==i & df_Polygons_SnowFraction_GAM_predictions$Polygon==j,]
              df_Polygon_SnowfractionDate <- df_Polygon_Snowfraction[df_Polygon_Snowfraction$NDSI_threshold==i & df_Polygon_Snowfraction$Polygon==j,]

              #Create plot for current NDSI_threshold and Polygon and store it in list_plots_snowfraction:
              list_plots_snowfraction[[i_index]][[j_index]] <- ggplot()+
                geom_point(data=df_Polygon_SnowFraction_GAM[df_Polygon_SnowFraction_GAM$outliers==FALSE,], aes(x=doy, y=SnowFraction))+
                geom_point(data=df_Polygon_SnowFraction_GAM[df_Polygon_SnowFraction_GAM$outliers==TRUE,], aes(x=doy, y=SnowFraction), col="black", pch=16, alpha=0.2)+
                geom_line(data=df_Polygon_SnowFraction_GAM_predictions, aes(x=doy, y=SnowFraction_gam_predict), col="#1620de", lwd=1.25)+
                geom_point(data=df_Polygon_SnowfractionDate[!is.na(df_Polygon_SnowfractionDate$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)+
                geom_vline(xintercept = 150, colour="grey", lty=2)+
                geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
                xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                ylab("Fraction of snow-covered pixels") +
                ggtitle(paste0("Year: ", year_ID, ", Polygon: ", j, ", NDSI: ", i))+
                theme_tom()

               }

            #Plot SnowFraction and model predictions in a separate plot per Polygon per NDSI_threshold
            # plots_snowfraction <- list_plots_snowfraction[i_index]
            # plots_per_page = 25
            # plots_snowfraction <- lapply(plots_snowfraction, function(x){split(x, ceiling(seq_along(plots_snowfraction[[1]])/plots_per_page))})
            # plots_snowfraction <- unname(unlist(plots_snowfraction, recursive = F))
            # pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_NDSI_threshold_", NDSI_threshold_char, "_Polygons_Plot_Snowfraction.pdf"), width=20, height=16, onefile = TRUE)
            # for (k in seq(length(plots_snowfraction))) { do.call("grid.arrange", plots_snowfraction[[k]]) }
            # dev.off()

          }

        }


########################################################################################################################################################################################

  #(II): AVERAGE BAND VALUES - Fit a Generalized Additive Model (GAM) through the average band values within each polygon
    if("avg_NDSI" %in% method){

########################################################################################################################################################################################

      #(A): FSC_Gascoin2020 - Fit a Generalized Additive Model (GAM) through the FSC_Gascoin2020 data within each polygon
      if("FSC_Gascoin2020" %in% method_output){

        #(A.0) Print message
          cat("\n")
          print("--------------------------------------------------------------------------------------------------------------------------")
          print(paste0("FIT A GAM THROUGH THE AVERAGE FSC_GASCOIN2020 DATA PER POLYGON"))
          print("--------------------------------------------------------------------------------------------------------------------------")

        #(A.1) Create an empty dataframe
          df_Polygons_FSC_Gascoin2020_GAM <- data.frame(FSC_Gascoin2020=numeric(),
                                                        Date=factor(),
                                                        doy=numeric(),
                                                        Polygon=factor(),
                                                        outliers=logical())

          df_Polygons_FSC_Gascoin2020_GAM_predictions <- data.frame(Polygon=character(),
                                                                    doy=numeric(),
                                                                    FSC_Gascoin2020_gam_predict=numeric(),
                                                                    stringsAsFactors=FALSE)

        #(A.2) Loop through all Polygons and fit a separate GAM (with sequential outlier removal) to the Polygon-specific FSC_Gascoin2020 data
          for(i in unique(df_Polygons_BandValues$Polygon)){

            #For debugging
            #i=unique(df_Polygons_BandValues$Polygon)[2]

            #Select Polygon-specific subset of data:
            df_Polygon_FSC_Gascoin2020_GAM_new <- df_Polygons_BandValues[df_Polygons_BandValues$Polygon==i &
                                                                         !is.na(df_Polygons_BandValues$FSC_Gascoin2020),
                                                                         c("FSC_Gascoin2020", "Date", "doy", "Polygon")]

            #If there are at least two datapoints, then continue to fit the GAM
            if(nrow(df_Polygon_FSC_Gascoin2020_GAM_new) > 1){

              #Fit a gam through the Polygon-specific FSC_Gascoin2020 ~ doy data and emply sequential outlier removal
              df_Polygon_FSC_Gascoin2020_GAM_new <- f_gam_SeqRemOutliers(data=df_Polygon_FSC_Gascoin2020_GAM_new, y="FSC_Gascoin2020", x="doy", outlier_removal=outlier_removal,
                                                                         outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                                         default_k=gam_k_outlier)

              #Sort df_Polygon_FSC_Gascoin2020_GAM_new by doy:
              df_Polygon_FSC_Gascoin2020_GAM_new <- df_Polygon_FSC_Gascoin2020_GAM_new[order(df_Polygon_FSC_Gascoin2020_GAM_new$doy),]

              #Bind the Polygon-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
               df_Polygons_FSC_Gascoin2020_GAM <- rbind(df_Polygons_FSC_Gascoin2020_GAM, df_Polygon_FSC_Gascoin2020_GAM_new)

              #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

                #Refit GAM through data
                 index <- which(df_Polygon_FSC_Gascoin2020_GAM_new$outliers==FALSE)
                 mod_gam <- with(df_Polygon_FSC_Gascoin2020_GAM_new[index,], mgcv::gam(FSC_Gascoin2020 ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

                #Use gam to make predictions on a more detailed interval
                 df_Polygon_FSC_Gascoin2020_GAM_predictions_new <- data.frame(Polygon=i, doy=seq(min(df_Polygon_FSC_Gascoin2020_GAM_new$doy), max(df_Polygon_FSC_Gascoin2020_GAM_new$doy), 0.01))
                 df_Polygon_FSC_Gascoin2020_GAM_predictions_new$FSC_Gascoin2020_gam_predict <- stats::predict(mod_gam, newdata=df_Polygon_FSC_Gascoin2020_GAM_predictions_new, type="response")
                 df_Polygon_FSC_Gascoin2020_GAM_predictions_new <- df_Polygon_FSC_Gascoin2020_GAM_predictions_new[order(df_Polygon_FSC_Gascoin2020_GAM_predictions_new$doy),]

                #Add predictions to df_Polygons_FSC_Gascoin2020_GAM_predictions dataframe:
                 df_Polygons_FSC_Gascoin2020_GAM_predictions <- rbind(df_Polygons_FSC_Gascoin2020_GAM_predictions, df_Polygon_FSC_Gascoin2020_GAM_predictions_new)

               }

            #If there are less than two datapoints, then append empty dataframes
            if(nrow(df_Polygon_FSC_Gascoin2020_GAM_new) < 2){

              #Bind the empty dataframe to the dataframe containing all dataframes from previous iterations
              if(nrow(df_Polygon_FSC_Gascoin2020_GAM_new) < 1){df_Polygon_FSC_Gascoin2020_GAM_new <- data.frame(FSC_Gascoin2020=NA, Date=NA, doy=NA, Polygon=i, outliers=NA)}
              if(nrow(df_Polygon_FSC_Gascoin2020_GAM_new) > 0){df_Polygon_FSC_Gascoin2020_GAM_new$outliers <- NA}
              df_Polygons_FSC_Gascoin2020_GAM <- rbind(df_Polygons_FSC_Gascoin2020_GAM, df_Polygon_FSC_Gascoin2020_GAM_new)

              #Add empty GAM predictions to df_Polygons_FSC_Gascoin2020_GAM_predictions dataframe:
              df_Polygon_FSC_Gascoin2020_GAM_predictions_new <- df_Polygon_FSC_Gascoin2020_GAM_new[,c("Polygon", "doy")]
              df_Polygon_FSC_Gascoin2020_GAM_predictions_new$FSC_Gascoin2020_gam_predict <- NA
              df_Polygons_FSC_Gascoin2020_GAM_predictions <- rbind(df_Polygons_FSC_Gascoin2020_GAM_predictions, df_Polygon_FSC_Gascoin2020_GAM_predictions_new)

            }

            }

          #Change Polygon column to factor
          df_Polygons_FSC_Gascoin2020_GAM$Polygon <- as.factor(as.character(df_Polygons_FSC_Gascoin2020_GAM$Polygon))
          df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon <- as.factor(as.character(df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon))

          #Add LocationID and year
          df_Polygons_FSC_Gascoin2020_GAM_predictions$location <- area_name
          df_Polygons_FSC_Gascoin2020_GAM_predictions$year <- year_ID

          #Save dataframe with GAM fits for FSC_Gascoin2020
          #write.csv(df_Polygons_FSC_Gascoin2020_GAM, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_FSCGascoin2020_GAM.csv"), row.names = FALSE)
          #write.csv(df_Polygons_FSC_Gascoin2020_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_GAM_Predictions_FSCGascoin2020.csv"), row.names = FALSE)

        #(A.3) Plot the raw FSC_Gascoin2020 datapoints and gam predictions for each Polygon:

           #Plot FSC_Gascoin2020 and model predictions for all Polygons in a single plot
            p_Polygon_FSC_Gascoin2020 = ggplot()+
              geom_point(data=df_Polygons_FSC_Gascoin2020_GAM, aes(x=doy, y=FSC_Gascoin2020, fill=Polygon, col=Polygon))+
              geom_line(data=df_Polygons_FSC_Gascoin2020_GAM_predictions, aes(x=doy, y=FSC_Gascoin2020_gam_predict, col=Polygon)) +
              geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
              xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
              ylab("FSC_Gascoin2020") +
              theme_tom()

        #(A.4): Calculate at which day of year the FSC_Gascoin2020 value is equal to Snowfraction_threshold for each Polygon using predictions from mod_gam

            #Setup parallel processing
            numCores <- detectCores()
            cl <- makePSOCKcluster(numCores)
            registerDoSNOW(cl)

            #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector
            results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                        pixelIDs_split = list(unique(df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon)), #levels of Polygon (input needs to be a list)
                                                        df_pixel_y = df_Polygons_FSC_Gascoin2020_GAM_predictions, #dataframe containing GAM predictions
                                                        pixel_ID_column="Polygon", #Grouping column
                                                        y="FSC_Gascoin2020_gam_predict", #response variable in GAM
                                                        x="doy", #predictor variable in GAM
                                                        pixel_gam_plots = FALSE, #Should GAM plots be created
                                                        y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated

            #Turn parallel processing off and run sequentially again after this point
            stopCluster(cl)
            registerDoSEQ()

            #Store dates of snowmelt per Polygon
            df_Polygon_FSC_Gascoin2020 <- results[[1]]
            df_Polygon_FSC_Gascoin2020 <- as.data.frame(do.call(rbind, df_Polygon_FSC_Gascoin2020))
            colnames(df_Polygon_FSC_Gascoin2020)[colnames(df_Polygon_FSC_Gascoin2020)=="pixel_ID"] <- "Polygon"
            colnames(df_Polygon_FSC_Gascoin2020)[colnames(df_Polygon_FSC_Gascoin2020)=="x_threshold"] <- "doy"
            colnames(df_Polygon_FSC_Gascoin2020)[colnames(df_Polygon_FSC_Gascoin2020)=="y_threshold"] <- "Snowfraction_threshold"

            #Change column Polygon to a factor:
            df_Polygon_FSC_Gascoin2020$Polygon <- as.factor(as.character(df_Polygon_FSC_Gascoin2020$Polygon))

            #Add LocationID and year
            df_Polygon_FSC_Gascoin2020$location <- area_name
            df_Polygon_FSC_Gascoin2020$year <- year_ID

            #Save dates of snowmelt per polygon per FSC_Gascoin2020 threshold as a .csv file
            #write.csv(df_Polygon_FSC_Gascoin2020, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Snowmelt_FSCGascoin2020.csv"), row.names = FALSE)

        #(A.5): Create a separate plot with GAM predictions per polygon

            #Create an empty list to store plots
            list_plots_fsc_gascoin2020 <- vector('list', length(unique(df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon)))

            #Loop through all Polygons
            for(i in unique(df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon)){

              #i=unique(df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon)[1]

              #Create an index variable for parameter i
              i_index <- which( unique(df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon) == i)

              #select datasets for current Snowfraction_threshold and Polygon:
              df_Polygon_FSC_Gascoin2020_GAM <- df_Polygons_FSC_Gascoin2020_GAM[df_Polygons_FSC_Gascoin2020_GAM$Polygon==i,]
              df_Polygon_FSC_Gascoin2020_GAM_predictions <- df_Polygons_FSC_Gascoin2020_GAM_predictions[df_Polygons_FSC_Gascoin2020_GAM_predictions$Polygon==i,]
              df_Polygon_FSC_Gascoin2020_snowmelt <- df_Polygon_FSC_Gascoin2020[df_Polygon_FSC_Gascoin2020$Polygon==i,]

              #create plot for current Polygon and store it in list_plots_fsc_gascoin2020:
              list_plots_fsc_gascoin2020[[i_index]] <- ggplot()+
                geom_point(data=df_Polygon_FSC_Gascoin2020_GAM[df_Polygon_FSC_Gascoin2020_GAM$outliers==FALSE,], aes(x=doy, y=FSC_Gascoin2020))+
                geom_point(data=df_Polygon_FSC_Gascoin2020_GAM[df_Polygon_FSC_Gascoin2020_GAM$outliers==TRUE,], aes(x=doy, y=FSC_Gascoin2020), col="black", pch=16, alpha=0.2)+
                geom_line(data=df_Polygon_FSC_Gascoin2020_GAM_predictions, aes(x=doy, y=FSC_Gascoin2020_gam_predict), col="#1620de", lwd=1.25)+
                geom_point(data=df_Polygon_FSC_Gascoin2020_snowmelt[!is.na(df_Polygon_FSC_Gascoin2020_snowmelt$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)+
                geom_vline(xintercept = 150, colour="grey", lty=2)+
                geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
                xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                ylab("FSC_Gascoin2020") +
                ggtitle(paste0("Year: ", year_ID, ", Polygon: ", i))+
                theme_tom()

              }

            #Plot FSC_Gascoin2020 and model predictions in a separate plot per Polygon
            # plots_fsc_gascoin2020 <- list(list_plots_fsc_gascoin2020)
            # plots_per_page = 25
            # plots_fsc_gascoin2020 <- lapply(plots_fsc_gascoin2020, function(x){split(x, ceiling(seq_along(plots_fsc_gascoin2020[[1]])/plots_per_page))})
            # plots_fsc_gascoin2020 <- unname(unlist(plots_fsc_gascoin2020, recursive = F))
            # pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Plot_FSCGascoin2020.pdf"), width=20, height=16, onefile = TRUE)
            # for (k in seq(length(plots_fsc_gascoin2020))) { do.call("grid.arrange", plots_fsc_gascoin2020[[k]]) }
            # dev.off()

      }

########################################################################################################################################################################################

      #(B): FSC_Aalstad2020 - Fit a Generalized Additive Model (GAM) through the FSC_Aalstad2020 data within each polygon
      if("FSC_Aalstad2020" %in% method_output){

        #(B.0) Print message
          cat("\n")
          print("--------------------------------------------------------------------------------------------------------------------------")
          print(paste0("FIT A GAM THROUGH THE AVERAGE FSC_AALSTAD2020 DATA PER POLYGON"))
          print("--------------------------------------------------------------------------------------------------------------------------")

        #(B.1) Create an empty dataframe
          df_Polygons_FSC_Aalstad2020_GAM <- data.frame(FSC_Aalstad2020=numeric(),
                                                        Date=factor(),
                                                        doy=numeric(),
                                                        Polygon=factor(),
                                                        outliers=logical())

          df_Polygons_FSC_Aalstad2020_GAM_predictions <- data.frame(Polygon=character(),
                                                                    doy=numeric(),
                                                                    FSC_Aalstad2020_gam_predict=numeric(),
                                                                    stringsAsFactors=FALSE)

        #(B.2) Loop through all Polygons and fit a separate GAM (with sequential outlier removal) to the Polygon-specific FSC_Aalstad2020 data
          for(i in unique(df_Polygons_BandValues$Polygon)){

            #For debugging
            #i=unique(df_Polygons_BandValues$Polygon)[1]

            #Select Polygon-specific subset of data:
            df_Polygon_FSC_Aalstad2020_GAM_new <- df_Polygons_BandValues[df_Polygons_BandValues$Polygon==i &
                                                                         !is.na(df_Polygons_BandValues$FSC_Aalstad2020),
                                                                         c("FSC_Aalstad2020", "Date", "doy", "Polygon")]

            #If there are at least two datapoints, then continue to fit the GAM
            if(nrow(df_Polygon_FSC_Aalstad2020_GAM_new) > 1){

              #Fit a gam through the Polygon-specific FSC_Aalstad2020 ~ doy data and emply sequential outlier removal
              df_Polygon_FSC_Aalstad2020_GAM_new <- f_gam_SeqRemOutliers(data=df_Polygon_FSC_Aalstad2020_GAM_new, y="FSC_Aalstad2020", x="doy", outlier_removal=outlier_removal,
                                                                         outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                                         default_k=gam_k_outlier)

              #Sort df_Polygon_FSC_Aalstad2020_GAM_new by doy:
              df_Polygon_FSC_Aalstad2020_GAM_new <- df_Polygon_FSC_Aalstad2020_GAM_new[order(df_Polygon_FSC_Aalstad2020_GAM_new$doy),]

              #Bind the Polygon-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
              df_Polygons_FSC_Aalstad2020_GAM <- rbind(df_Polygons_FSC_Aalstad2020_GAM, df_Polygon_FSC_Aalstad2020_GAM_new)

              #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

                #Refit GAM through data
                 index <- which(df_Polygon_FSC_Aalstad2020_GAM_new$outliers==FALSE)
                 mod_gam <- with(df_Polygon_FSC_Aalstad2020_GAM_new[index,], mgcv::gam(FSC_Aalstad2020 ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

                #Use gam to make predictions on a more detailed interval
                 df_Polygon_FSC_Aalstad2020_GAM_predictions_new <- data.frame(Polygon=i, doy=seq(min(df_Polygon_FSC_Aalstad2020_GAM_new$doy), max(df_Polygon_FSC_Aalstad2020_GAM_new$doy), 0.01))
                 df_Polygon_FSC_Aalstad2020_GAM_predictions_new$FSC_Aalstad2020_gam_predict <- stats::predict(mod_gam, newdata=df_Polygon_FSC_Aalstad2020_GAM_predictions_new, type="response")
                 df_Polygon_FSC_Aalstad2020_GAM_predictions_new <- df_Polygon_FSC_Aalstad2020_GAM_predictions_new[order(df_Polygon_FSC_Aalstad2020_GAM_predictions_new$doy),]

                #Add predictions to df_Polygons_FSC_Aalstad2020_GAM_predictions dataframe:
                 df_Polygons_FSC_Aalstad2020_GAM_predictions <- rbind(df_Polygons_FSC_Aalstad2020_GAM_predictions, df_Polygon_FSC_Aalstad2020_GAM_predictions_new)

              }

            #If there are less than two datapoints, then append empty dataframes
            if(nrow(df_Polygon_FSC_Aalstad2020_GAM_new) < 2){

              #Bind the empty dataframe to the dataframe containing all dataframes from previous iterations
              if(nrow(df_Polygon_FSC_Aalstad2020_GAM_new) < 1){df_Polygon_FSC_Aalstad2020_GAM_new <- data.frame(FSC_Aalstad2020=NA, Date=NA, doy=NA, Polygon=i, outliers=NA)}
              if(nrow(df_Polygon_FSC_Aalstad2020_GAM_new) > 0){df_Polygon_FSC_Aalstad2020_GAM_new$outliers <- NA}
              df_Polygons_FSC_Aalstad2020_GAM <- rbind(df_Polygons_FSC_Aalstad2020_GAM, df_Polygon_FSC_Aalstad2020_GAM_new)

              #Add empty GAM predictions to df_Polygons_FSC_Aalstad2020_GAM_predictions dataframe:
              df_Polygon_FSC_Aalstad2020_GAM_predictions_new <- df_Polygon_FSC_Aalstad2020_GAM_new[,c("Polygon", "doy")]
              df_Polygon_FSC_Aalstad2020_GAM_predictions_new$FSC_Aalstad2020_gam_predict <- NA
              df_Polygons_FSC_Aalstad2020_GAM_predictions <- rbind(df_Polygons_FSC_Aalstad2020_GAM_predictions, df_Polygon_FSC_Aalstad2020_GAM_predictions_new)

            }

            }

          #Change Polygon column to factor
          df_Polygons_FSC_Aalstad2020_GAM$Polygon <- as.factor(as.character(df_Polygons_FSC_Aalstad2020_GAM$Polygon))
          df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon <- as.factor(as.character(df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon))

          #Add LocationID and year
          df_Polygons_FSC_Aalstad2020_GAM_predictions$location <- area_name
          df_Polygons_FSC_Aalstad2020_GAM_predictions$year <- year_ID

          #Save dataframe with GAM fits for FSC_Aalstad2020
          #write.csv(df_Polygons_FSC_Aalstad2020_GAM, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_FSCAalstad2020_GAM.csv"), row.names = FALSE)
          #write.csv(df_Polygons_FSC_Aalstad2020_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_GAM_Predictions_FSCAalstad2020.csv"), row.names = FALSE)

        #(B.3) Plot the raw FSC_Aalstad2020 datapoints and gam predictions for each Polygon:

           #Plot FSC_Aalstad2020 and model predictions for all Polygons in a single plot
            p_Polygon_FSC_Aalstad2020 = ggplot()+
              geom_point(data=df_Polygons_FSC_Aalstad2020_GAM, aes(x=doy, y=FSC_Aalstad2020, fill=Polygon, col=Polygon))+
              geom_line(data=df_Polygons_FSC_Aalstad2020_GAM_predictions, aes(x=doy, y=FSC_Aalstad2020_gam_predict, col=Polygon)) +
              geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
              xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
              ylab("FSC_Aalstad2020") +
              theme_tom()

        #(B.4): Calculate at which day of year the FSC_Aalstad2020 value is equal to Snowfraction_threshold for each Polygon using predictions from mod_gam

            #Setup parallel processing
            numCores <- detectCores()
            cl <- makePSOCKcluster(numCores)
            registerDoSNOW(cl)

            #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in Snowfraction_threshold_vector
            results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                        pixelIDs_split = list(unique(df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon)), #levels of Polygon (input needs to be a list)
                                                        df_pixel_y = df_Polygons_FSC_Aalstad2020_GAM_predictions, #dataframe containing GAM predictions
                                                        pixel_ID_column="Polygon", #Grouping column
                                                        y="FSC_Aalstad2020_gam_predict", #response variable in GAM
                                                        x="doy", #predictor variable in GAM
                                                        pixel_gam_plots = FALSE, #Should GAM plots be created
                                                        y_threshold = Snowfraction_threshold_vector) #Which threshold values for 'y' should be calculated

            #Turn parallel processing off and run sequentially again after this point
            stopCluster(cl)
            registerDoSEQ()

            #Store dates of snowmelt per Polygon
            df_Polygon_FSC_Aalstad2020 <- results[[1]]
            df_Polygon_FSC_Aalstad2020 <- as.data.frame(do.call(rbind, df_Polygon_FSC_Aalstad2020))
            colnames(df_Polygon_FSC_Aalstad2020)[colnames(df_Polygon_FSC_Aalstad2020)=="pixel_ID"] <- "Polygon"
            colnames(df_Polygon_FSC_Aalstad2020)[colnames(df_Polygon_FSC_Aalstad2020)=="x_threshold"] <- "doy"
            colnames(df_Polygon_FSC_Aalstad2020)[colnames(df_Polygon_FSC_Aalstad2020)=="y_threshold"] <- "Snowfraction_threshold"

            #Change column Polygon to a factor:
            df_Polygon_FSC_Aalstad2020$Polygon <- as.factor(as.character(df_Polygon_FSC_Aalstad2020$Polygon))

            #Add LocationID and year
            df_Polygon_FSC_Aalstad2020$location <- area_name
            df_Polygon_FSC_Aalstad2020$year <- year_ID

            #Save dates of snowmelt per polygon per FSC_Aalstad2020 threshold as a .csv file
            #write.csv(df_Polygon_FSC_Aalstad2020, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Snowmelt_FSCAalstad2020.csv"), row.names = FALSE)

        #(B.5): Create a separate plot with GAM predictions per polygon

            #Create an empty list to store plots
            list_plots_fsc_aalstad2020 <- vector('list', length(unique(df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon)))

            #Loop through all Polygons
            for(i in unique(df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon)){

              #i=unique(df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon)[1]

              #Create an index variable for parameter i
              i_index <- which( unique(df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon) == i)

              #select datasets for current Snowfraction_threshold and Polygon:
              df_Polygon_FSC_Aalstad2020_GAM <- df_Polygons_FSC_Aalstad2020_GAM[df_Polygons_FSC_Aalstad2020_GAM$Polygon==i,]
              df_Polygon_FSC_Aalstad2020_GAM_predictions <- df_Polygons_FSC_Aalstad2020_GAM_predictions[df_Polygons_FSC_Aalstad2020_GAM_predictions$Polygon==i,]
              df_Polygon_FSC_Aalstad2020_snowmelt <- df_Polygon_FSC_Aalstad2020[df_Polygon_FSC_Aalstad2020$Polygon==i,]

              #create plot for current Polygon and store it in list_plots_fsc_aalstad2020:
              list_plots_fsc_aalstad2020[[i_index]] <- ggplot()+
                geom_point(data=df_Polygon_FSC_Aalstad2020_GAM[df_Polygon_FSC_Aalstad2020_GAM$outliers==FALSE,], aes(x=doy, y=FSC_Aalstad2020))+
                geom_point(data=df_Polygon_FSC_Aalstad2020_GAM[df_Polygon_FSC_Aalstad2020_GAM$outliers==TRUE,], aes(x=doy, y=FSC_Aalstad2020), col="black", pch=16, alpha=0.2)+
                geom_line(data=df_Polygon_FSC_Aalstad2020_GAM_predictions, aes(x=doy, y=FSC_Aalstad2020_gam_predict), col="#1620de", lwd=1.25)+
                geom_point(data=df_Polygon_FSC_Aalstad2020_snowmelt[!is.na(df_Polygon_FSC_Aalstad2020_snowmelt$doy),], aes(x=doy, y=Snowfraction_threshold), col="red", size=3)+
                geom_vline(xintercept = 150, colour="grey", lty=2)+
                geom_hline(yintercept = Snowfraction_threshold_vector, colour="grey", lty=2)+
                xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                ylab("FSC_Aalstad2020") +
                ggtitle(paste0("Year: ", year_ID, ", Polygon: ", i))+
                theme_tom()

              }

            #Plot FSC_Aalstad2020 and model predictions in a separate plot per Polygon
            # plots_fsc_aalstad2020 <- list(list_plots_fsc_aalstad2020)
            # plots_per_page = 25
            # plots_fsc_aalstad2020 <- lapply(plots_fsc_aalstad2020, function(x){split(x, ceiling(seq_along(plots_fsc_aalstad2020[[1]])/plots_per_page))})
            # plots_fsc_aalstad2020 <- unname(unlist(plots_fsc_aalstad2020, recursive = F))
            # pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Plot_FSCAalstad2020.pdf"), width=20, height=16, onefile = TRUE)
            # for (k in seq(length(plots_fsc_aalstad2020))) { do.call("grid.arrange", plots_fsc_aalstad2020[[k]]) }
            # dev.off()

      }

########################################################################################################################################################################################

      #(C): NDSI - Fit a Generalized Additive Model (GAM) through the NDSI data within each polygon
      if("NDSI" %in% method_output){

          #(C.0) Print message
            cat("\n")
            print("--------------------------------------------------------------------------------------------------------------------------")
            print(paste0("FIT A GAM THROUGH THE AVERAGE NDSI DATA PER POLYGON"))
            print("--------------------------------------------------------------------------------------------------------------------------")

          #(C.1) Create an empty dataframe
            df_Polygons_NDSI_GAM <- data.frame(NDSI=numeric(),
                                               Date=factor(),
                                               doy=numeric(),
                                               Polygon=factor(),
                                               outliers=logical())

            df_Polygons_NDSI_GAM_predictions <- data.frame(Polygon=character(),
                                                           doy=numeric(),
                                                           NDSI_gam_predict=numeric(),
                                                           stringsAsFactors=FALSE)

          #(C.2) Loop through all Polygons and fit a separate GAM (with sequential outlier removal) to the Polygon-specific NDSI data
            for(i in unique(df_Polygons_BandValues$Polygon)){

              #For debugging
              #i=unique(df_Polygons_BandValues$Polygon)[2]

              #Select Polygon-specific subset of data:
              df_Polygon_NDSI_GAM_new <- df_Polygons_BandValues[df_Polygons_BandValues$Polygon==i &
                                                                !is.na(df_Polygons_BandValues$NDSI),
                                                                c("NDSI", "Date", "doy", "Polygon")]

              #If there are at least two datapoints, then continue to fit the GAM
              if(nrow(df_Polygon_NDSI_GAM_new) > 1){

                #Fit a gam through the Polygon-specific NDSI ~ doy data and emply sequential outlier removal
                df_Polygon_NDSI_GAM_new <- f_gam_SeqRemOutliers(data=df_Polygon_NDSI_GAM_new, y="NDSI", x="doy", outlier_removal=outlier_removal,
                                                                 outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                                 default_k=gam_k_outlier)

                #Sort df_Polygon_NDSI_GAM_new by doy:
                df_Polygon_NDSI_GAM_new <- df_Polygon_NDSI_GAM_new[order(df_Polygon_NDSI_GAM_new$doy),]

                #Bind the Polygon-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
                 df_Polygons_NDSI_GAM <- rbind(df_Polygons_NDSI_GAM, df_Polygon_NDSI_GAM_new)

                #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

                  #Refit GAM through data
                   index <- which(df_Polygon_NDSI_GAM_new$outliers==FALSE)
                   mod_gam <- with(df_Polygon_NDSI_GAM_new[index,], mgcv::gam(NDSI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

                  #Use gam to make predictions on a more detailed interval
                   df_Polygon_NDSI_GAM_predictions_new <- data.frame(Polygon=i, doy=seq(min(df_Polygon_NDSI_GAM_new$doy), max(df_Polygon_NDSI_GAM_new$doy), 0.01))
                   df_Polygon_NDSI_GAM_predictions_new$NDSI_gam_predict <- stats::predict(mod_gam, newdata=df_Polygon_NDSI_GAM_predictions_new, type="response")
                   df_Polygon_NDSI_GAM_predictions_new <- df_Polygon_NDSI_GAM_predictions_new[order(df_Polygon_NDSI_GAM_predictions_new$doy),]

                  #Add predictions to df_Polygons_NDSI_GAM_predictions dataframe:
                   df_Polygons_NDSI_GAM_predictions <- rbind(df_Polygons_NDSI_GAM_predictions, df_Polygon_NDSI_GAM_predictions_new)

                }

              #If there are less than two datapoints, then append empty dataframes
              if(nrow(df_Polygon_NDSI_GAM_new) < 2){

                #Bind the empty dataframe to the dataframe containing all dataframes from previous iterations
                if(nrow(df_Polygon_NDSI_GAM_new) < 1){df_Polygon_NDSI_GAM_new <- data.frame(NDSI=NA, Date=NA, doy=NA, Polygon=i, outliers=NA)}
                if(nrow(df_Polygon_NDSI_GAM_new) > 0){df_Polygon_NDSI_GAM_new$outliers <- NA}
                df_Polygons_NDSI_GAM <- rbind(df_Polygons_NDSI_GAM, df_Polygon_NDSI_GAM_new)

                #Add empty GAM predictions to df_Polygons_NDSI_GAM_predictions dataframe:
                df_Polygon_NDSI_GAM_predictions_new <- df_Polygon_NDSI_GAM_new[,c("Polygon", "doy")]
                df_Polygon_NDSI_GAM_predictions_new$NDSI_gam_predict <- NA
                df_Polygons_NDSI_GAM_predictions <- rbind(df_Polygons_NDSI_GAM_predictions, df_Polygon_NDSI_GAM_predictions_new)

                }

              }

            #Change Polygon column to factor
            df_Polygons_NDSI_GAM$Polygon <- as.factor(as.character(df_Polygons_NDSI_GAM$Polygon))
            df_Polygons_NDSI_GAM_predictions$Polygon <- as.factor(as.character(df_Polygons_NDSI_GAM_predictions$Polygon))

            #Add LocationID and year
            df_Polygons_NDSI_GAM_predictions$location <- area_name
            df_Polygons_NDSI_GAM_predictions$year <- year_ID

            #Save dataframe with GAM fits for NDSI
            #write.csv(df_Polygons_NDSI_GAM, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_NDSI_GAM.csv"), row.names = FALSE)
            #write.csv(df_Polygons_NDSI_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_GAM_Predictions_NDSI.csv"), row.names = FALSE)

          #(C.3) Plot the raw NDSI datapoints and gam predictions for each Polygon:

             #Plot NDSI and model predictions for all Polygons in a single plot
              p_Polygon_NDSI = ggplot()+
                geom_point(data=df_Polygons_NDSI_GAM, aes(x=doy, y=NDSI, fill=Polygon, col=Polygon))+
                geom_line(data=df_Polygons_NDSI_GAM_predictions, aes(x=doy, y=NDSI_gam_predict, col=Polygon)) +
                geom_hline(yintercept = NDSI_threshold_vector, colour="grey", lty=2)+
                xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                ylab("NDSI") +
                theme_tom()

          #(C.4): Calculate at which day of year the NDSI value is equal to NDSI_threshold for each Polygon using predictions from mod_gam

              #Setup parallel processing
              numCores <- detectCores()
              cl <- makePSOCKcluster(numCores)
              registerDoSNOW(cl)

              #Use the function f_detect_threshold_date_parallel to extract the moments the GAM predictions cross the thresholds in NDSI_threshold_vector
              results <- f_detect_threshold_date_parallel(subset=1, #data subset (not relevant here, set to 1)
                                                          pixelIDs_split = list(unique(df_Polygons_NDSI_GAM_predictions$Polygon)), #levels of NDSI_threshold (input needs to be a list)
                                                          df_pixel_y = df_Polygons_NDSI_GAM_predictions, #dataframe containing GAM predictions
                                                          pixel_ID_column="Polygon", #Grouping column
                                                          y="NDSI_gam_predict", #response variable in GAM
                                                          x="doy", #predictor variable in GAM
                                                          pixel_gam_plots = FALSE, #Should GAM plots be created
                                                          y_threshold = NDSI_threshold_vector) #Which threshold values for 'y' should be calculated

              #Turn parallel processing off and run sequentially again after this point
              stopCluster(cl)
              registerDoSEQ()

              #Store dates of snowmelt per Polygon
              df_Polygon_NDSI <- results[[1]]
              df_Polygon_NDSI <- as.data.frame(do.call(rbind, df_Polygon_NDSI))
              colnames(df_Polygon_NDSI)[colnames(df_Polygon_NDSI)=="pixel_ID"] <- "Polygon"
              colnames(df_Polygon_NDSI)[colnames(df_Polygon_NDSI)=="x_threshold"] <- "doy"
              colnames(df_Polygon_NDSI)[colnames(df_Polygon_NDSI)=="y_threshold"] <- "NDSI_threshold"

              #Change column Polygon to a factor:
              df_Polygon_NDSI$Polygon <- as.factor(as.character(df_Polygon_NDSI$Polygon))

              #Add LocationID and year
              df_Polygon_NDSI$location <- area_name
              df_Polygon_NDSI$year <- year_ID

              #Save dates of snowmelt per polygon per NDSI threshold as a .csv file
              #write.csv(df_Polygon_NDSI, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Snowmelt_NDSI.csv"), row.names = FALSE)

          #(C.5): Create a separate plot with GAM predictions per polygon

              #Create an empty list to store plots
              list_plots_ndsi <- vector('list', length(unique(df_Polygons_NDSI_GAM_predictions$Polygon)))

              #Loop through all Polygons
              for(i in unique(df_Polygons_NDSI_GAM_predictions$Polygon)){

                #i=unique(df_Polygons_NDSI_GAM_predictions$Polygon)[1]

                #Create an index variable for parameter i
                i_index <- which( unique(df_Polygons_NDSI_GAM_predictions$Polygon) == i)

                #select datasets for current NDSI_threshold and Polygon:
                df_Polygon_NDSI_GAM <- df_Polygons_NDSI_GAM[df_Polygons_NDSI_GAM$Polygon==i,]
                df_Polygon_NDSI_GAM_predictions <- df_Polygons_NDSI_GAM_predictions[df_Polygons_NDSI_GAM_predictions$Polygon==i,]
                df_Polygon_NDSI_snowmelt <- df_Polygon_NDSI[df_Polygon_NDSI$Polygon==i,]

                #create plot for current Polygon and store it in list_plots_ndsi:
                list_plots_ndsi[[i_index]] <- ggplot()+
                  geom_point(data=df_Polygon_NDSI_GAM[df_Polygon_NDSI_GAM$outliers==FALSE,], aes(x=doy, y=NDSI))+
                  geom_point(data=df_Polygon_NDSI_GAM[df_Polygon_NDSI_GAM$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
                  geom_line(data=df_Polygon_NDSI_GAM_predictions, aes(x=doy, y=NDSI_gam_predict), col="#1620de", lwd=1.25)+
                  geom_point(data=df_Polygon_NDSI_snowmelt[!is.na(df_Polygon_NDSI_snowmelt$doy),], aes(x=doy, y=NDSI_threshold), col="red", size=3)+
                  geom_vline(xintercept = 150, colour="grey", lty=2)+
                  geom_hline(yintercept = NDSI_threshold_vector, colour="grey", lty=2)+
                  xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                  ylab("NDSI") +
                  ggtitle(paste0("Year: ", year_ID, ", Polygon: ", i))+
                  theme_tom()

                }

              #Plot NDSI and model predictions in a separate plot per Polygon
              # plots_ndsi <- list(list_plots_ndsi)
              # plots_per_page = 25
              # plots_ndsi <- lapply(plots_ndsi, function(x){split(x, ceiling(seq_along(plots_ndsi[[1]])/plots_per_page))})
              # plots_ndsi <- unname(unlist(plots_ndsi, recursive = F))
              # pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Plot_NDSI.pdf"), width=20, height=16, onefile = TRUE)
              # for (k in seq(length(plots_ndsi))) { do.call("grid.arrange", plots_ndsi[[k]]) }
              # dev.off()

      }

########################################################################################################################################################################################

      #(D): NDVI - Fit a Generalized Additive Model (GAM) through the NDVI data within each polygon
      if("NDVI" %in% method_output){

          #(D.0) Print message
            cat("\n")
            print("--------------------------------------------------------------------------------------------------------------------------")
            print(paste0("FIT A GAM THROUGH THE AVERAGE NDVI DATA PER POLYGON"))
            print("--------------------------------------------------------------------------------------------------------------------------")

          #(D.1) Create an empty dataframe
            df_Polygons_NDVI_GAM <- data.frame(NDVI=numeric(),
                                               Date=factor(),
                                               doy=numeric(),
                                               Polygon=factor(),
                                               outliers=logical())

            df_Polygons_NDVI_GAM_predictions <- data.frame(Polygon=character(),
                                                           doy=numeric(),
                                                           NDVI_gam_predict=numeric(),
                                                           stringsAsFactors=FALSE)

          #(D.2) Loop through all Polygons and fit a separate GAM (with sequential outlier removal) to the Polygon-specific NDVI data
            for(i in unique(df_Polygons_BandValues$Polygon)){

              #For debugging
              #i=unique(df_Polygons_BandValues$Polygon)[1]

              #Select Polygon-specific subset of data:
              df_Polygon_NDVI_GAM_new <- df_Polygons_BandValues[df_Polygons_BandValues$Polygon==i &
                                                                !is.na(df_Polygons_BandValues$NDVI),
                                                                c("NDVI", "Date", "doy", "Polygon")]

              #If there are at least two datapoints, then continue to fit the GAM
              if(nrow(df_Polygon_NDVI_GAM_new) > 1){

                #Fit a gam through the Polygon-specific NDVI ~ doy data and emply sequential outlier removal
                df_Polygon_NDVI_GAM_new <- f_gam_SeqRemOutliers(data=df_Polygon_NDVI_GAM_new, y="NDVI", x="doy", outlier_removal=outlier_removal,
                                                                 outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                                 default_k=gam_k_outlier)

                #Sort df_Polygon_NDVI_GAM_new by doy:
                df_Polygon_NDVI_GAM_new <- df_Polygon_NDVI_GAM_new[order(df_Polygon_NDVI_GAM_new$doy),]

                #Bind the Polygon-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
                 df_Polygons_NDVI_GAM <- rbind(df_Polygons_NDVI_GAM, df_Polygon_NDVI_GAM_new)

                #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

                  #Refit GAM through data
                   index <- which(df_Polygon_NDVI_GAM_new$outliers==FALSE)
                   mod_gam <- with(df_Polygon_NDVI_GAM_new[index,], mgcv::gam(NDVI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

                  #Use gam to make predictions on a more detailed (1-day) day of year interval
                   df_Polygon_NDVI_GAM_predictions_new <- data.frame(Polygon=i, doy=seq(min(df_Polygon_NDVI_GAM_new$doy), max(df_Polygon_NDVI_GAM_new$doy), 0.01))
                   df_Polygon_NDVI_GAM_predictions_new$NDVI_gam_predict <- stats::predict(mod_gam, newdata=df_Polygon_NDVI_GAM_predictions_new, type="response")
                   df_Polygon_NDVI_GAM_predictions_new <- df_Polygon_NDVI_GAM_predictions_new[order(df_Polygon_NDVI_GAM_predictions_new$doy),]

                  #Add predictions to df_Polygons_NDVI_GAM_predictions dataframe:
                   df_Polygons_NDVI_GAM_predictions <- rbind(df_Polygons_NDVI_GAM_predictions, df_Polygon_NDVI_GAM_predictions_new)

                }

              #If there are less than two datapoints, then append empty dataframes
              if(nrow(df_Polygon_NDVI_GAM_new) < 2){

                #Bind the empty dataframe to the dataframe containing all dataframes from previous iterations
                if(nrow(df_Polygon_NDVI_GAM_new) < 1){df_Polygon_NDVI_GAM_new <- data.frame(NDVI=NA, Date=NA, doy=NA, Polygon=i, outliers=NA)}
                if(nrow(df_Polygon_NDVI_GAM_new) > 0){df_Polygon_NDVI_GAM_new$outliers <- NA}
                df_Polygons_NDVI_GAM <- rbind(df_Polygons_NDVI_GAM, df_Polygon_NDVI_GAM_new)

                #Add empty GAM predictions to df_Polygons_NDVI_GAM_predictions dataframe:
                df_Polygon_NDVI_GAM_predictions_new <- df_Polygon_NDVI_GAM_new[,c("Polygon", "doy")]
                df_Polygon_NDVI_GAM_predictions_new$NDVI_gam_predict <- NA
                df_Polygons_NDVI_GAM_predictions <- rbind(df_Polygons_NDVI_GAM_predictions, df_Polygon_NDVI_GAM_predictions_new)

                }

              }

            #Change Polygon column to factor
            df_Polygons_NDVI_GAM$Polygon <- as.factor(as.character(df_Polygons_NDVI_GAM$Polygon))
            df_Polygons_NDVI_GAM_predictions$Polygon <- as.factor(as.character(df_Polygons_NDVI_GAM_predictions$Polygon))

            #Add LocationID and year
            df_Polygons_NDVI_GAM_predictions$location <- area_name
            df_Polygons_NDVI_GAM_predictions$year <- year_ID

            #Save dataframe with GAM fits for NDVI
            #write.csv(df_Polygons_NDVI_GAM, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_NDVI_GAM.csv"), row.names = FALSE)
            #write.csv(df_Polygons_NDVI_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_GAM_Predictions_NDVI.csv"), row.names = FALSE)

          #(D.3) Plot the raw NDVI datapoints and gam predictions for each Polygon:

             #Plot NDVI and model predictions for all Polygons in a single plot
              p_Polygon_NDVI = ggplot()+
                geom_point(data=df_Polygons_NDVI_GAM, aes(x=doy, y=NDVI, fill=Polygon, col=Polygon))+
                geom_line(data=df_Polygons_NDVI_GAM_predictions, aes(x=doy, y=NDVI_gam_predict, col=Polygon)) +
                xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                ylab("NDVI") +
                theme_tom()

          #(D.4): Create a separate plot with GAM predictions per polygon

              #Create an empty list to store plots
              list_plots_ndvi <- vector('list', length(unique(df_Polygons_NDVI_GAM_predictions$Polygon)))

              #Loop through all Polygons
              for(i in unique(df_Polygons_NDVI_GAM_predictions$Polygon)){

                #i=unique(df_Polygons_NDVI_GAM_predictions$Polygon)[1]

                #Create an index variable for parameter i
                i_index <- which( unique(df_Polygons_NDVI_GAM_predictions$Polygon) == i)

                #select datasets for current NDVI_threshold and Polygon:
                df_Polygon_NDVI_GAM <- df_Polygons_NDVI_GAM[df_Polygons_NDVI_GAM$Polygon==i,]
                df_Polygon_NDVI_GAM_predictions <- df_Polygons_NDVI_GAM_predictions[df_Polygons_NDVI_GAM_predictions$Polygon==i,]

                #create plot for current Polygon and store it in list_plots_ndvi:
                list_plots_ndvi[[i_index]] <- ggplot()+
                  geom_point(data=df_Polygon_NDVI_GAM[df_Polygon_NDVI_GAM$outliers==FALSE,], aes(x=doy, y=NDVI))+
                  geom_point(data=df_Polygon_NDVI_GAM[df_Polygon_NDVI_GAM$outliers==TRUE,], aes(x=doy, y=NDVI), col="black", pch=16, alpha=0.2)+
                  geom_line(data=df_Polygon_NDVI_GAM_predictions, aes(x=doy, y=NDVI_gam_predict), col="#08a31a", lwd=1.25)+
                  geom_vline(xintercept = 150, colour="grey", lty=2)+
                  xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                  ylab("NDVI") +
                  ggtitle(paste0("Year: ", year_ID, ", Polygon: ", i))+
                  theme_tom()

                }

              #Plot NDVI and model predictions in a separate plot per Polygon
              # plots_ndvi <- list(list_plots_ndvi)
              # plots_per_page = 25
              # plots_ndvi <- lapply(plots_ndvi, function(x){split(x, ceiling(seq_along(plots_ndvi[[1]])/plots_per_page))})
              # plots_ndvi <- unname(unlist(plots_ndvi, recursive = F))
              # pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Plot_NDVI.pdf"), width=20, height=16, onefile = TRUE)
              # for (k in seq(length(plots_ndvi))) { do.call("grid.arrange", plots_ndvi[[k]]) }
              # dev.off()

      }

########################################################################################################################################################################################

      #(E): NDMI - Fit a Generalized Additive Model (GAM) through the NDMI data within each polygon
      if("NDMI" %in% method_output){

          #(E.0) Print message
            cat("\n")
            print("--------------------------------------------------------------------------------------------------------------------------")
            print(paste0("FIT A GAM THROUGH THE AVERAGE NDMI DATA PER POLYGON"))
            print("--------------------------------------------------------------------------------------------------------------------------")

          #(E.1) Create an empty dataframe
            df_Polygons_NDMI_GAM <- data.frame(NDMI=numeric(),
                                               Date=factor(),
                                               doy=numeric(),
                                               Polygon=factor(),
                                               outliers=logical())

            df_Polygons_NDMI_GAM_predictions <- data.frame(Polygon=character(),
                                                           doy=numeric(),
                                                           NDMI_gam_predict=numeric(),
                                                           stringsAsFactors=FALSE)

          #(E.2) Loop through all Polygons and fit a separate GAM (with sequential outlier removal) to the Polygon-specific NDMI data
            for(i in unique(df_Polygons_BandValues$Polygon)){

              #For debugging
              #i=unique(df_Polygons_BandValues$Polygon)[1]

              #Select Polygon-specific subset of data:
              df_Polygon_NDMI_GAM_new <- df_Polygons_BandValues[df_Polygons_BandValues$Polygon==i &
                                                                !is.na(df_Polygons_BandValues$NDMI),
                                                                c("NDMI", "Date", "doy", "Polygon")]

              #If there are at least two datapoints, then continue to fit the GAM
              if(nrow(df_Polygon_NDMI_GAM_new) > 1){

                #Fit a gam through the Polygon-specific NDMI ~ doy data and emply sequential outlier removal
                df_Polygon_NDMI_GAM_new <- f_gam_SeqRemOutliers(data=df_Polygon_NDMI_GAM_new, y="NDMI", x="doy", outlier_removal=outlier_removal,
                                                                 outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                                 default_k=gam_k_outlier)

                #Sort df_Polygon_NDMI_GAM_new by doy:
                df_Polygon_NDMI_GAM_new <- df_Polygon_NDMI_GAM_new[order(df_Polygon_NDMI_GAM_new$doy),]

                #Bind the Polygon-specific dataframe with outlier-filtered GAM estimates to the dataframe containing all dataframes from previous iterations
                 df_Polygons_NDMI_GAM <- rbind(df_Polygons_NDMI_GAM, df_Polygon_NDMI_GAM_new)

                #Create more detailed predictions (not only at the doy present in the datframe) at a 1-day interval to plot more smooth curves

                  #Refit GAM through data
                   index <- which(df_Polygon_NDMI_GAM_new$outliers==FALSE)
                   mod_gam <- with(df_Polygon_NDMI_GAM_new[index,], mgcv::gam(NDMI ~ s(doy, k=min(gam_k, length(index)-1)), method="REML"))

                  #Use gam to make predictions on a more detailed (1-day) day of year interval
                   df_Polygon_NDMI_GAM_predictions_new <- data.frame(Polygon=i, doy=seq(min(df_Polygon_NDMI_GAM_new$doy), max(df_Polygon_NDMI_GAM_new$doy), 0.01))
                   df_Polygon_NDMI_GAM_predictions_new$NDMI_gam_predict <- stats::predict(mod_gam, newdata=df_Polygon_NDMI_GAM_predictions_new, type="response")
                   df_Polygon_NDMI_GAM_predictions_new <- df_Polygon_NDMI_GAM_predictions_new[order(df_Polygon_NDMI_GAM_predictions_new$doy),]

                  #Add predictions to df_Polygons_NDMI_GAM_predictions dataframe:
                   df_Polygons_NDMI_GAM_predictions <- rbind(df_Polygons_NDMI_GAM_predictions, df_Polygon_NDMI_GAM_predictions_new)

                }

              #If there are less than two datapoints, then append empty dataframes
              if(nrow(df_Polygon_NDMI_GAM_new) < 2){

                #Bind the empty dataframe to the dataframe containing all dataframes from previous iterations
                if(nrow(df_Polygon_NDMI_GAM_new) < 1){df_Polygon_NDMI_GAM_new <- data.frame(NDMI=NA, Date=NA, doy=NA, Polygon=i, outliers=NA)}
                if(nrow(df_Polygon_NDMI_GAM_new) > 0){df_Polygon_NDMI_GAM_new$outliers <- NA}
                df_Polygons_NDMI_GAM <- rbind(df_Polygons_NDMI_GAM, df_Polygon_NDMI_GAM_new)

                #Add empty GAM predictions to df_Polygons_NDMI_GAM_predictions dataframe:
                df_Polygon_NDMI_GAM_predictions_new <- df_Polygon_NDMI_GAM_new[,c("Polygon", "doy")]
                df_Polygon_NDMI_GAM_predictions_new$NDMI_gam_predict <- NA
                df_Polygons_NDMI_GAM_predictions <- rbind(df_Polygons_NDMI_GAM_predictions, df_Polygon_NDMI_GAM_predictions_new)

                }

              }

            #Change Polygon column to factor
            df_Polygons_NDMI_GAM$Polygon <- as.factor(as.character(df_Polygons_NDMI_GAM$Polygon))
            df_Polygons_NDMI_GAM_predictions$Polygon <- as.factor(as.character(df_Polygons_NDMI_GAM_predictions$Polygon))

            #Add LocationID and year
            df_Polygons_NDMI_GAM_predictions$location <- area_name
            df_Polygons_NDMI_GAM_predictions$year <- year_ID

            #Save dataframe with GAM fits for NDMI
            #write.csv(df_Polygons_NDMI_GAM, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Data_NDMI_GAM.csv"), row.names = FALSE)
            #write.csv(df_Polygons_NDMI_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_GAM_Predictions_NDMI.csv"), row.names = FALSE)

          #(E.3) Plot the raw NDMI datapoints and gam predictions for each Polygon:

             #Plot NDMI and model predictions for all Polygons in a single plot
              p_Polygon_NDMI = ggplot()+
                geom_point(data=df_Polygons_NDMI_GAM, aes(x=doy, y=NDMI, fill=Polygon, col=Polygon))+
                geom_line(data=df_Polygons_NDMI_GAM_predictions, aes(x=doy, y=NDMI_gam_predict, col=Polygon)) +
                xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                ylab("NDMI") +
                theme_tom()

          #(E.4): Create a separate plot with GAM predictions per polygon

              #Create an empty list to store plots
              list_plots_ndmi <- vector('list', length(unique(df_Polygons_NDMI_GAM_predictions$Polygon)))

              #Loop through all Polygons
              for(i in unique(df_Polygons_NDMI_GAM_predictions$Polygon)){

                #i=unique(df_Polygons_NDMI_GAM_predictions$Polygon)[1]

                #Create an index variable for parameter i
                i_index <- which( unique(df_Polygons_NDMI_GAM_predictions$Polygon) == i)

                #select datasets for current NDMI_threshold and Polygon:
                df_Polygon_NDMI_GAM <- df_Polygons_NDMI_GAM[df_Polygons_NDMI_GAM$Polygon==i,]
                df_Polygon_NDMI_GAM_predictions <- df_Polygons_NDMI_GAM_predictions[df_Polygons_NDMI_GAM_predictions$Polygon==i,]

                #create plot for current Polygon and store it in list_plots_ndmi:
                list_plots_ndmi[[i_index]] <- ggplot()+
                  geom_point(data=df_Polygon_NDMI_GAM[df_Polygon_NDMI_GAM$outliers==FALSE,], aes(x=doy, y=NDMI))+
                  geom_point(data=df_Polygon_NDMI_GAM[df_Polygon_NDMI_GAM$outliers==TRUE,], aes(x=doy, y=NDMI), col="black", pch=16, alpha=0.2)+
                  geom_line(data=df_Polygon_NDMI_GAM_predictions, aes(x=doy, y=NDMI_gam_predict), col="#16acde", lwd=1.25)+
                  geom_vline(xintercept = 150, colour="grey", lty=2)+
                  xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                  ylab("NDMI") +
                  ggtitle(paste0("Year: ", year_ID, ", Polygon: ", i))+
                  theme_tom()

              }

              #Plot NDMI and model predictions in a separate plot per Polygon
              # plots_ndmi <- list(list_plots_ndmi)
              # plots_per_page = 25
              # plots_ndmi <- lapply(plots_ndmi, function(x){split(x, ceiling(seq_along(plots_ndmi[[1]])/plots_per_page))})
              # plots_ndmi <- unname(unlist(plots_ndmi, recursive = F))
              # pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Plot_NDMI.pdf"), width=20, height=16, onefile = TRUE)
              # for (k in seq(length(plots_ndmi))) { do.call("grid.arrange", plots_ndmi[[k]]) }
              # dev.off()
      }

########################################################################################################################################################################################

      #(F) Plot NDSI, NDVI and NDMI together in a plot per Polygon

        # #Plot:
        #  p_Polygon_BANDS_grid <- ggplot()+
        #    geom_point(data=df_Polygons_NDSI_GAM[df_Polygons_NDSI_GAM$outliers==FALSE,], aes(x=doy, y=NDSI))+
        #    geom_point(data=df_Polygons_NDSI_GAM[df_Polygons_NDSI_GAM$outliers==TRUE,], aes(x=doy, y=NDSI), col="black", pch=16, alpha=0.2)+
        #    geom_point(data=df_Polygons_NDVI_GAM[df_Polygons_NDVI_GAM$outliers==FALSE,], aes(x=doy, y=NDVI), col="black")+
        #    geom_point(data=df_Polygons_NDVI_GAM[df_Polygons_NDVI_GAM$outliers==TRUE,], aes(x=doy, y=NDVI), col="black", pch=16, alpha=0.2)+
        #    geom_point(data=df_Polygons_NDMI_GAM[df_Polygons_NDMI_GAM$outliers==FALSE,], aes(x=doy, y=NDMI), col="black")+
        #    geom_point(data=df_Polygons_NDMI_GAM[df_Polygons_NDMI_GAM$outliers==TRUE,], aes(x=doy, y=NDMI), col="black", pch=16, alpha=0.2)+
        #    geom_line(data=df_Polygons_NDSI_GAM_predictions, aes(x=doy, y=NDSI_gam_predict), col="#1620de", lwd=1.25)+
        #    geom_line(data=df_Polygons_NDVI_GAM_predictions, aes(x=doy, y=NDVI_gam_predict), col="#08a31a", lwd=1.25)+
        #    geom_line(data=df_Polygons_NDMI_GAM_predictions, aes(x=doy, y=NDMI_gam_predict), col="#16acde", lwd=1.25)+
        #    geom_vline(xintercept = 150, colour="grey", lty=2)+
        #    facet_wrap(~Polygon, ncol=ceiling(length(unique(df_Polygons_NDMI_GAM$Polygon))^0.5))+
        #    xlab("Day of year") +
        #    ylab("Normalized Difference Band Index") +
        #    theme_tom()
        #
        #    ggsave(plot=p_Polygon_BANDS_grid, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", data_ID, "_Res", resolution, "_Polygons_Plot_AllBands.pdf"), width=14, height = 12)

}


########################################################################################################################################################################################

  #Store dataframes and plots in a separate list element per year

    #(1): Pixelcounts
    if(pixel_counts==TRUE){
      list_df_pixelcount[[which(year_ID==year_ID_vector)]] <- df_pixelcount
      list_p_pixelcounts[[which(year_ID==year_ID_vector)]] <- p_pixelcounts
      }

    #(2): Raw data
    if("avg_NDSI" %in% method){
      list_df_Polygons_BandValues[[which(year_ID==year_ID_vector)]] <- df_Polygons_BandValues
      }
    if("snowfraction" %in% method){
      list_df_Polygons_SnowFraction[[which(year_ID==year_ID_vector)]] <- df_Polygons_SnowFraction
      }

    #(3): Snowfraction
    if("snowfraction" %in% method){
      list_df_Polygons_SnowFraction_GAM_predictions[[which(year_ID==year_ID_vector)]] <- df_Polygons_SnowFraction_GAM_predictions
      list_df_Polygon_Snowfraction[[which(year_ID==year_ID_vector)]] <- df_Polygon_Snowfraction
      list_p_snowfraction[[which(year_ID==year_ID_vector)]] <- list_plots_snowfraction
      }

    #(4): FSC Gascoin2020
    if(("avg_NDSI" %in% method) & ("FSC_Gascoin2020" %in% method_output)){
      list_df_Polygons_FSC_Gascoin2020_GAM_predictions[[which(year_ID==year_ID_vector)]] <- df_Polygons_FSC_Gascoin2020_GAM_predictions
      list_df_Polygon_FSC_Gascoin2020[[which(year_ID==year_ID_vector)]] <- df_Polygon_FSC_Gascoin2020
      list_p_fsc_gascoin2020[[which(year_ID==year_ID_vector)]] <- list(list_plots_fsc_gascoin2020)
      }

    #(5): FSC Aalstad2020
    if(("avg_NDSI" %in% method) & ("FSC_Aalstad2020" %in% method_output)){
      list_df_Polygons_FSC_Aalstad2020_GAM_predictions[[which(year_ID==year_ID_vector)]] <- df_Polygons_FSC_Aalstad2020_GAM_predictions
      list_df_Polygon_FSC_Aalstad2020[[which(year_ID==year_ID_vector)]] <- df_Polygon_FSC_Aalstad2020
      list_p_fsc_aalstad2020[[which(year_ID==year_ID_vector)]] <- list(list_plots_fsc_aalstad2020)
      }

    #(6): NDSI
    if(("avg_NDSI" %in% method) & ("NDSI" %in% method_output)){
      list_df_Polygons_NDSI_GAM_predictions[[which(year_ID==year_ID_vector)]] <- df_Polygons_NDSI_GAM_predictions
      list_df_Polygon_NDSI[[which(year_ID==year_ID_vector)]] <- df_Polygon_NDSI
      list_p_ndsi[[which(year_ID==year_ID_vector)]] <- list(list_plots_ndsi)
      }

    #(7): NDVI
    if(("avg_NDSI" %in% method) & ("NDVI" %in% method_output)){
      list_df_Polygons_NDVI_GAM_predictions[[which(year_ID==year_ID_vector)]] <- df_Polygons_NDVI_GAM_predictions
      list_p_ndvi[[which(year_ID==year_ID_vector)]] <- list(list_plots_ndvi)
      }

    #(8): NDMI
    if(("avg_NDSI" %in% method) & ("NDMI" %in% method_output)){
      list_df_Polygons_NDMI_GAM_predictions[[which(year_ID==year_ID_vector)]] <- df_Polygons_NDMI_GAM_predictions
      list_p_ndmi[[which(year_ID==year_ID_vector)]] <- list(list_plots_ndmi)
      }

}

  #Store the results for all years for the current shapefile

    #(1): Pixel counts:
    if(pixel_counts==TRUE){

       #Dataframe
       df_pixelcount <- do.call("rbind", list_df_pixelcount)
       write.csv(df_pixelcount, file=paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Data_Pixel_Counts_polygon.csv"), quote=FALSE, row.names=FALSE)

       #Plot:
       list_p_pixelcounts <- lapply(list_p_pixelcounts, function(x) {do.call(list, list(x))})
       pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Plot_Pixel_Counts_polygon.pdf"), width=12, height=8, onefile = TRUE)
       for(k in seq(length(list_p_pixelcounts))) {do.call("grid.arrange", c(list_p_pixelcounts[[k]], ncol=1, nrow=1))}
       dev.off()

       }

    #(2): Raw data

       #Timeseries of Mean bandvalues:
       if("avg_NDSI" %in% method){
        df_Polygons_BandValues <- do.call("rbind", list_df_Polygons_BandValues)
        write.csv(df_Polygons_BandValues, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Data_MeanBandValues.csv"), row.names = FALSE)
        }

       #Timeseries of fraction of snowcover
       if("snowfraction" %in% method){
        df_Polygons_SnowFraction <- do.call("rbind", list_df_Polygons_SnowFraction)
        write.csv(df_Polygons_SnowFraction, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Data_SnowFraction.csv"), row.names = FALSE)
        }

    #(3): Snowfraction
    if("snowfraction" %in% method){

       #Predictions of GAM fit
       df_Polygons_SnowFraction_GAM_predictions <- do.call("rbind", list_df_Polygons_SnowFraction_GAM_predictions)
       write.csv(df_Polygons_SnowFraction_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_GAM_Predictions_SnowFraction.csv"), row.names = FALSE)

       #Date of snowmelt
       df_Polygon_Snowfraction <- do.call("rbind", list_df_Polygon_Snowfraction)
       write.csv(df_Polygon_Snowfraction, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Snowmelt_Snowfraction.csv"), row.names = FALSE)

       #Plot (separate page per year, separate plot per NDSI-threshold and per year):

         #Load plots
         plots_snowfraction <- list_p_snowfraction
         plots_per_page = 25

         #Collapse list within each year
         plots_snowfraction <- lapply(plots_snowfraction, function(inner_list) {do.call(c, inner_list)})

         #Split list of plots to only contain maximum 25 plots per page
         plots_snowfraction <- lapply(plots_snowfraction, function(x){split(x, ceiling(seq_along(plots_snowfraction[[1]])/plots_per_page))})
         plots_snowfraction <- unname(unlist(plots_snowfraction, recursive = F))

         #Determine number of columns based on number of plots per page
         max_plots <- max(unlist(lapply(plots_snowfraction, function(x){length(x)})))
         if(max_plots<2){var_ncol<-1}
         if(max_plots>1 & max_plots<5){var_ncol<-2}
         if(max_plots>4 & max_plots<10){var_ncol<-3}
         if(max_plots>9 & max_plots<17){var_ncol<-4}
         if(max_plots>16){var_ncol<-5}

         #Create plot
         pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Plot_Snowfraction.pdf"), width=20, height=16, onefile = TRUE)
         for(k in seq(length(plots_snowfraction))) {do.call("grid.arrange", c(plots_snowfraction[[k]], ncol=var_ncol, nrow=var_ncol))}
         dev.off()

        }

    #(4): FSC Gascoin2020
    if(("avg_NDSI" %in% method) & ("FSC_Gascoin2020" %in% method_output)){

       #Predictions of GAM fit
       df_Polygons_FSC_Gascoin2020_GAM_predictions <- do.call("rbind", list_df_Polygons_FSC_Gascoin2020_GAM_predictions)
       write.csv(df_Polygons_FSC_Gascoin2020_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_GAM_Predictions_FSC_Gascoin2020.csv"), row.names = FALSE)

       #Date of snowmelt
       df_Polygon_FSC_Gascoin2020 <- do.call("rbind", list_df_Polygon_FSC_Gascoin2020)
       write.csv(df_Polygon_FSC_Gascoin2020, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Snowmelt_FSC_Gascoin2020.csv"), row.names = FALSE)

       #Plot (separate page per year, separate plot per polygon:

         #Load plots
         plots_fsc_gascoin2020 <- list_p_fsc_gascoin2020
         plots_per_page = 25

         #Collapse list within each year
         plots_fsc_gascoin2020 <- lapply(plots_fsc_gascoin2020, function(inner_list) {do.call(c, inner_list)})

         #Split list of plots to only contain maximum 25 plots per page
         plots_fsc_gascoin2020 <- lapply(plots_fsc_gascoin2020, function(x){split(x, ceiling(seq_along(plots_fsc_gascoin2020[[1]])/plots_per_page))})
         plots_fsc_gascoin2020 <- unname(unlist(plots_fsc_gascoin2020, recursive = F))

         #Determine number of columns based on number of plots per page
         max_plots <- max(unlist(lapply(plots_fsc_gascoin2020, function(x){length(x)})))
         if(max_plots<2){var_ncol<-1}
         if(max_plots>1 & max_plots<5){var_ncol<-2}
         if(max_plots>4 & max_plots<10){var_ncol<-3}
         if(max_plots>9 & max_plots<17){var_ncol<-4}
         if(max_plots>16){var_ncol<-5}

         #Create plot
         pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Plot_FSC_Gascoin2020.pdf"), width=20, height=16, onefile = TRUE)
         for(k in seq(length(plots_fsc_gascoin2020))) {do.call("grid.arrange", c(plots_fsc_gascoin2020[[k]], ncol=var_ncol, nrow=var_ncol))}
         dev.off()

        }

    #(5): FSC Aalstad2020
    if(("avg_NDSI" %in% method) & ("FSC_Aalstad2020" %in% method_output)){

       #Predictions of GAM fit
       df_Polygons_FSC_Aalstad2020_GAM_predictions <- do.call("rbind", list_df_Polygons_FSC_Aalstad2020_GAM_predictions)
       write.csv(df_Polygons_FSC_Aalstad2020_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_GAM_Predictions_FSC_Aalstad2020.csv"), row.names = FALSE)

       #Date of snowmelt
       df_Polygon_FSC_Aalstad2020 <- do.call("rbind", list_df_Polygon_FSC_Aalstad2020)
       write.csv(df_Polygon_FSC_Aalstad2020, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Snowmelt_FSC_Aalstad2020.csv"), row.names = FALSE)

       #Plot (separate page per year, separate plot per polygon:

         #Load plots
         plots_fsc_aalstad2020 <- list_p_fsc_aalstad2020
         plots_per_page = 25

         #Collapse list within each year
         plots_fsc_aalstad2020 <- lapply(plots_fsc_aalstad2020, function(inner_list) {do.call(c, inner_list)})

         #Split list of plots to only contain maximum 25 plots per page
         plots_fsc_aalstad2020 <- lapply(plots_fsc_aalstad2020, function(x){split(x, ceiling(seq_along(plots_fsc_aalstad2020[[1]])/plots_per_page))})
         plots_fsc_aalstad2020 <- unname(unlist(plots_fsc_aalstad2020, recursive = F))

         #Determine number of columns based on number of plots per page
         max_plots <- max(unlist(lapply(plots_fsc_aalstad2020, function(x){length(x)})))
         if(max_plots<2){var_ncol<-1}
         if(max_plots>1 & max_plots<5){var_ncol<-2}
         if(max_plots>4 & max_plots<10){var_ncol<-3}
         if(max_plots>9 & max_plots<17){var_ncol<-4}
         if(max_plots>16){var_ncol<-5}

         #Create plot
         pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Plot_FSC_Aalstad2020.pdf"), width=20, height=16, onefile = TRUE)
         for(k in seq(length(plots_fsc_aalstad2020))) {do.call("grid.arrange", c(plots_fsc_aalstad2020[[k]], ncol=var_ncol, nrow=var_ncol))}
         dev.off()

        }

    #(6): NDSI
    if(("avg_NDSI" %in% method) & ("NDSI" %in% method_output)){

       #Predictions of GAM fit
       df_Polygons_NDSI_GAM_predictions <- do.call("rbind", list_df_Polygons_NDSI_GAM_predictions)
       write.csv(df_Polygons_NDSI_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_GAM_Predictions_NDSI.csv"), row.names = FALSE)

       #Date of snowmelt
       df_Polygon_NDSI <- do.call("rbind", list_df_Polygon_NDSI)
       write.csv(df_Polygon_NDSI, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Snowmelt_NDSI.csv"), row.names = FALSE)

       #Plot (separate page per year, separate plot per polygon:

         #Load plots
         plots_ndsi <- list_p_ndsi
         plots_per_page = 25

         #Collapse list within each year
         plots_ndsi <- lapply(plots_ndsi, function(inner_list) {do.call(c, inner_list)})

         #Split list of plots to only contain maximum 25 plots per page
         plots_ndsi <- lapply(plots_ndsi, function(x){split(x, ceiling(seq_along(plots_ndsi[[1]])/plots_per_page))})
         plots_ndsi <- unname(unlist(plots_ndsi, recursive = F))

         #Determine number of columns based on number of plots per page
         max_plots <- max(unlist(lapply(plots_ndsi, function(x){length(x)})))
         if(max_plots<2){var_ncol<-1}
         if(max_plots>1 & max_plots<5){var_ncol<-2}
         if(max_plots>4 & max_plots<10){var_ncol<-3}
         if(max_plots>9 & max_plots<17){var_ncol<-4}
         if(max_plots>16){var_ncol<-5}

         #Create plot
         pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Plot_NDSI.pdf"), width=20, height=16, onefile = TRUE)
         for(k in seq(length(plots_ndsi))) {do.call("grid.arrange", c(plots_ndsi[[k]], ncol=var_ncol, nrow=var_ncol))}
         dev.off()

        }

    #(7): NDVI
    if(("avg_NDSI" %in% method) & ("NDVI" %in% method_output)){

      #Predictions of GAM fit
      df_Polygons_NDVI_GAM_predictions <- do.call("rbind", list_df_Polygons_NDVI_GAM_predictions)
      write.csv(df_Polygons_NDVI_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_GAM_Predictions_NDVI.csv"), row.names = FALSE)

      #Plot (separate page per year, separate plot per polygon:

        #Load plots
        plots_ndvi <- list_p_ndvi
        plots_per_page = 25

        #Collapse list within each year
        plots_ndvi <- lapply(plots_ndvi, function(inner_list) {do.call(c, inner_list)})

        #Split list of plots to only contain maximum 25 plots per page
        plots_ndvi <- lapply(plots_ndvi, function(x){split(x, ceiling(seq_along(plots_ndvi[[1]])/plots_per_page))})
        plots_ndvi <- unname(unlist(plots_ndvi, recursive = F))

        #Determine number of columns based on number of plots per page
        max_plots <- max(unlist(lapply(plots_ndvi, function(x){length(x)})))
        if(max_plots<2){var_ncol<-1}
        if(max_plots>1 & max_plots<5){var_ncol<-2}
        if(max_plots>4 & max_plots<10){var_ncol<-3}
        if(max_plots>9 & max_plots<17){var_ncol<-4}
        if(max_plots>16){var_ncol<-5}

        #Create plot
        pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Plot_NDVI.pdf"), width=20, height=16, onefile = TRUE)
        for(k in seq(length(plots_ndvi))) {do.call("grid.arrange", c(plots_ndvi[[k]], ncol=var_ncol, nrow=var_ncol))}
        dev.off()

        }

    #(8): NDMI
    if(("avg_NDSI" %in% method) & ("NDMI" %in% method_output)){

      #Predictions of GAM fit
      df_Polygons_NDMI_GAM_predictions <- do.call("rbind", list_df_Polygons_NDMI_GAM_predictions)
      write.csv(df_Polygons_NDMI_GAM_predictions, paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_GAM_Predictions_NDMI.csv"), row.names = FALSE)

      #Plot (separate page per year, separate plot per polygon:

        #Load plots
        plots_ndmi <- list_p_ndmi
        plots_per_page = 25

        #Collapse list within each year
        plots_ndmi <- lapply(plots_ndmi, function(inner_list) {do.call(c, inner_list)})

        #Split list of plots to only contain maximum 25 plots per page
        plots_ndmi <- lapply(plots_ndmi, function(x){split(x, ceiling(seq_along(plots_ndmi[[1]])/plots_per_page))})
        plots_ndmi <- unname(unlist(plots_ndmi, recursive = F))

        #Determine number of columns based on number of plots per page
        max_plots <- max(unlist(lapply(plots_ndmi, function(x){length(x)})))
        if(max_plots<2){var_ncol<-1}
        if(max_plots>1 & max_plots<5){var_ncol<-2}
        if(max_plots>4 & max_plots<10){var_ncol<-3}
        if(max_plots>9 & max_plots<17){var_ncol<-4}
        if(max_plots>16){var_ncol<-5}

        #Create plot
        pdf(paste0(here(), "/Output/MODIS/03_Polygons_Snowmelt/", timestamp, "_", area_name, "_MODIS_", MODIS_cloud_masking_algorithm, "_Res", resolution, "_Polygons_Plot_NDMI.pdf"), width=20, height=16, onefile = TRUE)
        for(k in seq(length(plots_ndmi))) {do.call("grid.arrange", c(plots_ndmi[[k]], ncol=var_ncol, nrow=var_ncol))}
        dev.off()

        }

 }

##################################################################################################################################

#The End    
     
 #Print concluding remarks
 {cat("\n") ; cat("\n")
  print("--------------------------------------------------------------------------------------------------------------------------")
  print(paste0("THE ANALYSIS HAS COMPLETED"))
  cat("\n")
  print(paste0("-Generated data and plots can be found at ", here(), "/Output/MODIS/03_Polygons_Snowmelt/"))
  print("--------------------------------------------------------------------------------------------------------------------------")
  cat("\n")}  
     
##################################################################################################################################

