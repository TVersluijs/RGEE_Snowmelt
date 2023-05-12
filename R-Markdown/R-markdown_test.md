RGEE MODIS - Calculate the date of snowmelt for all pixels within an
area of interest (shapefile)
================
Tom Versluijs
2023-05-12

Extract MODIS satellite data and Calculate the date of snowmelt for
every 500mx500m pixel in an area of interest (shapefile).

## Copyright Tom Versluijs 2023-02-08. Do not use this code without permission. Contact information: <tom.versluijs@gmail.com>

# I. Setup Workspace

------------------------------------------------------------------------

### (1): Clear workspace and load packages

``` r
  #Clear workspace and set python environment
   rm(list=ls())
   library(here)
   rgee_environment_dir <- readRDS(paste0(here::here(), "/Input/rgee_environment_dir.rds"))
   reticulate::use_python(rgee_environment_dir, required=T)
```

    ## Warning: The request to
    ## `use_python("C:\Users\tomve\miniconda3\envs\rgee_py/python.exe")` will be
    ## ignored because the environment variable RETICULATE_PYTHON is set to
    ## "C:\Users\tomve\miniconda3\envs\rgee_py"

``` r
   reticulate::py_config()
```

    ## python:         C:/Users/tomve/miniconda3/envs/rgee_py/python.exe
    ## libpython:      C:/Users/tomve/miniconda3/envs/rgee_py/python311.dll
    ## pythonhome:     C:/Users/tomve/miniconda3/envs/rgee_py
    ## version:        3.11.3 | packaged by Anaconda, Inc. | (main, Apr 19 2023, 23:46:34) [MSC v.1916 64 bit (AMD64)]
    ## Architecture:   64bit
    ## numpy:          C:/Users/tomve/miniconda3/envs/rgee_py/Lib/site-packages/numpy
    ## numpy_version:  1.24.3
    ## 
    ## NOTE: Python version was forced by RETICULATE_PYTHON

``` r
  #Load packages
   library(pacman)
   p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr, foreach, parallel, doSNOW, gridExtra, webshot)
   webshot::install_phantomjs()
```

    ## Warning in utils::download.file(url, method = method, ...): the 'wininet'
    ## method is deprecated for http:// and https:// URLs

### (2): Define ggplot2 plotting theme

``` r
    theme_tom <- function(){
       theme_classic() %+replace%
       theme(axis.title = element_text(size=18),
           axis.text = element_text(size=16),
           legend.position = "none",
           strip.background = element_rect(fill = "white", colour = "black", size = rel(2)),
           complete = TRUE)}
```

### (3): Load auxiliary functions

``` r
   source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "MODIS_AuxiliaryFunctions")
   sapply(source_files, source, chdir = TRUE) ; rm(source_files)
```

### (4): Initialize earth engine and google drive

``` r
   rgee::ee_Initialize(user="tom.versluijs@gmail.com", drive = T)
```

    ## ── rgee 1.1.6.9999 ────────────────────────────────── earthengine-api 0.1.353 ── 
    ##  ✔ user: tom.versluijs@gmail.com 
    ##  ✔ Google Drive credentials: ✔ Google Drive credentials:  FOUND
    ##  ✔ Initializing Google Earth Engine: ✔ Initializing Google Earth Engine:  DONE!
    ##  ✔ Earth Engine account: users/escape 
    ## ────────────────────────────────────────────────────────────────────────────────

------------------------------------------------------------------------

# II. Specify parameters of interest

------------------------------------------------------------------------

### (5): Specify parameters used in analysis

``` r
   #(a): MODIS satellite

     #Specify resolution of images in meters
     resolution=500 # native resolution of satellite image (MODIS=500m)

   #(b) Area of interest

     #Specify name of study area
     area_name <- "ZAC"

     #Approximate central point of study area
     coordinates_point <- c(-20.49, 74.49)

     #Bounding box for study area (general area of interest)
     aoi <- ee$Geometry$Polygon(list(c(-20.34, 74.542),
                                     c(-20.7,  74.542),
                                     c(-20.7,  74.442),
                                     c(-20.34, 74.442)))

     #Name of Shapefile located in Input folder (specific area of interest)
     shapefile <- "ZAC_Outline_EPSG32627.shp"

     #Coordinate reference system
     #UTM: use "epsg:326XX" for northern hemisphere, "epsg:327XX" for southern hemisphere, where XX is the specific zone
     crs <- "epsg:32627"

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
```

------------------------------------------------------------------------

# III. Create a unique asset folder

------------------------------------------------------------------------

### (6) Create a unique dataID and asset folder for storing generated datafiles

``` r
    #Create a unique data_ID
     data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
     data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)

    #Create a unique Asset folder (delete this folder if already present)
     path_asset <- paste0(ee_get_assethome(), "/", data_ID)
     #ee_manage_assetlist(path_asset)
     ee_manage_delete(path_asset)
```

    ## EE object deleted: users/escape/ZAC22_MODIS_PGE11

``` r
     ee_manage_create(path_asset=path_asset, asset_type="Folder")
```

    ## GEE asset: users/escape/ZAC22_MODIS_PGE11 created

``` r
     #ee_manage_assetlist()
```

------------------------------------------------------------------------

# IV. Read and display the unfiltered data

------------------------------------------------------------------------

### (7) Extract MODIS satellite Surface Reflectance images for specified daterange and general area of interest

Note that specifying a region of interest using a polygon does not
function properly for the MODIS data. Instead, providing an initial
point and then afterwards clipping the image to a desired area of
interest works better.

``` r
    #Specify starting and ending date
     start_date_doy <- as.numeric(strftime(start_date, format = "%j"))
     end_date_doy <- as.numeric(strftime(end_date, format = "%j"))

    #Extract Sentinel-2 Surface Reflectance satellite data
     MODIS_col <- ee$ImageCollection('MODIS/006/MOD09GA')
     point <- ee$Geometry$Point(coordinates_point[1], coordinates_point[2])
     MODIS_col <- MODIS_col$
       filterBounds(point)$
       filterDate(start_date, end_date)
```

### (8) Read study area shapefile as a feature collection and clip image collection to area of interest

Make sure that the shapefile is in EPSG:326XX format (UTM). In addition
the ‘id’ column of the shapefile should not be empty (i.e. add 0 instead
of NULL or NA in QGIS layer attribute table)

``` r
    #Read aoi_Shapefile shapefile in root folder
     root_fldr <- here()
     aoi_Shapefile <- st_read(paste0(root_fldr, "/Input/Shapefiles/", shapefile), quiet=T)
     aoi_Shapefile <- st_transform(aoi_Shapefile, crs=crs)
     aoi_Shapefile <- sf_as_ee(aoi_Shapefile)
```

    ## Registered S3 method overwritten by 'geojsonsf':
    ##   method        from   
    ##   print.geojson geojson

``` r
     aoi_Shapefile <- ee$FeatureCollection(aoi_Shapefile)

    #Calculate the size of the study area in km2:
     img <- ee$Image$pixelArea()$divide(1000000)
     area2 <- img$reduceRegion(
       reducer= ee$Reducer$sum(),
       geometry= aoi_Shapefile,
       crs= crs,
       scale= 10,
       maxPixels= 1E13)
     paste0('Size of study area calculated using the pixel area method: ', round(ee$Number(area2$get('area'))$getInfo(),3), ' km2')
```

    ## [1] "Size of study area calculated using the pixel area method: 50.008 km2"

``` r
    #Clip all images to area of interest (aoi). Clipping to shapefile will occur at step 16B of this script
     MODIS_col <- MODIS_col$map(function(img){return(img$clip(aoi))})

    #Plot a single clipped image:

      #Select a single image for initial plot, clipped to aoi_Shapefile:
       image <- MODIS_col$filterDate(paste0(year_ID, "-06-10"), paste0(year_ID, "-08-15"))$first()$clipToCollection(aoi_Shapefile)

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
       Map$setCenter(coordinates_point[1], coordinates_point[2], 10)
       Map$addLayer(image,list(bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')+
       Map$addLayer(ndsi_MODIS,list(min=-1, max=1.5, palette=c('black', '0dffff', '0524ff', 'ffffff')), 'NDSI')+
       Map$addLayer(ndvi_MODIS,list(min=-1, max=1, palette=c('#FF0000','#00FF00')), 'NDVI')+
       Map$addLayer(ndwi_MODIS,list(min=0, max=1, palette=c('000000', '0dffff', '0524ff', 'ffffff')), 'NDWI')
```

![](R-markdown_test_files/figure-gfm/read%20shapefile%20and%20plot%20initial%20data-1.png)<!-- -->

Built with R-version 4.3.0
