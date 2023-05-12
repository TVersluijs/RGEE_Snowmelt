---
title: "RGEE MODIS - Calculate the date of snowmelt for all pixels within an area of interest (shapefile)"
author: "Tom Versluijs"
date: "2023-05-12"
output: 
  html_document:
    keep_md: TRUE
---



Extract MODIS satellite data and Calculate the date of snowmelt for every 500mx500m pixel in an area of interest (shapefile).

## Copyright Tom Versluijs 2023-02-08. Do not use this code without permission. Contact information: [tom.versluijs\@gmail.com](mailto:tom.versluijs@gmail.com){.email}

# I. Setup Workspace

------------------------------------------------------------------------

### (1): Clear workspace and load packages


```r
#Clear workspace and set python environment
 rm(list=ls())
 library(here)
 rgee_environment_dir <- readRDS(paste0(here::here(), "/Input/rgee_environment_dir.rds"))
 reticulate::use_python(rgee_environment_dir, required=T)
 reticulate::py_config()
```

```
## python:         C:/Users/tomve/miniconda3/envs/rgee_py/python.exe
## libpython:      C:/Users/tomve/miniconda3/envs/rgee_py/python311.dll
## pythonhome:     C:/Users/tomve/miniconda3/envs/rgee_py
## version:        3.11.3 | packaged by Anaconda, Inc. | (main, Apr 19 2023, 23:46:34) [MSC v.1916 64 bit (AMD64)]
## Architecture:   64bit
## numpy:          C:/Users/tomve/miniconda3/envs/rgee_py/Lib/site-packages/numpy
## numpy_version:  1.24.3
## 
## NOTE: Python version was forced by use_python function
```

```r
#Load packages
 library(pacman)
 p_load(sf, rgee, ggplot2, mgcv, googledrive, dplyr, foreach, parallel, doSNOW, gridExtra, webshot)
 #webshot::install_phantomjs(force=TRUE)
```

### (2): Define ggplot2 plotting theme


```r
theme_tom <- function(){
   theme_classic() %+replace%
   theme(axis.title = element_text(size=18),
       axis.text = element_text(size=16),
       legend.position = "none",
       strip.background = element_rect(fill = "white", colour = "black", size = rel(2)),
       complete = TRUE)}
```

### (3): Load auxiliary functions


```r
source_files <- list.files(path=paste0(here(), "/Input"), full.names=TRUE, recursive = TRUE, pattern = "MODIS_AuxiliaryFunctions")
sapply(source_files, source, chdir = TRUE) ; rm(source_files)
```

### (4): Initialize earth engine and google drive


```r
rgee::ee_Initialize(user="tom.versluijs@gmail.com", drive = T)
```

```
## ── rgee 1.1.6.9999 ────────────────────────────────── earthengine-api 0.1.353 ── 
##  ✔ user: tom.versluijs@gmail.com 
##  ✔ Google Drive credentials:
 ✔ Google Drive credentials:  FOUND
##  ✔ Initializing Google Earth Engine:
 ✔ Initializing Google Earth Engine:  DONE!
## 
 ✔ Earth Engine account: users/escape 
## ────────────────────────────────────────────────────────────────────────────────
```

------------------------------------------------------------------------

# II. Specify parameters of interest

------------------------------------------------------------------------

### (5): Specify parameters used in analysis


```r
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


```r
#Create a unique data_ID
 data_ID <- paste0(area_name, substr(year_ID,(nchar(year_ID)+1)-2,nchar(year_ID)), "_MODIS")
 data_ID <- paste0(data_ID, "_", MODIS_cloud_masking_algorithm)

#Create a unique Asset folder (delete this folder if already present)
 path_asset <- paste0(ee_get_assethome(), "/", data_ID)
 #tryCatch(ee_manage_assetlist(path_asset), error=function(error_message) {message("path_asset does not yet exist")})
 tryCatch(ee_manage_delete(path_asset), error=function(error_message){message("path_asset does not yet exist")})
```

```
## EE object deleted: users/escape/ZAC22_MODIS_PGE11
```

```r
 ee_manage_create(path_asset=path_asset, asset_type="Folder")
```

```
## GEE asset: users/escape/ZAC22_MODIS_PGE11 created
```

```r
 #ee_manage_assetlist()
```

------------------------------------------------------------------------

# IV. Read and display the unfiltered data

------------------------------------------------------------------------

### (7) Extract MODIS satellite Surface Reflectance images for specified daterange and general area of interest

Note that specifying a region of interest using a polygon does not function properly for the MODIS data. Instead, providing an initial point and then afterwards clipping the image to a desired area of interest works better.


```r
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

Make sure that the shapefile is in EPSG:326XX format (UTM). In addition the 'id' column of the shapefile should not be empty (i.e. add 0 instead of NULL or NA in QGIS layer attribute table)


```r
#Read aoi_Shapefile shapefile in root folder
 root_fldr <- here()
 aoi_Shapefile <- st_read(paste0(root_fldr, "/Input/Shapefiles/", shapefile), quiet=T)
 aoi_Shapefile <- st_transform(aoi_Shapefile, crs=crs)
 aoi_Shapefile <- sf_as_ee(aoi_Shapefile)
```

```
## Registered S3 method overwritten by 'geojsonsf':
##   method        from   
##   print.geojson geojson
```

```r
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

```
## [1] "Size of study area calculated using the pixel area method: 50.008 km2"
```

```r
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

```{=html}
<div class="leaflet html-widget html-fill-item-overflow-hidden html-fill-item" id="htmlwidget-521b2c698f0f0582fe12" style="width:672px;height:480px;"></div>
<script type="application/json" data-for="htmlwidget-521b2c698f0f0582fe12">{"x":{"options":{"minZoom":1,"maxZoom":24,"crs":{"crsClass":"L.CRS.EPSG3857","code":null,"proj4def":null,"projectedBounds":null,"options":{}},"preferCanvas":false,"bounceAtZoomLimits":false,"maxBounds":[[[-90,-370]],[[90,370]]]},"calls":[{"method":"addProviderTiles","args":["CartoDB.Positron","CartoDB.Positron","CartoDB.Positron",{"errorTileUrl":"","noWrap":false,"detectRetina":false,"pane":"tilePane","maxZoom":24}]},{"method":"addProviderTiles","args":["OpenStreetMap","OpenStreetMap","OpenStreetMap",{"errorTileUrl":"","noWrap":false,"detectRetina":false,"pane":"tilePane","maxZoom":24}]},{"method":"addProviderTiles","args":["CartoDB.DarkMatter","CartoDB.DarkMatter","CartoDB.DarkMatter",{"errorTileUrl":"","noWrap":false,"detectRetina":false,"pane":"tilePane","maxZoom":24}]},{"method":"addProviderTiles","args":["Esri.WorldImagery","Esri.WorldImagery","Esri.WorldImagery",{"errorTileUrl":"","noWrap":false,"detectRetina":false,"pane":"tilePane","maxZoom":24}]},{"method":"addProviderTiles","args":["OpenTopoMap","OpenTopoMap","OpenTopoMap",{"errorTileUrl":"","noWrap":false,"detectRetina":false,"pane":"tilePane","maxZoom":24}]},{"method":"addLayersControl","args":[["CartoDB.Positron","OpenStreetMap","CartoDB.DarkMatter","Esri.WorldImagery","OpenTopoMap"],[],{"collapsed":true,"autoZIndex":true,"position":"topleft"}]},{"method":"addScaleBar","args":[{"maxWidth":100,"metric":true,"imperial":true,"updateWhenIdle":true,"position":"bottomleft"}]},{"method":"addTiles","args":["https://earthengine.googleapis.com/v1alpha/projects/earthengine-legacy/maps/581e05a7d3f792b0f8ce13e8fb1a8732-bb1ddf938027073eedaf84b632e5f4e3/tiles/{z}/{x}/{y}","TRUE COLOR","TRUE COLOR",{"minZoom":0,"maxZoom":24,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false}]},{"method":"addLayersControl","args":[["CartoDB.Positron","OpenStreetMap","CartoDB.DarkMatter","Esri.WorldImagery","OpenTopoMap"],"TRUE COLOR",{"collapsed":true,"autoZIndex":true,"position":"topleft"}]},{"method":"hideGroup","args":[null]},{"method":"addTiles","args":["https://earthengine.googleapis.com/v1alpha/projects/earthengine-legacy/maps/5d6de8dd7c83fbf1b648270225721286-cd2a6107920090067cfbf242c86d1199/tiles/{z}/{x}/{y}","NDSI","NDSI",{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false}]},{"method":"addLayersControl","args":[["CartoDB.Positron","OpenStreetMap","CartoDB.DarkMatter","Esri.WorldImagery","OpenTopoMap"],["TRUE COLOR","NDSI"],{"collapsed":true,"autoZIndex":true,"position":"topleft"}]},{"method":"hideGroup","args":[null]},{"method":"addTiles","args":["https://earthengine.googleapis.com/v1alpha/projects/earthengine-legacy/maps/9747b8425ed6fe87650fa29e471a96ba-0451028497be24614f33299efacadbf7/tiles/{z}/{x}/{y}","NDVI","NDVI",{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false}]},{"method":"addLayersControl","args":[["CartoDB.Positron","OpenStreetMap","CartoDB.DarkMatter","Esri.WorldImagery","OpenTopoMap"],["TRUE COLOR","NDSI","NDVI"],{"collapsed":true,"autoZIndex":true,"position":"topleft"}]},{"method":"hideGroup","args":[null]},{"method":"addTiles","args":["https://earthengine.googleapis.com/v1alpha/projects/earthengine-legacy/maps/300626e7ce7e21a22e6a869ba42b0546-2a6882c3000d958d347c17297d07fa16/tiles/{z}/{x}/{y}","NDWI","NDWI",{"minZoom":0,"maxZoom":18,"tileSize":256,"subdomains":"abc","errorTileUrl":"","tms":false,"noWrap":false,"zoomOffset":0,"zoomReverse":false,"opacity":1,"zIndex":1,"detectRetina":false}]},{"method":"addLayersControl","args":[["CartoDB.Positron","OpenStreetMap","CartoDB.DarkMatter","Esri.WorldImagery","OpenTopoMap"],["TRUE COLOR","NDSI","NDVI","NDWI"],{"collapsed":true,"autoZIndex":true,"position":"topleft"}]},{"method":"hideGroup","args":[null]}],"setView":[[74.49,-20.49],10,[]]},"evals":[],"jsHooks":{"render":[{"code":"function(el, x, data) {\n  return (\n      function(el, x, data) {\n      // get the leaflet map\n      var map = this; //HTMLWidgets.find('#' + el.id);\n      // we need a new div element because we have to handle\n      // the mouseover output separately\n      // debugger;\n      function addElement () {\n      // generate new div Element\n      var newDiv = $(document.createElement('div'));\n      // append at end of leaflet htmlwidget container\n      $(el).append(newDiv);\n      //provide ID and style\n      newDiv.addClass('lnlt');\n      newDiv.css({\n      'position': 'relative',\n      'bottomleft':  '0px',\n      'background-color': 'rgba(255, 255, 255, 0.7)',\n      'box-shadow': '0 0 2px #bbb',\n      'background-clip': 'padding-box',\n      'margin': '0',\n      'padding-left': '5px',\n      'color': '#333',\n      'font': '9px/1.5 \"Helvetica Neue\", Arial, Helvetica, sans-serif',\n      'z-index': '700',\n      });\n      return newDiv;\n      }\n\n\n      // check for already existing lnlt class to not duplicate\n      var lnlt = $(el).find('.lnlt');\n\n      if(!lnlt.length) {\n      lnlt = addElement();\n\n      // grab the special div we generated in the beginning\n      // and put the mousmove output there\n\n      map.on('mousemove', function (e) {\n      if (e.originalEvent.ctrlKey) {\n      if (document.querySelector('.lnlt') === null) lnlt = addElement();\n      lnlt.text(\n                           ' lon: ' + (e.latlng.lng).toFixed(5) +\n                           ' | lat: ' + (e.latlng.lat).toFixed(5) +\n                           ' | zoom: ' + map.getZoom() +\n                           ' | x: ' + L.CRS.EPSG3857.project(e.latlng).x.toFixed(0) +\n                           ' | y: ' + L.CRS.EPSG3857.project(e.latlng).y.toFixed(0) +\n                           ' | epsg: 3857 ' +\n                           ' | proj4: +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs ');\n      } else {\n      if (document.querySelector('.lnlt') === null) lnlt = addElement();\n      lnlt.text(\n                      ' lon: ' + (e.latlng.lng).toFixed(5) +\n                      ' | lat: ' + (e.latlng.lat).toFixed(5) +\n                      ' | zoom: ' + map.getZoom() + ' ');\n      }\n      });\n\n      // remove the lnlt div when mouse leaves map\n      map.on('mouseout', function (e) {\n      var strip = document.querySelector('.lnlt');\n      if( strip !==null) strip.remove();\n      });\n\n      };\n\n      //$(el).keypress(67, function(e) {\n      map.on('preclick', function(e) {\n      if (e.originalEvent.ctrlKey) {\n      if (document.querySelector('.lnlt') === null) lnlt = addElement();\n      lnlt.text(\n                      ' lon: ' + (e.latlng.lng).toFixed(5) +\n                      ' | lat: ' + (e.latlng.lat).toFixed(5) +\n                      ' | zoom: ' + map.getZoom() + ' ');\n      var txt = document.querySelector('.lnlt').textContent;\n      console.log(txt);\n      //txt.innerText.focus();\n      //txt.select();\n      setClipboardText('\"' + txt + '\"');\n      }\n      });\n\n      }\n      ).call(this.getMap(), el, x, data);\n}","data":null},{"code":"function(el, x, data) {\n  return (function(el,x,data){\n           var map = this;\n\n           map.on('keypress', function(e) {\n               console.log(e.originalEvent.code);\n               var key = e.originalEvent.code;\n               if (key === 'KeyE') {\n                   var bb = this.getBounds();\n                   var txt = JSON.stringify(bb);\n                   console.log(txt);\n\n                   setClipboardText('\\'' + txt + '\\'');\n               }\n           })\n        }).call(this.getMap(), el, x, data);\n}","data":null}]}}</script>
```

### (9) Add NDSI, NDVI, NDMI and NDWI to the clipped image collection


```r
#Map normalized Difference functions over image collection
 MODIS_col <- MODIS_col$
   map(getNDSI)$
   map(getNDVI)$
   map(getNDMI)$
   map(getNDWI)
```

### (10) Create a timeseries gif of RGB images for the specified shapefile (Optional, for debugging)


```r
# #Check number of images in collection
#  MODIS_col$size()$getInfo()
# 
# #Create a timelapse video of RGB band
#  videoArgs <- list(dimensions=380, region=aoi,framesPerSecond=5, crs='EPSG:3857', bands=c("sur_refl_b01", "sur_refl_b04", "sur_refl_b03"),
#                    min=0, max=12000, gamma=c(1.9, 1.7, 1.7))
#  browseURL(MODIS_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$getVideoThumbURL(videoArgs))
# 
# #Create a timelapse video of NDSI band
#  palette=c('black', '0dffff', '0524ff', 'ffffff')
#  visFun_NDSI <-  function(img) {
#    return(img$visualize(bands='NDSI', min=-1, max=1.5, palette=palette)$copyProperties(img, img$propertyNames()))
#    }
#  MODIS_snow_RGB <- MODIS_col$map(function(img){return(img$clipToCollection(aoi_Shapefile))})$map(visFun_NDSI)
#  videoArgs <- list(dimensions=380, region=aoi, framesPerSecond=5, crs='EPSG:3857', bands=c('vis-red', 'vis-green', 'vis-blue'), min=0, max=255)
#  browseURL(MODIS_snow_RGB$getVideoThumbURL(videoArgs))
```




Built with R-version 4.3.0
