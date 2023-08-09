#######################################################################################

# Automated RGEE scripts to analyze the date of snowmelt, and timeseries of NDSI, NDVI and NDMI based on MODIS or Sentinel-2 satellite data.

All scripts in this Github respository rely on the R-package RGEE by Cesar Aybar et al. My scripts provide an automated workflow
to extract the timing of snowmelt at the level of pixels, points, subareas, or larger polygons based on either MODIS or Sentinel-2 
satellite data. While the user can specify the required parameters (resolution, spatial extent, date ranges, cloud/water filtering
etc) the rest of the script is automated and generates the required output data. In this ReadMe file you can find (I) Usage notes,
(II) Which script to use? and (III) Detailed description of RGEE scripts.

Copyright: Tom S.L. Versluijs 2023 (tom.versluijs@gmail.com)

When using any of these scripts for scientific publications, please cite the original Github repository by following the provided
information in the file 'CITATION.dff' found in the root folder of this Github respository.
<br />
<br />

#######################################################################################

# (I) USAGE NOTES

#######################################################################################

-(1) Make sure to download the complete github folder 'RGEE' to make sure all dependencies between input files, scripts and output files function properly. The complete repository can be cloned
     to your own Github repository using '<>Code/Clone' or '<>Code/Download ZIP'. 

-(2) Update to the newest versions of R and R-Studio before running the scripts.

-(3) If RGEE has not yet been installed, or does not function properly, make sure to install RGEE and its python dependencies using the script "00-RGEE_TomVersluijs_Installation.R". Note that
     getting RGEE to function properly can be quite a frustrating first hurdle to take. Things will get easier once everything is up and running!
	 
-(4) All scripts in this Github folder need to be run from within the R-project "RGEE.Rproj" in the root directory. Thus, first open the 'RGEE.Rproj' file in RStudio, and then select
     the required script under 'Files/RGEE/R-scripts', or "/File/Open file..."
	 
-(5) All scripts rely on auxiliary scripts in the /RGEE/Input folder that specify commonly used functions. These auxiliary scripts are automatically sourced at the beginning of all main scripts. 

-(6) Several scripts rely on a shapefile of the study area as input. To construct these shapefiles follow the QGIS manual "Manual_QGIS_CreateShapefilePolygons" found in /RGEE/Manuals. After 
     creating your own shapefile(s), place them in the /RGEE/Input/Shapefiles directory.
	 
-(7) Each script generates its own output folder in the /RGEE/Output folder in the root directory.

-(8) All scripts are completely automated and require the user to only alter parameters in the "#Specify parameters of interest" section at the start of each script. After adjusting these 
     parameters the complete code can be selected (CTRL-A) and can be run at once (CTRL + ENTER, or CTRL + R, or by clicking Run in R-Studio).
	 
-(9) Sentinel-2 scripts can be computationally demanding and depending on the occupancy of the Google Earth Engine servers might take up to 16 hours to complete! The servers can sometimes be
     so busy that the script returns a 'computation time out' or 'computation error'. If this happens, please try to run the script again at some other time. If the code still results in an
	 computation error then your spatial extend is probably to large for the specified resolution. Either decrease the spatial resolution of the analysis by setting the parameter 'resolution',
	 or reduce the spatial extent of the analysis.
	 
-(10) For questions regarding potential errors and bugs please contact tom.versluijs@gmail.com.
<br />
<br />

#######################################################################################

# (II) WHICH SCRIPT TO USE?

#######################################################################################

For a more thorough description of all RGEE scripts, see the section 'DETAILED DESCRIPTION OF RGEE SCRIPTS' below.

### General comparison of MODIS and Sentinel-2
The R-scripts provided below analyze data using either MODIS or Sentinel-2 satellite data. The choice between either satellite depends
mostly on the required spatial resolution (MODIS=500m, sentinel-2=10m), and how far back in time you want to analyze these data
(MODIS since 2000, Sentinel-2 since 2016). Moreover, if you aim to analyze very large spatial extents (>50-200km2) then Sentinel-2
scripts often result in a computation error (even when manually specifying a lower spatial resolution). In that case MODIS data
is a better choice. Below, you can find a few queries to aid in picking the right script for your research question.

### Which years do you want to analyse?
* MODIS: 	 	
	* years: 2000 - current
	* scripts: 01 - 03

* Sentinel-2: 	
	* years: 2016 - current
	* scripts: 04 - 09

### What is the required spatial resolution of your analysis?
* MODIS: 	 	
	* resolution: 500 meter			
	* scripts: 01 - 03

* Sentinel-2:  	
	* resolution: 10 meter
	* scripts: 04 - 09

### MODIS scripts: what kind of data do you want to analyse?
* Date of snowmelt, and average NDSI, NDVI and NDMI at multiple point locations (with or without buffer)				
	* Script 01

* Date of snowmelt at all pixels within a shapefile area																
	* Script 02

* Change in date of snowmelt per pixel within a shapefile area 															
	* Script 03


### Sentinel-2 scripts: what kind of data do you want to analyse?
* Date of snowmelt, and average NDSI, NDVI and NDMI at multiple point locations (with or without buffer)									
	* Script 04

* Date of snowmelt, and average NDSI, NDVI and NDMI within a small shapefile (<50km2)									
	* Script 05

* Date of snowmelt, and average NDSI, NDVI and NDMI at subareas located within a small shapefile area (<50km2)			
	* Script 06

* Date of snowmelt, and average NDSI, NDVI and NDMI at multiple point locations (with or without buffer) located within a small shapefile area (<50km2)	
	* Script 07

* Date of snowmelt for all pixels within a small shapefile (area < 50km2)												
	* Script 08
	
* Date of snowmelt for all pixels within a large shapefile (area < 250km2)												
	* Script 09
<br />
<br />

#######################################################################################

# (III) DETAILED DESCRIPTION OF RGEE SCRIPTS

#######################################################################################

## GENERAL SCRIPTS:

### *00-RGEE_TomVersluijs_Installation.R
Script to install RGEE by manually setting the Python path. This script is adapted from a tutorial by Ricardo Dal'Agnol da Silva.

## MODIS SCRIPTS:

### *01-RGEE_TomVersluijs_MODIS_Points_Snowmelt.R
In this script the date of snow melt is calculated based on MODIS data for all locations specified in an input file. The user 
can specify a bufferzone (radius) to depict the area in which snow melt will be analysed per location. All locations are analysed 
consecutively (using a loop). First, a location specific bounding box is drawn per point location taking into account the specified 
buffer zone) and MODIS satellite data is extracted within this bounding box. Second, clouds and permanent water bodies are filtered 
within this bounding box. Finally, snow melt is analysed within each locations's buffer zone based on one of the following methods 
(specified by the user by setting the parameter 'method'): (1) fitting a GAM through the average NDSI, NDVI and NDMI values over 
time, (2) Fitting a GAM through the fraction of snow-covered pixels over time, or (3) Fitting a GAM through the NDSI values per 
pixel. No shapefile is required as input for this script, which allows for input locations to be spaced far apart.					

<p float="left">
  <img align="top" src="./_pictures/01A-MODIS_Wrangel_Locations.png" width="49%" title="Five locations at Wrangel Island, each with a buffer of 5000m" />
  <img align="top" src="./_pictures/01B-MODIS_Wrangel_Location4_Pixels_NDSI_GAMS.png" width="49%" title="Method 3: GAMS fitted through NDSI values per pixel for a single location" /> 
 </p>
 <p float="left">
  <img align="top" src="./_pictures/01C-MODIS_Wrangel_Location4_Image_Snowmelt.png" width="49%" title="Method 3: Image of the pixel-specific date of snowmelt for a single location" />
  <img align="top" src="./_pictures/01D-MODIS_Wrangel_Locations_Snowmelt_GAMS.png" width="49%" title="Method 3: GAMS fitted through the fraction of snow-covered pixels for each location" />
</p>
<br />
<br />
				
### *02-RGEE_TomVersluijs_MODIS_Shapefile_Pixel_Snowmelt.R
Create pixel-level maps (500m resolution) of the date of snowmelt in a study area based on MODIS data. Snowmelt is calculated 
per pixel by fitting a GAM through the average NDSI data and extracting the moment this GAM crosses a user specified NDSI threshold. 
This script requires a shapefile of the study area as input. The user can specify whether clouds and permanent waterbodies need to
be masked from the analysis.

<p float="left">
  <img align="top" src="./_pictures/05A-Sentinel2_Zackenberg_Shapefile.png" width="49%" title="Shapefile area for Zackenberg (c.a. 50km2)" />
  <img align="top" src="./_pictures/02A-MODIS_Zackenberg_Pixels_NDSI_GAMS.png" width="49%" title="GAMS fitted through NDSI values for each pixel within the shapefile" />
</p>
<p float="left">
  <img align="top" src="./_pictures/02B-MODIS_Zackenberg_Image_Snowmelt.png" width="49%" title="Image of the date of snowmelt for all pixels within the shapefile (500m resolution)" /> 
</p>
<br />
<br />	
													   
### *03-RGEE_TomVersluijs_MODIS_Shapefile_Pixel_ChangeInSnowmelt.R
This script requires MODIS snowmelt maps generated using script "02-RGEE_TomVersluijs_MODIS_Shapefile_Pixel_Snowmelt.R" as input.
It imports the pixel-level snowmelt images for all analyzed years and transforms them into an image with the change in the date 
of snowmelt over the years for each pixel (i.e. slope of linear regression) and another image with the average date of snowmelt 
over the years for each pixel (i.e. intercept of linear regression).

<p float="left">
  <img align="top" src="./_pictures/05A-Sentinel2_Zackenberg_Shapefile.png" width="49%" title="Shapefile area for Zackenberg (c.a. 50km2)" />
  <img align="top" src="./_pictures/03A-MODIS_Zackenberg_Pixels_Regressions_ChangeInSnowmelt.png" width="49%" title="Linear regressions fitted to pixel-specific timeseries of the date of snowmelt" />
</p>
<p float="left">
  <img align="top" src="./_pictures/03B-MODIS_Zackenberg_Image_ChangeInSnowmelt.png" width="49%" title="Image of the change in date of snowmelt (slope) for all pixels within the shapefile" /> 
  <img align="top" src="./_pictures/03C-MODIS_Zackenberg_Image_SnowmeltIntercept.png" width="49%" title="Image of the average date of snowmelt (intercept at average year) for all pixels within the shapefile" /> 
</p>
<br />
<br />

## SENTINEL-2 SCRIPTS:

### *04-RGEE_TomVersluijs_S2_Points_Snowmelt.R
In this script the date of snow melt is calculated based on Sentinel-2 data for all locations specified in an input file. The 
user can specify a bufferzone (radius) to depict the area in which snow melt will be analysed per location. All locations are 
analysed consecutively (using a loop). First, a location specific bounding box is drawn per point location (taking into account 
the specified buffer zone) and Sentinel-2 satellite data is extracted within this bounding box. Second, clouds and permanent 
water bodies are filtered within this bounding box. Third, if this bounding box overlaps with multiple satellite tiles for a 
certain day, a composite image is created (picking the pixel with least cloud cover). Finally, snow melt is analysed within each 
locations's buffer zone based on one of the following methods (specified by the user by setting the parameter 'method'): (1) 
fitting a GAM through the average NDSI, NDVI and NDMI values over time, (2) Fitting a GAM through the fraction of snow-covered 
pixels over time. Note that calculations of the date of snowmelt per pixel are not conducted in this script because these are
too memory intensive for this high resolution sentinel-2 data and will result in a computation error. No shapefile is required 
as input for this script, which allows for input locations to be spaced far apart.										

<p float="left">
  <img align="top" src="./_pictures/04A-Sentinel2_Zackenberg_Locations.png" width="49%" title="Five locations at Zackenberg, each with a buffer of 250" />
  <img align="top" src="./_pictures/04B-Sentinel2_Zackenberg_Locations_NDSI_Snowmelt.png" width="49%" title="Method 1: GAMS fitted through average NDSI values per location" /> 
</p>
<p float="left">
  <img align="top" src="./_pictures/04C-Sentinel2_Zackenberg_Locations_SnowFraction_Snowmelt.png" width="49%" title="Method 2: GAMS fitted through the fraction of snow-covered pixels per location" />
</p>
<br />
<br />	
								
### *05-RGEE_TomVersluijs_S2_Shapefile_Snowmelt.R 
Use Sentinel-2 data to extract time series of the average NDVI, NDMI, and NDSI and the fraction of snowcover for a single area of 
interest (i.e. a single shapefile). The user can specify whether clouds and permanent waterbodies need to be masked, and whether
a composite image per day of year needs to be generated (merging multiple satellite photos for that day). The current script 
requires a shapefile of the study area as input. It only works for small areas of c.a. 50 km2 (larger areas might result in computation 
errors).

<p float="left">
  <img align="top" src="./_pictures/05A-Sentinel2_Zackenberg_Shapefile.png" width="49%" title="Shapefile area for Zackenberg (c.a. 50km2)" />
  <img align="top" src="./_pictures/05B-Sentinel2_Zackenberg_Shapefile_SnowFraction.png" width="49%" title="GAM fitted through fraction of snowcover over time within the shapefile area" /> 
</p>
<br />
<br />

### *06-RGEE_TomVersluijs_S2_Shapefile_SubAreas_Snowmelt.R
Use Sentinel-2 data to extract time series of the average NDVI, NDMI, and NDSI and the fraction of snowcover for all sub areas located
within a shapefile (these sub areas can be specified by creating a multipolygon in e.g. QGIS). The user can specify whether clouds and 
permanent waterbodies need to be masked, and whether a composite image per day of year needs to be generated (merging multiple satellite 
photos for that day). The current script requires a shapefile including subareas within the study area as input. It only works for small 
areas of c.a. 50 km2 (larger areas might result in computation errors).

<p float="left">
  <img align="top" src="./_pictures/06A-Sentinel2_Zackenberg_Shapefile_SubAreas.png" width="49%" title="Shapefile area for 10 subareas in Zackenberg (each c.a. 5km2)" />
  <img align="top" src="./_pictures/06B-Sentinel2_Zackenberg_Shapefile_SubAreas_FractionSnowCover.png" width="49%" title="GAM fitted through fraction of snowcover over time within each subarea" /> 
</p>
<br />
<br />

### *07-RGEE_TomVersluijs_S2_Shapefile_Points_Snowmelt.R
In this script the date of snow melt is calculated based on Sentinel-2 data for all locations specified in an input file. The
user can specify a bufferzone (radius) to depict the area in which snow melt will be analysed per location. All locations are
required to be located within a single shapefile and are then analysed simultaneously. First, clouds and permanent water bodies 
are filtered within the shapefile. Second, if the shapefile overlaps with multiple satellite tiles for a certain day, a composite
image is created (picking the pixel with least cloudcover). Finally, snow melt is analysed for each location's buffer zone based 
on one of the following methods (specified by the user by setting the parameter 'method'): (1) fitting a GAM through the average 
NDSI, NDVI and NDMI values over time, (2) Fitting a GAM through the fraction of snow-covered pixels over time. Note that this 
script is similar to script "*04-RGEE_TomVersluijs_S2_Points_Snowmelt.R". However, in the latter script all points are analysed 
consecutively using a loop, which makes that script significantly slower to run. It does however not rely on a shapefile and 
thus works for points spaced much further apart (e.g. tracking data of birds). The current script only works for small areas
of c.a. 50 km2 (larger areas might result in computation errors).
	
<p float="left">
  <img align="top" src="./_pictures/07A-Sentinel2_Zackenberg_Shapefile_Points.png" width="49%" title="Five point locations with a buffer zone of 250m located within a shapefile area in Zackenberg" />
  <img align="top" src="./_pictures/07B-Sentinel2_Zackenberg_Shapefile_Points_FractionSnowCover.png" width="49%" title="GAM fitted through fraction of snowcover over time within each location's bufferzone" /> 
</p>
<br />
<br />	

### *08-RGEE_TomVersluijs_S2_Shapefile_Pixel_Snowmelt.R
Use Sentinel-2 data to create pixel-level maps (10m resolution) of the date of snowmelt for an area up to c.a. 50km2. Snowmelt
is calculated per pixel by fitting a GAM through the average NDSI data and extracting the moment this GAM crosses a user specified
NDSI threshold. This script requires a single shapefile of the study area as input. The user can specify whether clouds and
permanent waterbodies need to be masked. No composite image can be generated because this will result in computation errors.

<p float="left">
  <img align="top" src="./_pictures/05A-Sentinel2_Zackenberg_Shapefile.png" width="49%" title="Shapefile area for Zackenberg (c.a. 50km2)" />
  <img align="top" src="./_pictures/08A_Sentinel2_Zackenberg_Shapefile_Image_Snowmelt.png" width="49%" title="Image of the date of snowmelt at Zackenberg for all pixels within the shapefile (10m resolution)" />
</p>								
<br />
<br />

### *09-RGEE_TomVersluijs_S2_Shapefile_Pixel_Snowmelt_LargeAreas.R
Use Sentinel-2 data to create pixel-level maps (10m resolution) of the date of snowmelt for an area up to c.a. 200km2. Snowmelt
is calculated per pixel by fitting a GAM through the average NDSI data and extracting the moment this GAM crosses a user specified
NDSI threshold. This script requires that the shapefile of the study area is split-up into several smaller shapefiles to prevent
memory issues on the GEE-server. No composite image can be generated because this will result in computation errors. This script
currently only works by splitting the main shapefile into four smaller shapefiles.

<p float="left">
  <img align="top" src="./_pictures/09A-Sentinel2_Taymir_Shapefile_Subareas.png" width="49%" title="Shapefile split up into 4 subareas for Taymir (c.a. 250km2)" />
  <img align="top" src="./_pictures/09B_Sentinel2_Taymir_Shapefile_Image_Snowmelt.png" width="49%" title="Image of the date of snowmelt at Taymir for all pixels within the shapefile (10m resolution)" />
</p>																	