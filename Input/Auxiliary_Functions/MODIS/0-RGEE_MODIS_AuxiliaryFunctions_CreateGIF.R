####################################################################################################################################
#
#                         Function to create a GIF animation from an image collection
#
#################################################################################################################################### 

#DESCRIPTION: This script contains the function 'f_img_col_to_gif' that transforms an image collection to a date-annotated
#GIF-animation and automatically saves the output to the output folder. In addition it contains a visualization function
#'f_band_to_RGB' to transform a single-band image (like NDSI) to an RGB-image as this is required for the GIF.

####################################################################################################################################

#(I): Transform an image collection to a text-annotated GIF and store it in the output folder
f_img_col_to_gif <- function(img_col=image_collection,  #image collection
                             RGB_bands=c("name_of_red_band", "name_of_green_band", "name_of_blue_band"), #names of RGB-bands in image collection
                             shapefile=aoi_Shapefile, #shapefile of area of interest
                             centroid_buffer_m=0, #buffer zone surrounding centroid of shapefile in meters (to enlarge area for GIF)
                             gif_dimensions=gif_max_pixels, #Maximum dimensions of GIF (pixels)
                             gif_fps=5, #Frames per second of GIF
                             gif_min=0, #Value to map to 0 (to improve contrast). A good rule of thumb is to set min to values that represent the 2nd percentile of the data
                             gif_max=8000, #Value to map to 255 (to improve contrast). A good rule of thumb is to set max to values that represent the 98th percentile of the data
                             gif_gamma=c(1.9, 1.7, 1.7), #Gamma correction factors (one for each band)
                             gif_crs='EPSG:3857', #CRS project of the output
                             gif_text_property="Date", #Name of property in each image which will be used to annotate each frame of the GIF
                             gif_text_position="northwest", #Location of text
                             gif_text_position_adjustment="+0+0", #Small scale adjustment of text in meters
                             gif_text_size=14, #Size of text
                             gif_text_col="#FFFFFF", #Colour of text
                             output_fldr="Path_of_output_folder", #Name of output folder
                             file_name="GIF_RGB"){ #Name of output file
                             #data_ID=data_ID, #previously created variable (required for naming of output file)
                             #resolution=resolution, #previously created variable (required for naming of output file)
                             #year_ID=year_ID, #previously created variable (required for naming of output file)
                             #timestamp=timestamp #previously created variable (required for naming of output file)
  
  #(0): Print progress
    print(paste0("Start creating GIF: ")) 
    print(paste0(output_fldr, timestamp, "_", data_ID, "_Res", resolution, "_", file_name, ".gif"))
    cat("\n")
  
  #(1): Clip images by the shapefile area
    
    #Select only the bands of interest
    img_col <- img_col$select(RGB_bands)
  
    #Add buffer of size 'centroid_buffer_m' to the shapefile area
    if(centroid_buffer_m>0){
      shapefile_large <-  shapefile$geometry()$centroid()$buffer(centroid_buffer_m)
      img_col_gif <- img_col$map(function(img){return(img$clipToCollection(ee$FeatureCollection(shapefile_large)))})
      }
    
    #Do not add a buffer to the shapefile area
    if(!centroid_buffer_m>0){
      shapefile_large <- shapefile$geometry()
      img_col_gif <- img_col$map(function(img){return(img$clipToCollection(ee$FeatureCollection(shapefile_large)))})
      }
  
  #(2): Overlay a polygon of the shapefile area on top of each image  
  
    #Transform shapefile into a polygon layer with styling (e.g., color, fillColor, line width and RGB bandnames)
    polygon_layer <- shapefile$style(color="red", fillColor="#00000000", width=2)$rename(RGB_bands)
    
    #Make sure the polygon layer has the same bands (and order) as the image collection. We therefore only
    #retain the RGB bands from the img_col_gif:
    img_col_gif <- img_col_gif$map(function(img){return(img$select(RGB_bands))})
    
    #Overlay polygon layer on each image in the image collection
    img_col_gif_with_polygon <- img_col_gif$map(function(img) {return(img$blend(polygon_layer))})
    
    # #Visualize a single image (for debugging)
    # image <- img_col_gif_with_polygon$filterDate(paste0(year_ID, "-07-23"), end_date)$first()
    # Map$setCenter(coordinates_point$getInfo()$coordinates[1], coordinates_point$getInfo()$coordinates[2], 10)
    # Map$addLayer(image,list(bands=c(RGB_bands), min=100, max=8000, gamma=c(1.9, 1.7, 1.7)), 'TRUE COLOR')
  
  #(3): Extract 'gif_text_property' of each image (to annotate each frame of the GIF)
  
    #create a current timestamp to prevent identical names on Google Drive
    current_timestamp_gif1 <- gsub('\\.', '', format(Sys.time(), "%y%m%d%H%M%OS2"))
    
    #Setup task
    task_vector_gif1 <- ee_table_to_drive(
      collection= img_col, #using img_col instead of img_col_gif_with_polygon to speed up code 
      description = paste0(current_timestamp_gif1, "_", data_ID, "_Res", resolution, "_GIF_Text"),
      folder="RGEE_tmp",
      fileFormat="CSV",
      selectors=gif_text_property
      )
  
    #Run and monitor task
    task_vector_gif1$start()
    ee_monitoring(task_vector_gif1, quiet=F, task_time=10, max_attempts=1000000)
    
    #Import results
    exported_stats <- ee_drive_to_local(task=task_vector_gif1, dsn=paste0(output_fldr, timestamp, "_", data_ID, "_Res", resolution, "_GIF_Text"))
    gif_text <- read.csv(exported_stats)  ;   unlink(exported_stats)
    gif_text <- gif_text[,gif_text_property]
    
  #(4): Create GIF
  
    #Catch potential errors
    
      #Create and annotate GIF
      tryCatch({
        videoArgs <- list(dimensions=gif_dimensions, region=shapefile_large, framesPerSecond=gif_fps, crs=gif_crs, bands=c(RGB_bands), min=gif_min, max=gif_max, gamma=gif_gamma)
        animation <- rgeeExtra::ee_utils_gif_creator(img_col_gif_with_polygon, quiet=F, videoArgs, mode="wb")
        animation <- animation %>% rgeeExtra::ee_utils_gif_annotate(gif_text, gravity=gif_text_position, location=gif_text_position_adjustment, size=gif_text_size, color=gif_text_col)
        gc(reset = TRUE)
        rgeeExtra::ee_utils_gif_save(animation, path = paste0(output_fldr, timestamp, "_", data_ID, "_Res", resolution, "_", file_name, ".gif"))
        }, error = function(cond){return({
          print(paste0("Unable to create: '", file_name, ".gif'"))
          print(paste0("  -Detailed error: ", cond))
          print("  -In case of error - 'Total request size (X pixels) must be less than or equal to X pixels' - reduce 'gif_max_pixels'.")
          cat("\n")
          })})
    
}

#(II): Transform each image in an image collection from a single band (e.g. NDSI) to an RGB image

#Map a visualization function over the image collection. This function converts each image to an RGB image (i.e. a red, green and blue
#colour band with min and max values of 0 and 255 respectively). In the function it can be specified which data band needs to be converted,
#what min and max values of the original band should correspond to 0 and 255 in the RGB image, and which colour palette should be used.
f_band_to_RGB <-  function(img, 
                           band='NDSI', #Name of the band that needs to be converted to RGB
                           min_value=-1, #Value of band to map to 0
                           max_value=1.5, #Value of band to map to 255
                           palette=c('#000000', '#0dffff', '#0524ff', '#ffffff')){ #colour palette for converted band
  return(img$
           visualize(bands=band, min=min_value, max=max_value, palette=palette)$
           copyProperties(img, img$propertyNames()))}
