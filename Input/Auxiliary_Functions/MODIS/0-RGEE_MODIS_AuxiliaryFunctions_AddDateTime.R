####################################################################################################################################
#
#                                              MODIS - Add image properties
#
####################################################################################################################################
 
#DESCRIPTION: This script contains functions that can be mapped over an MODIS image collection to add datetime and day of year as 
#image properties to all images in the collection.

####################################################################################################################################

#(I): Add datetime and day of year (doy) to an image
 add_Date = function(img){
  date <- img$date()$format("YYYY-MM-dd hh:mm:ss")
  doy <- ee$Date(img$get('system:time_start'))$getRelative('day', 'year')
  return(img$
           set('Date', date)$
           set("doy", doy))
  }

#(II): Add 'seconds since 1st of january' as a new band to each image 
 add_Time <- function(img){
   seconds <- ee$Number(ee$Date(img$get('system:time_start'))$getRelative('second', 'year'))
   seconds_img <- ee$Image$constant(ee$Number(seconds))$toUint32()$rename('seconds') 
   return(img$addBands(seconds_img))
 } 
 
 
####################################################################################################################################

