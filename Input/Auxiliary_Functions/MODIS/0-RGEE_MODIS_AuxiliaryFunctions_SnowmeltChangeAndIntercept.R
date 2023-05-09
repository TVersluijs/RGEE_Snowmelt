#This function calculates the change in a response variable over time for all pixels in an image (e.g. change in the date of 
#snowmelt over years). It does so by fitting a linear model for each pixel (y~year) and extracting the slope and intercept.
#To speed up processing of all pixels, the script runs the code in parallel (on multiple cores) for each defined subset of  
#pixels stored in the variable 'pixelIDs_split'. The output of each parallel run is stored in a separate position in a list. 
#Within each of these list elements the function returns two outputs, which are each stored in a separate sublist within the 
#main list. The first sublist contains a dataframe with the slope and intercept of the linear model for each pixel (of that 
#subset). The second sublist contains the plot of the response variable vs year and the fitted linear model as a red line.
f_ChangeInSnowmelt_FitLinearModel_parallel_2 <- function(subset, MODIS_pixel_y=MODIS_pixel_snowmelt, 
                                                         pixel_ID_column="pixel_ID", y="doy_snowmelt", x="Year"){
  
  #Select data for the current data subset (ranging from 1:length(pixel_IDs_split))
  MODIS_pixel_y_subset <- MODIS_pixel_y[MODIS_pixel_y[,pixel_ID_column] %in% pixel_IDs_split[[subset]],]
  print(paste0("Data subset ", subset, " out of ", length(pixel_IDs_split)))
  
  #Specify parameters of progress bar (normal printing of progress does not work within a foreach function)
  library(progress)
  pb <- progress_bar$new(
    format = "Pixel = :pixel [:bar] :elapsed | eta: :eta",
    total = length(unique(MODIS_pixel_y_subset[,pixel_ID_column])),
    width = 100)
  pb_func <- function(n){ pb$tick(tokens = list(pixel = unique(MODIS_pixel_y_subset[,pixel_ID_column])[n])) }
  
  #Specify a function to combine pixel-level dataframes and plots in a separate list
  comb <- function(x, ...){ lapply(seq_along(x), function(i){c(x[[i]], lapply(list(...), function(y){y[[i]]}))}) }
  
  #Specify a foreach (%dopar%) parallel loop through that (within the current subset) iterates through the data for all pixels
  results <- foreach(pixel=unique(MODIS_pixel_y_subset[,pixel_ID_column]), 
                     .combine = 'comb', .multicombine=TRUE,
                     .init=list(list(), list()),
                     .packages="ggplot2",
                     .options.snow=list(progress = pb_func)) %dopar% {
                       
                       # #for debugging:
                       # i=1
                       # pixel=unique(MODIS_pixel_y_subset[,pixel_ID_column])[i]
                       
                       #Select dataframe for current pixel
                       df_tmp <- MODIS_pixel_y_subset[MODIS_pixel_y_subset[,pixel_ID_column]==pixel,]
                       
                       #Change colnames of response and predictor variables to y and x respectively
                       colnames(df_tmp)[which(colnames(df_tmp)==y)] <- "y"
                       colnames(df_tmp)[which(colnames(df_tmp)==x)] <- "x"
                       
                       #Sort df_tmp by x
                       df_tmp <- df_tmp[order(df_tmp$x),]
                       
                       #If there are at least two doy_snowmelt measures available
                       if(sum(!is.na(df_tmp$y)) > 1){
                         
                         #Fit a linear model through the filtered pixel-specific y~Year data
                         mod_lm <- lm(y ~ I(x - 2010), data=df_tmp)
                         y_change <- mod_lm$coefficients['I(x - 2010)']
                         y_intercept <- mod_lm$coefficients['(Intercept)']
                         
                         #Use the fitted LM to make predictions of y based on a 1 year interval in Year
                         newdata <- base::data.frame(x=seq(min(df_tmp$x), max(df_tmp$x), 1))
                         newdata$lm_predict <- stats::predict(mod_lm, newdata=newdata, type="response")
                         
                         #Construct a plot with model predictions
                         p_tmp <- ggplot2::ggplot()+
                           geom_point(data=df_tmp, aes(x=x, y=y))+
                           geom_line(data=newdata, aes(x=x, y=lm_predict), col = "red")+
                           xlab(x)+
                           ylab(y)+
                           theme_classic()+
                           ggtitle(pixel)
                         
                       }
                       
                       #If there are less than two datapoints available
                       if(sum(!is.na(df_tmp$y)) < 2){
                         
                         y_change <- NA
                         y_intercept <- NA
                         
                         p_tmp <- ggplot() + 
                           annotate("text", x = NA, y = NA, size=8, label = "NA") + 
                           theme_void()+
                           ggtitle(pixel)
                         
                       }
                       
                       #Store the change in doy snowmelt (slope of linear regression) and intercept for each pixel in dataframe:
                       df_pixel_ChangeIn_Y <- base::data.frame(pixel_ID=pixel, 
                                                               y_change=y_change, 
                                                               y_intercept=y_intercept)
                       
                       #Store both the dataframe and plot in the multiresultclass object  
                       result <- list(df_pixel_ChangeIn_Y, p_tmp)
                       
                       #Output multiresult object  
                       return(result)
                     }
  
}
