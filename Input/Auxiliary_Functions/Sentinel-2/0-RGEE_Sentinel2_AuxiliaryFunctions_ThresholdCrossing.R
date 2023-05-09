####################################################################################################################################
#
#                  Sentinel 2 - Function for detecting the moment of crossing of threshold per pixel using GAMs
#
#################################################################################################################################### 

#DESCRIPTION: This script contains the function 'f_detect_threshold_date_parallel_1' that calculates per pixel for which day of 
#year a threshold value in the response variable is crossed (e.g. the date of snowmelt when NDSI is the response variable). This
#is done by fitting a GAM through the data (response~doy) and interpolating the moment this threshold is crossed in the direction 
#from above to below the threshold value (for snowmelt, using parameter direction="down").To speed up processing of all pixels, 
#the script runs the code in parallel (on multiple cores) for each defined subset of pixels stored in the variable 'pixelIDs_
#split'. The output is a dataframe with the moment the threshold was crossed for each pixel (within the respective pixel 
#subset), which is stored in a consecutive position in a list. 

#The function f_detect_threshold_date_parallel_1()' takes as input the number of the subset (provided using an external lapply 
#function), dataframe 'df_pixel_y', the column in df_pixel_y where the pixelID is specified (pixel_ID_column), the name of the 
#dependent  variable in df_pixel_y (y), the name of the independent variable in df_pixel_y (x), the threshold used for the 
#dependent variable (y_threshold).

####################################################################################################################################

#(I): Specify function to define the moment of crossing of the threshold
f_detect_threshold_date_parallel_1 <- function(subset, df_pixel_y=df_pixel_ndsi, pixel_ID_column="pixel_ID", y="NDSI", 
                                               x="doy", y_threshold=NDSI_threshold){
  
  #Select data for the current data subset (ranging from 1:length(pixelIDs_split))
  df_subset <- df_pixel_y[df_pixel_y[,pixel_ID_column] %in% pixelIDs_split[[subset]],]
  print(paste0("Data subset ", subset, " out of ", length(pixelIDs_split)))
  
  #Specify parameters of progress bar (normal printing of progress does not work within a foreach function)
  library(progress)
  pb <- progress_bar$new(
    format = "Pixel = :pixel [:bar] :elapsed | eta: :eta",
    total = length(unique(df_subset[,pixel_ID_column])),
    width = 100)
  pb_func <- function(n){
    pb$tick(tokens = list(pixel = unique(df_subset[,pixel_ID_column])[n]))
  }
  
  #Specify a foreach (%dopar%) parallel loop through that (within the current subset) iterates through the data for all pixels
  results <- foreach(pixel=unique(df_subset[,pixel_ID_column]), 
                               .export=c("year_ID", "outlier_thresh_1", "outlier_thresh_2", "outlier_removal", "gam_k_outlier", "gam_k", "f_gam_SeqRemOutliers"), #"y_threshold"
                               #.export specifies all variables needed from outside this function
                               .combine=rbind, 
                               .options.snow=list(progress = pb_func)) %dopar% {
                                 
                                 # #for debugging:
                                 # i=1
                                 # pixel=unique(df_subset[,pixel_ID_column])[i]
                                 # outlier_removal=TRUE                  
                                 
                                 #Select dataframe for current pixel
                                 df_tmp <- df_subset[df_subset[,pixel_ID_column]==pixel,]
                                 
                                 #Filter dataset by sequentially removing outliers using the function 'f_gam_SeqRemOutliers'
                                 df_tmp <- f_gam_SeqRemOutliers(data=df_tmp, y=y, x=x, outlier_removal=outlier_removal, 
                                                                outlier_thresh_1=outlier_thresh_1, outlier_thresh_2=outlier_thresh_2,
                                                                default_k=gam_k_outlier)
                                 
                                 #Sort df_tmp by "x" (i.e. by day of year)
                                 df_tmp <- df_tmp[order(df_tmp[,x]),]
                                 
                                 #Change column names of response and predictor variable to y and x respectively
                                 colnames(df_tmp)[colnames(df_tmp)==y] <- "y"
                                 colnames(df_tmp)[colnames(df_tmp)==x] <- "x"
                                 
                                 #Fit a GAM through the filtered pixel-specific y~x data (exclude outliers)
                                 index <- which(df_tmp$outliers==FALSE)
                                 mod_gam <- mgcv::gam(y ~ s(x, k=min(gam_k, length(index)-1)), data = df_tmp[index, ], method="REML")
                                 
                                 #Use the fitted GAM to make predictions of y based on a 1 day interval in x
                                 newdata <- base::data.frame(x=seq(min(df_tmp$x), max(df_tmp$x), 1))
                                 newdata$gam_predict <- stats::predict(mod_gam, newdata=newdata, type="response")
                                 
                                 #Detect cutoff points where curve goes from above y_threshold to below y_threshold
                                 newdata$cutoff <- ifelse(newdata$gam_predict >= y_threshold, 1, 0)
                                 newdata$dif <- c(0, diff(newdata$cutoff))
                                 #the column 'cutoff' indicates whether the gam prediction is above (1) or below (0) the y_threshold
                                 #the column 'dif' indicates when there is a change from 1 to 0 (-1) or 0 to 1 (1) in the column cutoff
                                 #Thus, those rows where the column 'dif' is equal to -1 indicate moments where the y value changes from above
                                 #the threshold to below the threshold. It might be possible that this happens multiple times within a season due to
                                 #measurement errors or cloud effects. We therefore need to determine which 'cutoff' most likely corresponds to the 
                                 #final moment of threshold crossing (e.g. the actual moment of snowmelt)
                                 
                                 #If y_threshold was at least crossed once:
                                 if(any(newdata$dif<0)){
                                   
                                   #Select all moments (cutoffs) where dif==-1
                                   cutoffs <- data.frame(index=which(newdata$dif<0))
                                   
                                   #For the period 30 days after each cutoff point, sum the number of days that have a y value larger than y_threshold. If a 
                                   #cutoff represents the final moment the threshold is crossed, then we do not expect any days after this moment with y > 
                                   #y_threshold. Thus, the closer this sum is to 0, the more likely this cutoff corresponds to the final moment  that this 
                                   #threshold is crossed (i.e. the actual moment of snowmelt).
                                   cutoffs$min <- cutoffs$index -30
                                   cutoffs$min[cutoffs$min<1] <- 1
                                   cutoffs$max <- cutoffs$index + 29
                                   cutoffs$max[cutoffs$max>nrow(newdata)] <- nrow(newdata)
                                   cutoffs$sum_cutoff_plus_30 <- apply(cutoffs, 1, function(x){sum(newdata$cutoff[x['index']:(x['max'])])})
                                   cutoff_best <- cutoffs[cutoffs$sum_cutoff_plus_30==min(cutoffs$sum_cutoff_plus_30),'index'][1]
                                   
                                   #Approximate day on which threshold is crossed in period from (cutoff_best-1 : cutoff_best) 
                                   newdata_subset <- newdata[max(0, cutoff_best-2) : min(cutoff_best+1, nrow(newdata)),]
                                   x_threshold <- stats::approx(x = newdata_subset$gam_predict, y = newdata_subset$x, xout = y_threshold)$y[1]
                                   
                                   #For debugging:
                                   # p_tmp <- ggplot2::ggplot()+
                                   #   geom_point(data=df_tmp[index,], aes(x=x, y=y))+
                                   #   geom_point(data=df_tmp[df_tmp$outliers==TRUE,], aes(x=x, y=y), col="black", pch=16, alpha=0.2)+
                                   #   geom_line(data=newdata, aes(x=x, y=gam_predict), col = "red") +
                                   #   geom_point(aes(x=x_threshold, y=y_threshold), col="blue", size=3)+
                                   #   geom_hline(yintercept=y_threshold, lty=2, col="grey")+
                                   #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                                   #   ylab(paste0(y, "-value at pixel")) +
                                   #   theme_classic()+
                                   #   ggtitle(pixel)
                                   
                                 }
                                 
                                 #If threshold was not crossed:
                                 if(!any(newdata$dif<0)){
                                   
                                   #There was no date at which the threshold was crossed (set x_threshold to NA)
                                   x_threshold <- NA
                                   
                                   #For debugging:
                                   # p_tmp <- ggplot2::ggplot()+
                                   #   geom_point(data=df_tmp[index,], aes(x=x, y=y))+
                                   #   geom_point(data=df_tmp[df_tmp$outliers==TRUE,], aes(x=x, y=y), col="black", pch=16, alpha=0.2)+
                                   #   geom_line(data=newdata, aes(x=x, y=gam_predict), col = "red") +
                                   #   geom_hline(yintercept=y_threshold, lty=2, col="grey")+
                                   #   xlab(paste0("Day of year (starting at 01-01-", year_ID, ")")) +
                                   #   ylab(paste0(y, "-value at pixel")) +
                                   #   theme_classic()+
                                   #   ggtitle(pixel)
                                   
                                 }
                                 
                                 #for debugging:
                                 #print(p_tmp)
                                 
                                 #Store date at which threshold was crossed (e.g. date of snowmelt in case of y=NDSI) for each pixel in dataframe:
                                 df_pixel_threshold <- base::data.frame(pixel_ID=pixel, x_threshold=x_threshold)
                                 return(df_pixel_threshold)
                               }
  
}
 
 
###############################################################################################################################################