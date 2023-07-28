####################################################################################################################################
#
#                            Sentinel 2 - Function for sequential outlier detection using GAMs
#
#################################################################################################################################### 

#DESCRIPTION: This script contains the function f_gam_SeqRemOutliers() that detects outliers in a dataset by fitting GAMS and 
#employing a sequential residual-filtering process. The initial step is to fit a GAM through the data, and calculate model 
#predictions and residuals. We then mark all datapoints as outliers of which residuals are >= (outlier_thresh_1 * the range of  
#the y-data). We then re-fit a GAM to the same dataset excluding these initial outliers. We again make predictions and calculate 
#residuals. We then mark all datapoints as additional outliers where residuals >= (outlier_thresh_2 * the range of the y-data).  
#This gives us a final dataframe in which outliers have been marked following this two-step approach. Note that a sequential step 
#is required because initially some datapoints might falsely be assigned a large residual because of some extreme outliers. 
#After the initial removal of these extreme outliers and refitting of the GAM, it can be better assessed which data points
#truly have a large residual and can thus be assigned as 'actual' outliers.

#The function f_gam_SeqRemOutliers() takes as input a dataframe (data), the name of the dependent variable (y), the name 
#of the independent variable (x), whether outlier removal should be employed (outlier_removal), an initial outlier threshold 
#(outlier_thresh_1), a secondary outlier threshold (outlier_threshold_2), and the parameter 'k' indicating the minimum number 
#of knots used in fitting of the GAMs.

####################################################################################################################################

#(I): Specify function to sequentially remove oultiers from a GAM
 f_gam_SeqRemOutliers <- function(data, y, x, outlier_removal=outlier_removal, 
                                  outlier_thresh_1=outlier_thresh_1, 
                                  outlier_thresh_2=outlier_thresh_2,
                                  default_k=10){
    
    #Remove all rows with NA in columns 'x' or 'y'
     data <- data[!is.na(data[,x]) & !is.na(data[,y]),] 
   
    #Find range of data
     data_range <- (max(data[,y]) - min(data[,y]))
    
    #Fit initial GAM through data
     mod_gam <- mgcv::gam(data[,y] ~ s(data[,x], k=min(default_k, nrow(data)-1)), method="REML")
     data$gam_residuals <- base::round(mod_gam$residuals, 4)
     data$outliers = FALSE
    
    #Sequentially detect outliers (first on outlier_thresh_1, and after refitting a GAM, on outlier_thresh_2)
     if(outlier_removal==TRUE){
       
       #Mark outliers where Residuals >= outlier_thresh_1 * data_range
        data$outliers[abs(data$gam_residuals) >= (outlier_thresh_1 * data_range)] <- TRUE
       
       #Refit a gam to the reduced data
        index <- which(data$outliers==FALSE)
        mod_gam <- mgcv::gam(data[index,y] ~ s(data[index,x], k=min(default_k, length(index)-1)), method="REML")
        data$gam_residuals[index] <- base::round(mod_gam$residuals, 4)
       
       #Mark outliers where Residuals >= outlier_thresh_2 * data_range
        data$outliers[abs(data$gam_residuals) >= (outlier_thresh_2 * data_range)] <- TRUE
        }
    
    #Return the dataframe with a new boolean column "outliers"
     return(data[,-which(colnames(data)=="gam_residuals")])
    }
 
 
###############################################################################################################################################