######################################################################################################################################

# Script: Script with installation procedures to use RGEE
# Author: Tom Versluijs adapted from Ricardo Dal'Agnol da Silva
# Date Created: 2024-03-116
# R version 4.3.2
# Adapted from the following video tutorial: https://www.youtube.com/watch?v=1-k6wNL2hlo&ab_channel=RicardoDalagnol

######################################################################################################################################

 # clean environment
  rm(list = ls()); gc()

 # general libraries
  #renv::restore() #revert to last version of R-packages used to successfully run this script (optional).
  utils::install.packages("pacman")
  library(pacman)
  p_load(here)

######################################################################################################################################

#(1): GEE account

######################################################################################################################################

  ## you need a GEE account before being able to use RGEE.
  ## log in the https://code.earthengine.google.com/ and register for one.

######################################################################################################################################

#(2): installing python through conda environment 

######################################################################################################################################

  ## the conda environment is where the GEE Python API will be located. The RGEE package uses it.

  ## first you need to install the Miniconda OUTSIDE of R
  ## install Miniconda3 at https://docs.conda.io/en/latest/miniconda.html (keep all default settings, i.e. install for local user only)
  ## open 'anaconda' in the command prompt (window button --> anaconda, you will see anaconda prompt)
  ## then type in the commands below one-by-one (without the #) to install the rgee_py environment and packages:
  # conda create -n rgee_py pip python     #or "conda create -n rgee_py pip python=3.11" for a specific version
  # conda activate rgee_py
  # pip install google-api-python-client
  # pip install earthengine-api==0.1.370            #or 'pip install earthengine-api' for the newest version
  # pip install numpy
  
  #Important: there is a bug in the newest versions of the earthengine-api that result in the error "credentials have expired" when
  #initializing google earth engine using RGEE. To prevent this we have installed an older version of the earthengine-api (0.1.370)
  #at line 44.

  ## ok conda should now be installed, now lets get the path to the environment, type inside anaconda:
  # conda env list

  ## copy the path to the rgee_py environment, you will need to set it in the variable below inside R:
  ## note the use of double backslashes \\ 
  ## this below is where is located in MY computer, you have to use the 'conda env list' command to find where it is located on yours
  rgee_environment_dir = "C:\\Users\\tomve\\miniconda3\\envs\\rgee_py"
  saveRDS(object = rgee_environment_dir, file = paste0(here(), "/Input/rgee_environment_dir.rds"))
  
######################################################################################################################################

#(3): Install the gcloud CLI

######################################################################################################################################

  # Install the gcloud CLI: follow the instructions of the web page link below. Pick all the default install options (i.e. 
  #install for local user only). https://cloud.google.com/sdk/docs/install

  # Can check that gcloud works within the rgee python environment : 
  #     Run python.exe found in folder 'rgee_environment_dir'
  #     In the python console :
  #         import ee
  #         ee.Authenticate()
  # 
  #     Credentials on my computer have been saved to C:/Users/tomve/AppData/Roaming/gcloud/application_default_credentials.json.
  #
  # RESTART THE R SESSION AFTER INSTALLING GCLOUD BEFORE RUNNING THE NEXT LINES

  #Check if gcloud is found in windows
  system("gcloud --version")
 
######################################################################################################################################

#(4): pre-requirements for R

######################################################################################################################################

  ## R: version at least 4.3.2 (this is the version that I tested so far and works)
  # Link: https://cran.r-project.org/bin/windows/base/

  ## RStudio: a recent version is recommended.
  ## Older versions do not show the GEE images in Viewer correctly.
  # Link: https://www.rstudio.com/products/rstudio/download/

  ## RTools: needed to build some packages in R from source
  # Link: https://cran.r-project.org/bin/windows/Rtools/
  ## after installing the Rtools, make sure to run this command line below inside RStudio:
  # writeLines('PATH="${RTOOLS42_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")


######################################################################################################################################

#(5): R packages

######################################################################################################################################

  ## if you installed everything above, you can now install the packages inside R

  # install/load general packages used in the scripts
  p_load(raster,
         rgdal,
         rgeos,
         sp,
         sf,
         leaflet,
         mapview,
         caret)

  ##Manually install some other packages
  #utils::install.packages("devtools")
  #devtools:::install_github("gearslaboratory/gdalUtils")

  # now some more specific packages related to using rgee
  p_load(geojsonio, remotes, reticulate, devtools, googledrive, gdalUtils)
  devtools:::install_github("r-spatial/rgee")

  ## IMPORTANT! Restart R before proceeding with the code below. Also try restarting if the installation
  ## did not finish properly and run the installation again after restart. In any case make sure to
  #rerun 'rgee_environment_dir' at line 53. 

######################################################################################################################################

#(6): set python

######################################################################################################################################

  #Specify which python version is used by reticulate
  reticulate::use_python(rgee_environment_dir, required=T)

  #Specify which python version is used by rgee
  rgee::ee_install_set_pyenv(
    py_path = rgee_environment_dir, # Change it for your own Python PATH
    py_env = "rgee_py" # Change it for your own Python ENV
    )

  #Manually override python paths for reticulate and rgee
  Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
  Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)
  
  #Specify installation folder ('google-cloud-sdk\\bin') of gcloud CLI manually:
  Sys.setenv(EARTHENGINE_GCLOUD = "C:\\Users\\tomve\\AppData\\Local\\Google\\Cloud SDK\\google-cloud-sdk\\bin")
  #old: "C:\\Program Files (x86)\\Google\\Cloud SDK\\google-cloud-sdk\\bin"
  
  # Initialize the Python Environment
  # to clean credentials: ee_clean_credentials(). Now no longer working
  # If everything fails try to run this code before ee_Initialize: ee_Authenticate(auth_mode='notebook')
  rgee::ee_Initialize(user="tom.versluijs@gmail.com", drive = T)

  ##The first time running ee_Initialize some text will appear about google drive credentials, and ask you to log in your GEE account.
  #If after running that there are four green checkmarks at 'user', 'Google Drive credentials', 'Initializing Google Earth Engine'
  #and 'Earth Engine account' then the installation and initialization has completed succesfully. Congratulations!

  #As a final step we create an asset folder on your local machine for the temporary storage of output files
  ee$data$createAssetHome("users/Tom_assets")
  #ee_get_assethome()
  
  # #If ee_initialize() results in an error that your credentials have expired, you should downgrade the earthengine API version to
  # #version 0.1.370 (https://github.com/r-spatial/rgee/issues/353#issuecomment-1983765552)
  # library(reticulate)
  # py_config() # see the name of your conda (python) environment, in my case "rgee_py" 
  # reticulate::py_install('earthengine-api==0.1.370', envname='rgee_py') 
  # 
  # # Check the installation of "earthengine-api" with 
  # pyl <- py_list_packages()
  # pyl[pyl$package == "earthengine-api", ]
  # 
  # # check python version with
  # py_run_string("import sys; print(sys.version)")
  
######################################################################################################################################

#(7): Some final checks

######################################################################################################################################
  
  #Check folder where credentials have been stored (C:\Users\tomve\.config\earthengine)
  rgee::ee_get_earthengine_path()

  #Check summary of python configuration
  reticulate::py_config()

  #Checks python and earthengine-api installation (might still result in an error, but as long as ee_Initialize is ok)
  rgee::ee_check_python()
  rgee::ee_check_python_packages()
  rgee::ee_check_credentials()
  rgee::ee_check_gcloud()
  #potential gcloud error might be due to misspecification of path in which is checked for gcloud

  #Make sure an .Renviron file is generated in ~//Documents with the following information (without hashtag)
   #PATH="${RTOOLS42_HOME}\usr\bin;${PATH}"
   #EARTHENGINE_INIT_MESSAGE="True"
   #EARTHENGINE_PYTHON="C:\Users\tomve\miniconda3\envs\rgee_py"
   #EARTHENGINE_ENV="rgee_py"
  
  #Make sure a environments.txt file is generated in C:\Users\tomve\.conda with the following information (without hashtag)
   #C:\Users\tomve\miniconda3\envs\rgee_py
  
  #Make sure an earthengine folder is present at 'C:\Users\tomve\.config' containing the same files as found in the folder
  #'earthengine_backup_working_credentials' also found at "C:\Users\tomve\.config"
  

######################################################################################################################################
