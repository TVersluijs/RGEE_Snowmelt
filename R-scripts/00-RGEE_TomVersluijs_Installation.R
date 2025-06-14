######################################################################################################################################

# Script: Script with installation procedures to use RGEE
# Author: Tom Versluijs, adapted from Ricardo Dal'Agnol da Silva
# Date Created: 2025-06-14
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

#(1): Setup a GEE account

######################################################################################################################################

  ## you need a GEE account before being able to use RGEE.
  ## log in the https://code.earthengine.google.com/ and register for one.

######################################################################################################################################

#(2): Install python through the conda environment

######################################################################################################################################

  ## The conda environment is where the GEE Python API will be located and from where the the RGEE package uses it.

  ## First you need to install the Miniconda OUTSIDE of R
  ## Install Miniconda3 at https://docs.conda.io/en/latest/miniconda.html (keep all default settings, i.e. install for local user only)
  ## Open 'anaconda' in the command prompt (window button --> anaconda, you will see anaconda prompt)
  ## Then type in the commands below one-by-one (without the #) to install the rgee_py environment and packages:
  # conda create -n rgee_py pip python=3.8.18     #or "conda create -n rgee_py pip python" for the newest python version
  # conda activate rgee_py
  # pip install google-api-python-client
  # pip install earthengine-api==0.1.370            #or 'pip install earthengine-api' for the newest earthengine-api version
  # pip install numpy

  #Important: there is a bug in the newest versions of the earthengine-api that result in the error "credentials have expired" when
  #initializing Google Earth Engine using RGEE. To prevent this we have installed an older version of the earthengine-api (0.1.370)
  #using an older version of python (3.8.18) at line 44 (see 'https://github.com/r-spatial/rgee/issues/353#issuecomment-1983765552' for more details).

  ## Conda should now be installed, now lets get the path to the environment. Type inside anaconda:
  # conda env list

  ## copy the path for the 'rgee_py' environment and set it in the variable below inside R:
  ## note the use of double backslashes \\
  ## This below is where is located in MY computer, you have to use the 'conda env list' command to find where it is located on your
  ## computer.
  rgee_environment_dir = "C:\\Users\\USERNAME\\miniconda3\\envs\\rgee_py"
  saveRDS(object = rgee_environment_dir, file = paste0(here(), "/Input/rgee_environment_dir.rds"))

######################################################################################################################################

#(3): Setup your Google Cloud Project

######################################################################################################################################

  #Starting 13 November 2024, users will need to use a Cloud project to access Google Earth Engine. The use
  #of Google Earth Engine will remain free for non-commercial and research usage. However, is mandatory to
  #set up a Google Cloud Project and register this for non-commercial use. More information can be found here:
  #https://developers.google.com/earth-engine/guides/transition_to_cloud_projects

  #(I): Install the Google Cloud CLI:

    #Download the Google Cloud Installer CLI from this website: (https://cloud.google.com/sdk/docs/install)
    #and follow the steps below to install gcloud.

    #Make sure to remember the folder where you install GCloud because we need this at line 211 of this script.

      #-(1): Pick all the default install options (i.e. install for Single user only, and keep 'Bundled Python' and 'Cloud Tools for Powershell' checked).
      #-(2): After the installation has completed, the installer gives you the option to create Start Menu and Desktop shortcuts,
      #      start the Google Cloud SDK shell, and configure the Google Cloud CLI. Make sure that you leave the options to start
      #      the shell and configure your installation selected. This will initiate the gcloud configuration.
      #-(3): In the command prompt you will be asked if you want to sign in. Press 'Y' to confirm and make sure you sign in with
      #      the gmail account that you plan to use to access Google Earth Engine.
      #-(4): After signing in you are asked to pick a cloud project. You get several options of which you need to pick 'Create a new
      #      project'. Select this option by entering the corresponding number.
      #-(5): You are asked to enter a Project ID. This is the name of your project on the Google Cloud Platform. I suggest to
      #      pick something like: YOURNAME_rgee_snowmelt. Press enter to confirm.
      #-(6): You should now get a confirmation that Google cloud is configured and is ready to be used. You can close the cmd window.

  #(II): Authenticate gcloud from within the rgee python environment:

    #The aim of the Earth Engine authentication flows is to get a security "token" from your signed-in account which can be
    #stored to give your scripts permission to access your data. Follow the steps below to complete the authentication process.
    #For more information, see: https://developers.google.com/earth-engine/guides/auth.

      #-(1): Run 'python.exe' found in folder specified in the variable 'rgee_environment_dir' (see line 58). Make sure that you
      #      find and run the executable in this exact location (i.e. in the envs/rgee_py subfolder).
      #-(2): Run the following two commands in the python console. This will retrieve your Earth Engine Credentials.
      #      Note that this commonly results in the following error: 'Invalid value for [--scopes]: https://www.google
      #      apis.com/auth/cloud-platform scope is required but is not requested.'. If this is the case, then you can
      #      run step 3 below instead.
      #         import ee
      #         ee.Authenticate()
      #-(3): Only run the following code if step 2 above returned an error. If the error at step 2 pertains to a
      #      different missing scope, then you can add this scope yourself in the code below by separating multiple
      #      scopes by a comma: scopes=["scope1", "scope2"]
      #         import ee
      #         ee.Authenticate(auth_mode="gcloud", scopes=['https://www.googleapis.com/auth/cloud-platform'])
      #-(4): A browser window will open to an account selection page. Select the gmail account that you want to use
      #      for authentication. Browse the consent screen and indicate that you are willing to grant the requested
      #      scopes by clicking 'Allow'. You should then get a confirmation that the authentication was successful.
      #      On the webpage it will state: "You are now authenticated with the gcloud CLI!". You can close the webpage
      #      Back in the python executable, it should say "Successfully saved authorization token." These credentials
      #      have now been stored to a .json file that is on my computer located at: c:/Users/USERNAME/AppData/Roaming/
      #      gcloud/application_default_credentials.json.

  #(III): Check if gcloud is found in windows (run code in R). This step can be skipped when running on a MAC.
    system("gcloud --version")

  #(IV): Enable the Earth Engine API for your project

    #Go to https://console.cloud.google.com/apis/library/earthengine.googleapis.com. Login with your gmail account. Select
    #your gcloud project from the dropdown menu. Then click ENABLE. You might get a message stating that you need credentials
    #to use this API. There is no need for this as this is not required when accessing Google Earth Engine through RGEE.

  #(V): Register project for noncommercial usage

    #(1): Go to https://code.earthengine.google.com/register and login with your gmail account.
    #(2): Click 'Register a Noncommercial or Commercial Cloud project'.
    #(3): Select Unpaid usage if you intend to use Google Earth Engine for research purposes.
    #(4): Pick 'Academia & Research' from the dropdown menu and click Next.
    #(5): Choose en existing Google Cloud Project and select the project ID that you have specified under step I.5 above.
    #(6): Click 'Continue to summary'.
    #(7): Your project is now registered for non-commercial purposes. This concludes the setup of your Google Cloud Project.

######################################################################################################################################

#(4): Requirements for R

######################################################################################################################################

  ## R: version 4.4.1 (this is the newest version that I tested successfully)
  # Link: https://cran.r-project.org/bin/windows/base/

  ## RStudio: version at least 2024.09.0 build 375 (this is the version that I tested successfully)
  ## Older versions do not show the GEE images in Viewer correctly.
  # Link: https://www.rstudio.com/products/rstudio/download/

  ## (OPTIONAL) RTools: used to build some packages in R from source
  # Note that RTools is likely not required for RGEE to function so its installation can be skipped.
  # Link: https://cran.r-project.org/bin/windows/Rtools/
  ## after installing the Rtools, make sure to run this command line below inside RStudio:
  # writeLines('PATH="${RTOOLS42_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

######################################################################################################################################

#(5): Install R packages

######################################################################################################################################

  ## if you installed everything above, you can now install the packages inside R

  # install/load general packages used in the scripts
  p_load(raster,
         #rgdal, #package not supported anymore
         #rgeos, #package not supported anymore
         sp,
         sf,
         leaflet,
         mapview,
         caret)

  ##Manually install some other packages
  #utils::install.packages("devtools")
  #devtools:::install_github("gearslaboratory/gdalUtils")
  #Note that this might result in an error due to the dependency on rgdal. This is however not essential for RGEE to function.

  # now some more specific packages related to using rgee
  p_load(geojsonio, remotes, reticulate, devtools, googledrive, magick)
  p_load(rgee)
  #devtools:::install_github("r-spatial/rgee") #Note: install Github version only when default rgee package does not work
  remotes::install_github("r-earthengine/rgeeExtra")

  #IMPORTANT! Save your R-script at this point and restart R before proceeding with the code below.
  #Also try restarting if the installation did not finish properly and run the installation again after
  #restart. After re-starting make sure to re-run all lines of R-code up to this point.

######################################################################################################################################

#(6): set python version and initialize Google Earth Engine

######################################################################################################################################

  #Specify which python version is used by reticulate
  reticulate::use_python(rgee_environment_dir, required=T)

  #Specify which python version is used by rgee
  rgee::ee_install_set_pyenv(
    py_path = rgee_environment_dir,
    py_env = "rgee_py"
    )

  #Manually override python paths for reticulate and rgee
  Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
  Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)

  #Important: reboot R-Studio at this point (no need to rerun any code after restarting).

  #Specify installation folder ('google-cloud-sdk\\bin') of gcloud CLI manually:
  Sys.setenv(EARTHENGINE_GCLOUD = "C:\\Users\\USERNAME\\AppData\\Local\\Google\\Cloud SDK\\google-cloud-sdk\\bin")

  #Authorize access to Earth Engine using OAuth2
  rgee::ee_Authenticate(user="ENTER YOUR OWN GOOGLE EMAIL (@gmail.com)", auth_mode='notebook')

  #Initialize Earth engine
  #Your 'project' is the ID of your Google Cloud Project which can be found when logging in at https://console.cloud.google.com/cloud-resource-manager
  rgee::ee_Initialize(user="ENTER YOUR OWN GOOGLE EMAIL (@gmail.com)",
                      project='ENTER THE PROJECT ID OF YOUR GOOGLE CLOUD PROJECT HERE (suggested format at line 88 was: YOURNAME_rgee_snowmelt)',
                      drive = T)

  #When running this for the first time you will need to conduct one or two authentication steps:
  #(1): You might need to allow Tidyverse to access and manage files on your Google Drive. IMPORTANT: make sure to cross the
  #     box to give this permission before confirming.
  #(2): Notebook Authenticator:
  #  -Verify that the correct user account is listed. If not, press 'SWITCH ACCOUNT' to select another gmail account.
  #  -Select your previously created Google Cloud Project to use for Authentication. Make sure that the project is listed in
  #   the dropdown menu under 'Earth Engine enabled Cloud Projects'. If not, redo step '3.V' from this installation script.
  #  -Keep 'Use read-only scopes' UNchecked.
  #  -Click 'Generate Token'.
  #  -You will be shown an account selection page. Click the user account that you want to grant access to.
  #  -A warning page is presented, indicating that Google did not create the app (i.e. the code in the notebook). Click Continue to acknowledge.
  #  -A consent screen pops up asking if you are willing to grant the requested scopes. Click 'Select all' and press 'Continue'.
  #  -The following screen will show an authorization verification code that can be copied and pasted back into the R-console.
  #  -You should get the confirmation: "Successfully saved authorization token."

  #If after running ee_Initialize there are five green checkmarks at 'user', 'Google Drive credentials', 'Initializing Google Earth Engine'
  #'Earth Engine account' and 'Python Path' then the installation and initialization has completed successfully. Congratulations!!

  #As a final step we create an asset folder on your local machine for the temporary storage of output files.
  #Note: This step can be skipped if it results in the error that an asset folder already exists.
  rgee::ee$data$createAssetHome("users/fldr_assets")

######################################################################################################################################

#(7): Some final checks

######################################################################################################################################

  #Note: in case errors are encountered in the steps below it is advisable to restart R and rerun these steps.

  #Check folder where credentials have been stored (on my pc this is C:\Users\USERNAME\.config\earthengine)
  rgee::ee_get_earthengine_path()

  #Check summary of python configuration
  reticulate::py_config()

  #Checks python and earthengine-api installation (this might result in an error, but this is no problem as long as ee_Initialize is ok)
  rgee::ee_check_python()
  rgee::ee_check_python_packages() #NOTE can be ignored
  rgee::ee_check_credentials()
  #rgee::ee_check_gcloud()
  #rgee::ee_check_gcloud() often resulted in a gcloud error. This is due to the misspecification of the path in which is checked for gcloud in this
  #function. This error can thus be ignored and the check has been disabled.

  #(OPTIONAL): Check that an .Renviron file is generated in ~//Documents with the following information
   #EARTHENGINE_PYTHON="C:\Users\USERNAME\miniconda3\envs\rgee_py"
   #EARTHENGINE_ENV="rgee_py"

  #(OPTIONAL): Check that an 'environments.txt' file is generated in C:\Users\USERNAME\.conda with the following information
   #C:\Users\USERNAME\miniconda3\envs\rgee_py

  #(OPTIONAL): Make sure an earthengine folder is present at 'C:\Users\USERNAME\.config'


######################################################################################################################################
