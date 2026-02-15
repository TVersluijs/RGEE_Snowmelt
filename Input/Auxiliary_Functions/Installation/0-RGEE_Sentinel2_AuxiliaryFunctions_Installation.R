
#Custom initialization of earth engine using reticulate based on code by Carl Higgs (https://github.com/r-spatial/rgee/issues/361).
#Code was expanded by adding the functionality of authorizing google drive and storing the credentials in the right location.
f_ee_Init <- function(user = NULL, project = NULL, fldr_asset = NULL, fldr_credentials_user = NULL, drive = TRUE, quiet = FALSE) {

  #Load packages
  library(rgee)
  library(reticulate)

  #Get credentials from environment if not provided
  if(is.null(user)){
    user <- Sys.getenv("EARTHENGINE_USER", unset = NA)
    if (is.na(user) || user == ""){stop("EARTHENGINE_USER not set in .Renviron")}
    }

  if(is.null(project)){
    project <- Sys.getenv("EARTHENGINE_PROJECT", unset = NA)
    if (is.na(project) || project == ""){stop("EARTHENGINE_PROJECT not set in .Renviron")}
    }

  #Get asset folder from environment (bypasses ee_Initialize()'s check)
  if(is.null(fldr_asset)){
    fldr_asset <- Sys.getenv("EARTHENGINE_ASSET", unset = NA)
    if(is.na(fldr_asset) || fldr_asset == ""){stop("EARTHENGINE_ASSET not set in .Renviron")}
    }

  #Set variables in environment
  Sys.setenv(EARTHENGINE_USER = user)
  Sys.setenv(EARTHENGINE_PROJECT = project)
  Sys.setenv(EARTHENGINE_ASSET = paste0("users/", fldr_asset))

  #Extract path to earthengine credentials
  fldr_credentials <- dirname(fldr_credentials_user)

  #Print information
  if(!quiet){
    message("Initializing Earth Engine:")
    message("  User: ", user)
    message("  Project: ", project)
    message("  Asset folder: ", paste0("users/", fldr_asset))
    }

  #Set fldr_asset in rgee's namespace (bypasses the broken ee_check_root_folder() check)
  #Use assignInNamespace to modify the locked namespace
  tryCatch({assignInNamespace("fldr_asset", fldr_asset, ns = "rgee")}, error = function(e) {
    # If that fails, just set a global option that we can check later
    options(rgee.ee_user = paste0("users/", fldr_asset))
    })

  #Google Drive authentication
  if(drive==TRUE){

    #(1): Generate GD credential in the ./earthengine folder
    repeat {

      #List all files in the credentials folder
      full_credentials <- list.files(path = fldr_credentials, full.names = TRUE)

      #Check if any files match the google drive credential name
      drive_condition <- grepl(".*_.*@.*", basename(full_credentials))

      #If the GDC does not exist in fldr_credentials then run the authentication
      if(!any(drive_condition)){

        #run googledrive::drive_auth
        suppressMessages(googledrive::drive_auth(email=FALSE, cache = fldr_credentials))
        if(!quiet){message("✓ Google Drive Credentials created ...")}

        }

      #If at least one GDC exists in fldr_credentials make sure only the most recent remains
      if(any(drive_condition)) {

        #Extract credentials file name
        drive_credentials <- full_credentials[drive_condition]

        #Run authentication for specific user email to renew credential
        suppressMessages(googledrive::drive_auth(email = user, cache = fldr_credentials))
        if(!quiet){message(paste0("✓ Google Drive Credentials saved at ", fldr_credentials))}

        #Avoid that multiple token files are found in fldr_credentials
        new_full_credentials <- list.files(path = fldr_credentials, full.names = TRUE)
        new_drive_condition <- grepl(".*_.*@.*", basename(new_full_credentials))
        if(sum(new_drive_condition) > 1){

          #Delete the oldest tokens if the system detects two different token files.
          tokens <- new_full_credentials[new_drive_condition]
          files_credentials_time <- file.info(tokens)$ctime
          drive_credential_to_remove <- tokens[files_credentials_time < max(files_credentials_time)]
          file.remove(drive_credential_to_remove)

        }
        break
      }

    }

    #(2): Copy GD credential to the user folder
    if (file.exists(fldr_credentials_user)) {

      #Remove previous GD credentials in fldr_credentials_user
      new_full_credentials_user <- list.files(path = fldr_credentials_user, full.names = TRUE)
      new_drive_condition_user <- grepl(".*_.*@.*", basename(new_full_credentials_user))
      if (sum(new_drive_condition_user) > 0) {
        tokens_user <- new_full_credentials_user[new_drive_condition_user]
        file.remove(tokens_user)
      }

      #Copy new GD credential to user folder
      file.copy(
        from = drive_credentials,
        to = sprintf("%s/%s", fldr_credentials_user, basename(drive_credentials)),
        overwrite = TRUE)

      #Print progress
      if(!quiet){message(paste0("✓ Google Drive Credentials copied to ", fldr_credentials_user))}

      #Store final path to GD credential
      path_GDC_user <- sprintf("%s/%s", fldr_credentials_user, basename(drive_credentials))

    }

  }
  if(drive==FALSE){path_GDC_user <- NA}

  #Create rgee_sessioninfo.txt if it doesn't exist (to prevent rgee bug)

    #Create directory if it doesn't exist
    if(!dir.exists(fldr_credentials)){dir.create(fldr_credentials, recursive = TRUE)}

    #Specify file path
    sessioninfo_file <- file.path(fldr_credentials, "rgee_sessioninfo.txt")

    #If file already exists first delete it
    if(file.exists(sessioninfo_file)){
      file.remove(sessioninfo_file)
    }

    #Create data.frame matching rgee's format:
    df <- data.frame(
      email = user,
      user = paste0("users/", fldr_asset),
      drive_cre = path_GDC_user,
      gcs_cre = NA,
      stringsAsFactors = FALSE
      )

    #Write using write.table like rgee does
    write.table(df, sessioninfo_file, row.names = FALSE)
    if(!quiet){message("✓ rgee_sessioninfo.txt file generated")}

  #Now call the Python Earth Engine initialization directly
  tryCatch({

    #Import earth engine
    ee <- reticulate::import("ee")

    #Authenticate using cached credentials
    ee$Initialize(project = project)

    #If no warning or error, print success
    if(!quiet){message("✓ Earth Engine initialized successfully")}

    invisible(TRUE)

    }, error = function(e) {stop("Failed to initialize Earth Engine:\n",
                                 "  ", e$message, "\n\n",
                                 "Try running: library(rgee); ee_Authenticate()")
                                 })


}

#Custom creation of user home asset folder. This function attempts to create a new asset home folder, but if a folder
#already exists it catches this error message, extracts the path of the already existing folder and renames fldr_asset
#to this existing folder.
f_CreateAssetHome <- function(fldr_asset = NULL){
  tryCatch({rgee::ee$data$createAssetHome(paste0("users/", fldr_asset))},
           error = function(e) {
             msg <- conditionMessage(e)
             if(grepl("User home folder already exists", msg)) {
                  path <- sub('.*User home folder already exists:\\s*"([^"]+)".*', "\\1", msg)
                  cat("Asset folder path already exists:", path, "\n")

                  #Change fldr_asset to the already existing asset folder name (both in local and global environment)
                  fldr_asset <- sub("^projects/earthengine-legacy/assets/users/", "", path)
                  assign("fldr_asset", sub("^projects/earthengine-legacy/assets/users/", "", path), envir = .GlobalEnv)
                  cat("  -fldr_asset changed to:", fldr_asset, "\n")

                  return(invisible(NULL))}
             stop(e)
           })
}
