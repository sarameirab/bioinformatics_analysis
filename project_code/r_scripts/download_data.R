RANDOM_FOREST_RESULTS_URL <- "https://drive.google.com/file/d/1n8hV1mfLbsqR-zABud4JX1LccxhVhRH9"

download_url <- function(url, type = "csv", export_type = TRUE) {
  
  # 1. Define the base directory
  base_dir <- "~/workspace/bioinformatics_analysis/data"
  
  # 2. Ensure the directory exists, create it if not
  if (!dir.exists(base_dir)) {
    message("Creating directory: ", base_dir)
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # 3. Create a hash from the URL
  file_hash <- digest(url, algo = "md5")
  
  # 4. Construct the filename with the hash and type (extension)
  file_name <- paste0(file_hash, ".", type)
  
  # 5. Construct the full output path
  output_file_path <- file.path(base_dir, file_name)
  
  # 6. Check if the file already exists
  if (file.exists(output_file_path)) {
    message("File already exists, skipping download: ", output_file_path)
  } else {
    message("File not found. Downloading to: ", output_file_path)
    
    # 7. Download the file (only if it doesn't exist)
    tryCatch({
      # Determine if 'type' should be used for export or just extension
      if (export_type) {
        googledrive::drive_download(
          file = googledrive::as_id(url), # Use as_id to robustly get the file ID
          path = output_file_path,
          type = type, # Use type for potential export (e.g., Sheets to CSV)
          overwrite = FALSE # Ensure we don't overwrite (though check handles this)
        )
      } else {
        googledrive::drive_download(
          file = googledrive::as_id(url),
          path = output_file_path,
          overwrite = FALSE
        )
      }
      message("Download complete.")
    }, error = function(e) {
      message("An error occurred during download: ", e$message)
      # Optionally, you might want to remove a potentially partially downloaded file
      if (file.exists(output_file_path)) {
        file.remove(output_file_path)
      }
      # Set output_file_path to NULL or NA to indicate failure
      output_file_path <<- NA # Using <<- to modify the outer scope variable
    })
  }
  
  # 8. Return the path of the file (or NA if download failed)
  return(output_file_path)
}