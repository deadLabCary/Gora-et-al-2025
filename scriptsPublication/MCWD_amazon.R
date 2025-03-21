library(terra)
library(sf)


#Create annual raster stacks from raw CHELSA data
#####

#Generate annual precipitation raster stacks from monthly precip rasters

# List all precip_files in the directory
precip_files <- list.precip_files(pattern = "precipitation_directory\\.2.1.tif$", full.names = TRUE)

precip_files


# Define the years
years <- unique(gsub(".*(19[0-9]{2}).*", "\\1", precip_files[grep("19[0-9]{2}", precip_files)]))
years2 <- unique(gsub(".*(20[0-9]{2}).*", "\\1", precip_files[grep("20[0-9]{2}", precip_files)])) # Extracts unique years

MCWD_years<-c(years, years2)

MCWD_years

#create list
precip_files_by_year <- list()

# Loop through each year and filter precip_files
for (year in MCWD_years) {
  precip_files_by_year[[year]] <- precip_files[grep(year, precip_files)]
}

# Print the list of precip_files categorized by year
print(precip_files_by_year)


# Loop through each year and filter precip_files
for (year in MCWD_years) {
  
  # Create a raster stack (collection of layers)
  raster_stack <- rast(precip_files_by_year[[year]])
  
  
  #write stack to file
  writeRaster(raster_stack, filename = paste("precip_stacks/", year, "precip_raster_stack.tif"),  overwrite=FALSE)
  
  
  
}


#Generate annual evapotranspiration raster stacks from monthly et rasters

# List all et_files in the directory
et_files <- list.et_files(pattern = "et_directory\\.2.1.tif$", full.names = TRUE)

et_files


# Define the years
years <- unique(gsub(".*(19[0-9]{2}).*", "\\1", et_files[grep("19[0-9]{2}", et_files)]))
years2 <- unique(gsub(".*(20[0-9]{2}).*", "\\1", et_files[grep("20[0-9]{2}", et_files)])) # Extracts unique years

MCWD_years<-c(years, years2)

MCWD_years

#create list
et_files_by_year <- list()

# Loop through each year and filter et_files
for (year in MCWD_years) {
  et_files_by_year[[year]] <- et_files[grep(year, et_files)]
}

# Print the list of et_files categorized by year
print(et_files_by_year)


# Loop through each year and filter et_files
for (year in MCWD_years) {
  
  # Create a raster stack (collection of layers)
  raster_stack <- rast(et_files_by_year[[year]])
  
  
  #write stack to file
  writeRaster(raster_stack, filename = paste("et_stacks/", year, "et_raster_stack.tif"),  overwrite=FALSE)
  
  
  
}



# Apply Amazon crop and mask, correct scale of precip data
#####

#ET

# List all files in the directory
et_files <- list.files(pattern = "et_stacks\\.tif$", full.names = TRUE)

#define amazon boundary
amazon<-vect("amazon_shapefile.shp")

#set output directory
output_path<-"amazon_et/"

#
for (file in et_files){ 
  a<-rast(file)
  et_raster_amazon<-crop(a, amazon)
  et_raster_amazon <-mask(et_raster_amazon, amazon)
  output<-file.path(output_path, paste0(basename(file), "et_stack.tif"))
  writeRaster(et_raster_amazon, output, overwrite = TRUE)
  
  
}



# Precip

# List all files in the directory
precip_files <- list.files(pattern = "precip_stacks\\.tif$", full.names = TRUE)
amazon<-vect("amazon.shp")



output_path<-"amazon_precip/"

for (file in precip_files){ 
  a<-rast(file)
  precip_raster_amazon<-crop(a, amazon)
  precip_raster_amazon <-mask(precip_raster_amazon, amazon)
  precip_raster_amazon<-precip_raster_amazon/100
  output<-file.path(output_path, paste0(basename(file), "_corrected.tif"))
  writeRaster(precip_raster_amazon, output, overwrite = TRUE)
  
  
}


# Create annual reset raster from precipitation data by identifying mean max precip per pixel across years
#####

generate_reset_mask <- function(precip_files) {
  precip_stacks <- lapply(precip_files, rast)  # Load precipitation stacks
  precip_stack <- do.call(c, precip_stacks)   # Combine into a single SpatRaster
  
  
  index <- rep(1:12, length.out = nlyr(precip_stack))
  avg_precip <- tapp(precip_stack, index, fun = mean)
  
  # Identify the highest precipitation month per pixel
  max_month <- which.max(avg_precip)
  
  # Create a 12-band binary raster stack
  reset_mask <- rast(avg_precip)
  for (m in 1:12) {
    reset_mask[[m]] <- ifel(max_month == m, 0, 1)
  }
  return(reset_mask)  # Each layer corresponds to a month (1 for reset month)
}


#list files
precip_files <- list.files("./amazon_precip/", full.names = TRUE)

#select 1990-2019
precip_files <- precip_files[grepl(paste0(1990:2019, collapse="|"), precip_files)]

#apply function to generate reset mask
reset_mask <- generate_reset_mask(precip_files)  

#write to file
writeRaster(reset_mask, filename = "reset_mask.tif", overwrite=TRUE)


#Calculate Monthly Water Deficit
#####


#read in directories for monthly raster stackes
dir1 <- "mcwd/precip_stacks_annual"   
dir2 <- "mcwd/et_stacks_annual"

# Folder to save results
output_dir <- "mcwd/CWD"           
if(!dir.exists(output_dir)) dir.create(output_dir)

#pull files and sort
sortList <- function(dirN){
  # List all raster stacks 
  f <- list.files(dirN, full.names = TRUE)
  
  # Ensure lists are sorted to match corresponding rasters
  
  years <- sapply(f, function(i){
    sub("^\\s*([0-9]{4}).*", "\\1", basename(i))
  })
  years <- as.numeric(as.vector(years))
  
  f <- f[order(years)]
  return(f)
}

fileList <- lapply(c(dir1, dir2), sortList)
precip_files <- fileList[[1]]
et_files <- fileList[[2]]



#loop for creating monthly WD raster stacks
for (i in seq_along(precip_files)) {
  # Load corresponding rasters
  r1 <- rast(precip_files[i])
  r2 <- rast(et_files[i])
  
  #get year
  year <- sub("^\\s*([0-9]{4}).*", "\\1", basename(precip_files[i]))
  
  # Perform subtraction
  r_result <- r1 - r2
  
  # Create output filename
  output_file <- file.path(output_dir, paste0(year, "_cwd_stack.tif"))
  
  # Save result
  writeRaster(r_result, output_file, overwrite = TRUE)
}



#Calculate Cumulative Water Deficit and Maximum Cumulative Water Deficit
#####

#list CWD files
deficit_files <- list.files("mcwd/CWD", full.names = TRUE)

#select 1990-2019
deficit_files <- deficit_files[grepl(paste0(1990:2019, collapse="|"), deficit_files)]




# Calculate CWD and MCWD

#initialize december carryover
dec_prev<-0

#define reset month
reset_mask<-rast("mcwd/reset_mask.tif")


#Function to Calculate CWD
compute_cwd <- function(deficit_stack) {
  n_layers <- nlyr(deficit_stack)
  cwd <- deficit_stack[[1]]# Initialize with first layer
  cwd<-cwd + dec_prev
  rclMat <- matrix(c(0, 999999, 0), ncol=1, nrow=3)
  test <- classify(cwd, rclMat, right=TRUE)
  cwd<-cwd*reset_mask[[1]]
  min_cwd <- cwd 
  reset_mask<-reset_mask
  
  for (i in 2:n_layers) {
    cwd <- cwd + deficit_stack[[i]]  # Accumulate deficit
    cwd[cwd > 0] <- 0 # Reset when balance turns positive
    min_cwd <- min(min_cwd, cwd)
    cwd<-cwd*reset_mask[[i]]
  }
  
  dec_prev<<-cwd
  return(min_cwd)
}

#function for computing annual MCWD
compute_mcwd <- function(deficit_stack, output_filename) {
  # Compute Maximum Cumulative Water Deficit
  mcwd_raster <- compute_cwd(deficit_stack)
  
  # Save output raster
  writeRaster(mcwd_raster, output_filename, overwrite=TRUE)
  
  return(mcwd_raster)
}

#assign directory for annual output
output_dir <- "mcwd/CWD/MCWD/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Loop through each deficit raster stack to get MCWD
for (i in seq_along(deficit_files)) {
  # Load deficit raster stack
  deficit_stack <- rast(deficit_files[i])
  
  #year
  year <- sub("^\\s*([0-9]{4}).*", "\\1", basename(deficit_files[i]))
  
  # Generate output filename
  output_filename <- file.path(output_dir, paste0(year, "_mcwd", ".tif"))
  
  # Compute MCWD and save output
  compute_mcwd(deficit_stack, output_filename)
  
  cat("Processed:", output_filename, "\n")
}

#create a raster stack from annual MCWD
annual_reset_files<-list.files("CWD/MCWD", full.names = TRUE)
annual_reset_stack<-rast(annual_reset_files)

#write MCWD raster stack
writeRaster(annual_reset_stack, filename = "mcwd/CWD/MCWD/MCWD_precip_reset_1990_2019.tif")

