##########################################################
## Purpose: Process ERA5 climate data to look at cape and vpd metrics
##          over different areas
## Output:
##      - CAPE thresholds
##      - Rdata files with formatted data from netcdf files
## Creator: Ian McGregor
## Contact: mcgregori@caryinstitute.org
##########################################################
library(terra)
library(data.table)
library(parallel)
source("scripts/publication/parseClimDataFuns.R")

##-------------------------------------------------------------##
## DESCRIPTION OF ANALYSES
## for feng-adjacent analysis
### - Step 3. get table of all cape values in afternoon over entire timeframe
### - Step 4. calculate mean of those values

## for cape above thresholds, we need to
### - Step 3. calculate number hours above threshold for each pixel for each year
### - Step 4. calculate the mean hours across the full timeframe (each threshold)

## for vpd, we have to
### - Step 3. get mean vpd per pixel per month
### - Step 4. per pixel, determine what are the 3 driest consec months per year
###           using those mean vpd values
### - Step 5. Keep the highest (driest) vpd mean per year, then calculate the mean 
###           of those driest quarters (mean of means) across full timeframe

## for mcwd, we have to
###   calculate monthly cumulative water deficit per pixel, annually resetting
###     to 0 in the month that experienced the most precipitation for that pixel 
###   then get the annual max CWD
### - See script `MCWD_amazon.R` for the actual code. In this script, we simply
###     load the output MCWD raster.


##-------------------------------------------------------------##
## Step 0: Define key parameters
loc <- "amazon" #aus, amazon, or southAmerica
yearStart <- 1990 # earliest startYear = 1971.
yearEnd <- 2019 # for ERA5, 2019. For worldclim, 2000

##-------------------------------------------------------------##
## Step 1: define file paths and necessary vars
dataList <- prepData(loc)
climVars <- dataList$climVars
polyClip <- dataList$spatPoly
sites <- dataList$sites

locPath <- paste0(loc, "/climVars/")
if(!dir.exists(locPath)) dir.create(locPath)
fVar <- list.files(paste0(loc, "/netcdf"), full.names=TRUE)
timeSeries <- getTimes(fVar, climVars)

if(yearStart != 1971){
  nTimeStart <- min(which(grepl(as.character(yearStart), timeSeries)))
  timeSeries <- timeSeries[nTimeStart:length(timeSeries)]
} else {
  nTimeStart <- 1
}

nCores <- 55
tz <- "America/Manaus"
thresholds <- c(1023, 1900, 3000) #from Step 2 below
getCapeThresh <- FALSE

##-------------------------------------------------------------##
## Step 2: Calculate CAPE thresholds for amazon
## Only need to do once! thresholds are defined above
if(getCapeThresh){
  capeVals <- getCapeVals(fVar, nCores, sites, nTimeStart, polyClip)

  vals <- unlist(capeVals, recursive = FALSE, use.names = FALSE)
  range(vals); median(vals)
  quantile(vals, c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.995, 0.999))
}

## for Amazon 1971-2019, results are
### range = 0, 466656
### quantile(vals, c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.995, 0.999))
###     25%      50%      75%      90%      95%      99%    99.5%    99.9% 
### 126.375  492.000 1051.375 1671.875 2088.000 2933.500 3261.250 3974.000 

## For analysis, get 3 cape values
### the number of hours that cape was above each of the following
### - total > 1023 per year y axis, x-axis=mean VPD during those hours
### - total > 1900
### - total > 3000

##-------------------------------------------------------------##
## Step 3: function to parse netcdf files.
## Output = individual Rdata files, only need to run once, then can do
##          different analyses without having to redo the data
cl <- makeCluster(nCores)
clusterEvalQ(cl, library(terra))
clusterEvalQ(cl, library(data.table))
clusterExport(cl, c("fVar", "timeSeries", "sites", "parseNetcdf", "polyClip", "tz",
                      "calcVPD", "extractData", "locPath", "identifyFocalCells"))

parLapply(cl, X=1:length(timeSeries), parseNetcdf, fVar, timeSeries, thresholds, 
              sites, loc, locPath)

stopCluster(cl)

##-------------------------------------------------------------##
## Step 4: Individual metrics
allVars <- c("capeThresh", "vpdBauman", "mcwd")
stormData <- analyzeVars(allVars, locPath, nCores, yearStart, sites)
fwrite(stormData, paste0(loc, "/processedClimVars", yearStart, "_", yearEnd, ".csv"))
