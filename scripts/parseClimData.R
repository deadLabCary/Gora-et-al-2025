##########################################################
## Purpose: Process ERA5 climate data to look at cape and vpd metrics
##          over different areas
##
## Creator: Ian McGregor, March 2024
## Contact: 
## System: R Version 4.1.3
##########################################################
library(terra)
library(data.table)
library(parallel)
source("parseClimDataFuns.R")

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


##-------------------------------------------------------------##
## Step 0: Define key parameters
loc <- "amazon" #aus, amazon, or southAmerica
yearStart <- 1979 # earliest startYear = 1971. Max timeframe is 2019
yearEnd <- 2019 # for ERA5, 2019. For worldclim, 2000

##-------------------------------------------------------------##
## Step 1: define file paths and necessary vars
dataList <- prepData(loc)
climVars <- dataList$climVars
polyClip <- dataList$spatPoly
sites <- dataList$sites

## there is a duplicated plot ID from Adriane's data. For now, I am adding a
## 0.1 to the end of one of them (so e.g. we will have 2080.1 and 2080.2).
if(loc=="southAmerica"){
  nr <- which(sites$Plot.ID==2080.1)
  sites[nr[2], Plot.ID := 2080.2]
}

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
tz <- ifelse(loc=="amazon|southAmerica", "America/Manaus", "Australia/Brisbane")
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

## For analysis, get 4 cape values
### the number of hours that cape was above each of the following
### - mean afternoon cape > 1023 (Feng paper)
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
allVars <- c("capeFeng", "capeThresh", "vpdBauman", "mcwd")
stormData <- analyzeVars(allVars, locPath, nCores, yearStart, sites)

# for adriane's sites, need to bring in the coordinates of each plot
if(loc=="southAmerica"){
  setkeyv(stormData, "site")
  setkeyv(sites, "Plot.ID")
  stormData <- stormData[sites]
  stormData <- stormData[order(site)]
}
fwrite(stormData, paste0(loc, "/processedClimVars", yearStart, "_", yearEnd, ".csv"))

##-------------------------------------------------------------##
# Step 5: PLOTS
library(data.table)
library(ggplot2)
library(ggpubr)
library(grid)

dt <- fread("amazon/processedClimVars1990.csv")

#-------------#
# For creating the 7-panel figure for storms perspective piece, need to run
## - Step 1 = individual plots
## - Step 2 = heatmap of vpd vs cape metrics
## - Step 3 = combine

#-------------#
# Plots Step 1: Create the individual rasters for each variable
library(terra)
library(RColorBrewer)
fl <- list.files("amazon/netcdf", full.names=TRUE)
r <- rast(paste0(loc, "/analysisRastTemplate.tif"))

## Pixel location in amazon basin
# values(r)[!is.na(values(r))] <- cells(r)
# plot(r, col=viridis::viridis(75), main="Site (pixel) number distribution")

# make raster stack of the data
varName <- colnames(dt)[2:ncol(dt)]
createStack(varName, r, dt, rastFile="amazon/processedClimVars1990_2019.tif")

#-------------#
# Plots Step 2: vpd (or mcwd) vs cape metrics heatmap (single and multiplot)
dt <- fread("amazon/processedClimVars1990.csv")

heatmaps <- lapply(c("vpd", "mcwd"), function(v){
  if(v=="vpd"){
    idVar <- "vpdDryMean_kPa"
    xlab <- "Mean VPD of Driest Quarter (kPa)"
  } else if(v=="mcwd"){
    idVar <- "mcwdAnnualMeanReset"
    xlab <- "Mean Annual MCWD (mm)"
  }
  
  dtSub <- melt(dt, id.vars=c("site", idVar))

  ## binned count distribution map
  allVars <- c("capeHoursWeakM", "capeHoursModM", "capeHoursStrongM", "capeAftMean")
  ylabs <- c("", "", "", "Mean afternoon CAPE")
  xlabs <- c("", "", "", xlab)
  titles <- c("Low CAPE threshold", "Moderate CAPE threshold", 
              "High CAPE threshold", "")

  # fig = "standalone" (red-or-yel) or "supplement" (greyscale)
  plotsHeatmap <- lapply(1:length(allVars), makePlotsHeatmap, dtSub, allVars, 
                  figType="supplement", idVar, ylabs, xlabs, titles) 

  names(plotsHeatmap) <- c("weak", "mod", "strong", "aft")
  return(plotsHeatmap)
})
heatmapsVPD <- heatmaps[[1]]
heatmapsMCWD <- heatmaps[[2]]

## CAPE thresholds figure
pl <- ggarrange(heatmapsVPD$weak, heatmapsMCWD$weak,
                heatmapsVPD$mod, heatmapsMCWD$mod,
                heatmapsVPD$strong, heatmapsMCWD$strong,
                nrow=3, ncol=2, labels=c("a", "d", "b", "e", "c", "f"), 
                label.x=0.03)
w <- annotate_figure(pl, left = textGrob("Time above CAPE threshold (mean hours per year)", 
                                    rot = 90, gp = gpar(cex = 1.2)),
                    bottom = textGrob(c("Mean VPD of Driest Quarter (kPa)",
                                        "Mean annual MCWD (mm)"),
                                      x=c(0.22, 0.72),
                                      y=c(0.54, 0.54),
                                      gp = gpar(cex = 1.2)))

png("amazon/figures/capeThresholds.png", width=32, height=24, res=350, units="cm")
print(w)
dev.off()

#-------------#
# Plots Step 3: Create / format individual maps for multi-panel figure
## make maps for cape thresholds, vpd, and mcwd
filePre <- "amazon/figures/"
plotList <- makePlots(r, dt)
names(plotList)[6] <- "mcwdAnnualMean"

#-------------#
## Make multiple figures using `makePlots()` function
# fig1 = VPD + 3 CAPE thresholds, 4-panel square
# fig2 = CAPE thresholds, 3-panel vertical
# fig3 = CAPE mean afternoon, 1 panel
# fig4 = CAPE thresholds + mean afternoon CAPE, 4-panel square

for(fig in 1:4){
  savePlots(figN=fig, plotList, plotsHeatmap, filePre, figType="standalone")
}

#-------------#
# 3-panel for paper
## figType argument doesn't matter for this
savePlots(figN=5, plotList, heatmapsVPD, filePre, figType="standalone")

#-------------#
# 7-panel for paper

## Step 1: Create `plotList` from the maps
plotPanel4 <- savePlots(figN=1, plotList, plotsHeatmap, filePre, figType="supplement")
## Step 2: Create `pl` from vpd vs cape metrics with fig="supplement" and add
##          the annotation with the figure text
## Step 3: Combine here
library(gridExtra)

png("amazon/figures/stormsPerspective7Pan.png", res=350, width=46, height=16,
    units="cm")
grid.arrange(w,                                    
             plotPanel4,                               
             ncol = 4, nrow = 1, 
             layout_matrix = rbind(c(1,1, 2, 2), c(1,1, 2, 2)))
dev.off()
