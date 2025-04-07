##########################################################
## Purpose: Plot the processed ERA5 climate data to look at metrics
##          over different areas
## Output:
##      - summary plots and figures for publication
## Creator: Ian McGregor
## Contact: mcgregori@caryinstitute.org
##########################################################
library(data.table)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(terra)
library(RColorBrewer)

source("scriptsPublication/plotsClimFuns.R")

dt <- fread("amazon/processedClimVars1990.csv")
loc <- "amazon"

#-------------#
# For creating the 7-panel figure for storms perspective piece, need to run
## - Step 1 = individual plots
## - Step 2 = heatmap of vpd vs cape metrics
## - Step 3 = combine

#-------------#
# Plots Step 1: Create the individual rasters for each variable
r <- rast(paste0(loc, "/analysisRastTemplate.tif"))

## Pixel location in amazon basin
# values(r)[!is.na(values(r))] <- cells(r)
# plot(r, col=viridis::viridis(75), main="Site (pixel) number distribution")

# make raster stack of the data
varName <- colnames(dt)[2:ncol(dt)]
createStack(varName, r, dt, rastFile="amazon/processedClimVars1990_2019.tif")

#-------------#
# Plots Step 2: vpd (or mcwd) vs cape metrics heatmap (single and multiplot)
heatmaps <- lapply(c("vpd", "mcwd", "both"), function(v){
  if(v=="vpd"){
    idVar <- "vpdDryMean_kPa"
    xlab <- "Mean VPD of Driest Quarter (kPa)"
  } else if(v=="mcwd"){
    idVar <- "mcwdAnnualMean"
    xlab <- "Mean Annual MCWD (mm)"
  } else {
    idVar <- "vpdDryMean_kPa"
    allVars <- "mcwdAnnualMean"
    ylabs <- "Mean Annual MCWD (mm)"
    xlabs <- "Mean VPD of Driest Quarter (kPa)"
  }
  
  dtSub <- melt(dt, id.vars=c("site", idVar))
  
  if(v!="both"){
    ## binned count distribution map
    allVars <- c("capeHoursWeakM", "capeHoursModM", "capeHoursStrongM", "capeAftMean")
    ylabs <- c("", "", "", "Mean afternoon CAPE")
    xlabs <- c("", "", "", xlab)
    titles <- c("Low CAPE threshold", "Moderate CAPE threshold", 
                "High CAPE threshold", "")
  } else {
    ylabs <- ""
    xlabs <- ""
    titles <- ""
  }
  
  # fig = "standalone" (red-or-yel) or "supplement" (greyscale)
  plotsHeatmap <- lapply(1:length(allVars), makePlotsHeatmap, dtSub, allVars, 
                  figType="supplement", idVar, ylabs, xlabs, titles) 
  
  if(v=="both"){
    names(plotsHeatmap) <- "vpdmcwd"
  } else {
    names(plotsHeatmap) <- c("weak", "mod", "strong", "aft")
  }
  return(plotsHeatmap)
})
heatmapsVPD <- heatmaps[[1]]
heatmapsMCWD <- heatmaps[[2]]
heatmapsBoth <- heatmaps[[3]]

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
                                      gp = gpar(cex = 1)))

png("amazon/figures/capeThresholds.png", width=32, height=24, res=350, units="cm")
print(w)
dev.off()

#-------------#
# Plots Step 3: Create / format individual maps for multi-panel figure
## make maps for cape thresholds, vpd, and mcwd
filePre <- "amazon/figures/"
plotListAll <- makePlots(r, dt)
plotList <- lapply(plotListAll, '[[', 1)
plotListHist <- lapply(plotListAll, '[[', 2)

#-------------#
## Make multiple figures using `makePlots()` function
# fig1 = VPD + 3 CAPE thresholds, 4-panel square
# fig2 = CAPE thresholds, 3-panel vertical
# fig3 = CAPE mean afternoon, 1 panel
# fig4 = CAPE thresholds + mean afternoon CAPE, 4-panel square

for(fig in 1:2){
  savePlots(figN=fig, plotList, plotListHist, heatmapsVPD, heatmapsMCWD, heatmapsBoth, 
          filePre, figType="standalone")
}

#-------------#
# 9-panel for paper (main figure in manuscript)
## figType argument doesn't matter for this
savePlots(figN=3, plotList, plotListHist, heatmapsVPD, heatmapsMCWD, heatmapsBoth, 
          filePre, figType="standalone")

#-------------#
# 11-panel for paper (SI)
## Step 1: Create `plotList` from the maps
plotPanel5 <- savePlots(figN=1, plotList, plotListHist, heatmapsVPD, heatmapsMCWD, 
                        heatmapsBoth, filePre, figType="supplement")
                        
## Step 2: Combine the heatmaps plot with the 5 panel plot
png("amazon/figures/stormsPerspective11Pan.png", res=350, width=46, height=16,
    units="cm")
grid.arrange(w,                                    
             plotPanel5,                               
             ncol = 5, nrow = 1, 
             layout_matrix = rbind(c(1,1, 2, 2,2), c(1,1, 2, 2, 2)))
dev.off()