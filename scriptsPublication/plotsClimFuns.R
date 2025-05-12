##########################################################
## Purpose: Functions to accompany plotsClim.R
##
## Creator: Ian McGregor
## Contact: mcgregori@caryinstitute.org
##########################################################

# --------------------------------------------------------------------#
## createStack = make a stack of the target variable rasters
## makePlotsHeatmap = make heatmaps of CAPE hours (different thresholds) vs 
##                    different combinations of VPD and MCWD
## formatTifs = arrange and format variable rasters for plotting
## makePlots = create plots of the different climate variables over the basin
## savePlots = arrange the plots and save them to different pngs
# --------------------------------------------------------------------#
createStack <- function(varName, r, dt, maskForest, rastFile){
  rasts <- lapply(varName, function(X){
    r1 <- r
    r1[dt$site] <- dt[, get(X)]
    names(r1) <- X
    
    return(r1)
  })
  
  w <- rast(rasts)
  
  writeRaster(w, rastFile, overwrite=TRUE)
}
makePlotsHeatmap <- function(i, dtSub, allVars, figType, idVar, xlabs, ylabs, titles){
  dtSub <- dtSub[variable == allVars[i]]

  if(i==1) breaks <- seq(0, 150, length.out=6)[2:6]
  if(i==2) breaks <- seq(0, 100, length.out=6)[2:6]
  if(i==3) breaks <- seq(0, 160, length.out=6)[2:6]
  if(i==4) breaks <- seq(0, 150, length.out=6)[2:6]

  ## ERA5 data is 27.8 km res (0.25 degrees)
  areaPix <- 27.8*27.8
  breaksArea <- round((breaks*areaPix)/1000)

  if(figType=="standalone"){
    colors <- c("#663300", "#CC0000", "orange", "#FFCC00",
                "#FFFF33", "#FFFFCC")
    colLoess <- "white"
  } else if(figType=="supplement"){
    colors <- brewer.pal(9, "Greens")[2:9]
    # colPal <- colorRampPalette(colors)
    # colLoess <- "#F0F0F0"
    colLoess <- "#333333"
  }
  
  if(grepl("vpd", idVar) & grepl("mcwd", allVars[i])){
    p <- ggplot(dtSub) +
      aes(y=dtSub[,get(idVar)], x=value) + 
      geom_bin2d(bins=40)
  } else {
    p <- ggplot(dtSub) +
      aes(x=dtSub[,get(idVar)], y=value) + 
      geom_bin2d(bins=40)
  }
  
  p <- p +
    geom_smooth(method = "loess", color=colLoess, linewidth=1) +
    theme_classic() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab(xlabs[i]) +
    ylab(ylabs[i]) +
    ggtitle(titles[i]) +
    theme(plot.margin = unit(c(0.3,0.2,0,0.3), 'lines'))
  
  if(grepl("vpd", idVar) & !grepl("mcwd", allVars[i])){
    p <- p +
          scale_fill_gradientn(colors=colors, breaks=breaks, labels=breaksArea, 
                        name="1000s of km2")
  } else {
    p <- p + 
        scale_fill_gradientn(colors=colors, name="1000s of km2") +
        scale_x_reverse()
  }
  return(p)
}
formatTifs <- function(var, r, dt){
  cellsAll <- 1:ncell(r)
  naCells <- setdiff(1:ncell(r), dt$site)

  cellsAll[dt$site] <- dt[, get(var)]
  cellsAll[naCells] <- NA
  values(r) <- cellsAll

  return(r)
}
makePlots <- function(r, dt){
  varPlots <- colnames(dt)
  varPlots <- varPlots[!grepl("site", varPlots)]
  plotVarList <- lapply(varPlots, formatTifs, r, dt)
  names(plotVarList) <- varPlots
  
  titles <- c(rep("Hours", 3), "VPD (kPa)", "MCWD (mm)")
  main <- c("Low CAPE threshold", "Moderate CAPE threshold", 
            "High CAPE threshold", 
            "Mean VPD in driest quarter", "Mean Annual MCWD")

  plotList <- lapply(1:length(plotVarList), function(i){
    if(grepl("cape", varPlots[i])){
      cols <- c(brewer.pal(9, "Blues")[3:9], "black")
      colHist <- cols[length(cols)/2]
      binwidth <- 25
    } else if(grepl("vpd", varPlots[i])){
      ## lighter = wetter, darker = drier
      # cols <- rev(c("#663300", "#CC0000", "orange", "#FFCC00", "#FFFF33", "#FFFFCC"))
      cols <- brewer.pal(9, "YlOrRd")
      colHist <- "red"
      binwidth <- 0.1
    } else if(grepl("mcwd", varPlots[i])){
      ## lighter = wetter, darker = drier
      # cols <- c("black", "#663300", "#7b532b", "#996633", "#CC9933", 
      #               "#CC9900", "#FFCC33", "#fde18d", "#FFFFCC")
      cols <- c("black", rev(brewer.pal(9, "YlOrBr")[c(1, 3, 5, 7:9)]))
      colHist <- "#683904fb"
      binwidth <- 50
    }

    colPal <- colorRampPalette(cols)

    rf <- plotVarList[[i]]
    df <- as.data.frame(rf, xy=TRUE)
    colnames(df)[3] <- "value"
    
    ## histograms
    ## total area for each pixel is 28.07809 km x 28.07809 km = 788.379138 km2
    gHist <- ggplot(df, aes(x=value)) + 
      geom_histogram(aes(y=after_stat(count)*0.788379138), position="identity",
                      color="black", fill=colHist, alpha=0.5, binwidth=binwidth,
                      closed="left", boundary=0)+
      xlab(titles[i]) +
      ylab("Area (1000s of km)") +
      ggtitle(main[i]) +
      # ylim(0, 1000) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.01)),
                          limits=c(0,1100), 
                          breaks=seq(0,1000,250),
                          labels=seq(0, 1000, 250)) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank())
    
    ## mcwd x axis is reversed in `savePlots()`, so we fix the axis there
    if(!grepl("mcwd", varPlots[i])){
      gHist <- gHist +
                scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
    }
    
    if(i %in% c(1,6)){
      dft <- as.data.table(df)
      ## quantiles for paper
      if(i==1) dft[value >= 300, .N] / nrow(dft) #total forest area above hours at cape threshold
      if(i==6) dft[value <= -500, .N] / nrow(dft) #total forest area below MCWD threshold
    }
    
    ## maps
    g <- ggplot(df, aes(x=x, y=y, fill=value)) + 
      geom_raster() +
      xlab("") +
      ylab("") +
      ggtitle(main[i]) +
      theme_bw() +
      theme(plot.margin = unit(c(0.1,0.2,0,0.1), 'lines'),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) #top, right, bottom, left
    
    if(grepl("mcwd", varPlots[i])){
      g <- g +
        scale_fill_gradientn(colours=cols, name=titles[i], 
                              breaks=c(0, -250, -500, -750, -1000),
                              labels=c(0, -250, -500, -750, -1000))
                  
    } else {
      g <- g +
        scale_fill_gradientn(colours=colPal(10), name=titles[i])
    }
    

    return(list(g=g, gHist=gHist))
  })
  
  names(plotList) <- names(plotVarList)
  return(plotList)
}
savePlots <- function(figN, plotList, plotListHist, heatmapsVPD, heatmapsMCWD, 
                      heatmapsBoth, filePre, figType){
  res <- 350
  units <- "cm"
  labelX <- -78
  labelY <- 7

  # VPD + 3 CAPE thresholds, 4-panel square
  if(figN==1){
    fileOut <- paste0(filePre, "vpd_mcwd_capeThresh.png")
    widthF <- 22
    heightF <- 14
    
    if(figType=="supplement"){
      labels <- c("", "g", "h", "", "j", "k", "l")
    } else if(figType=="standalone"){
      labels <- c("A", "B", "C", "D")
    }
    
    topRow <- ggarrange(
                NULL, 
                plotList$vpdDryMean_kPa,
                plotList$mcwdAnnualMean,
                NULL,
                nrow=1, widths = c(0.35, 1, 1, 0.35),
                labels=labels[1:4],
                label.x=0.03, label.y=1)
    bottomRow <- ggarrange(
                plotList$capeHoursWeakM,
                plotList$capeHoursModM,
                plotList$capeHoursStrongM,
                ncol=3, labels=labels[5:7],
                label.x=0.03, label.y=1)
    
    pl <- ggarrange(topRow, bottomRow, nrow=2)
  }

  # CAPE thresholds, 3-panel vertical
  if(figN==2){
    fileOut <- paste0(filePre, "capeThresh.png")
    widthF <- 11
    heightF <- 21

    pl <- ggarrange(plotList$capeHoursWeakM,
                    plotList$capeHoursModM,
                    plotList$capeHoursStrongM,
                    ncol=1, nrow=3,
                    labels=c("A", "B", "C"))
  }

  if(figN==3){
    fileOut <- paste0(filePre, "stormsPerspective9Pan.png")
    widthF <- 50
    heightF <- 36
    
    pHeatVPD <- heatmapsVPD$weak +
            ggtitle("") +
            xlab("Mean VPD of driest quarter (kPa)") +
            ylab("Mean hours per year above CAPE threshold") +
            theme(legend.position="inside", legend.position.inside = c(0.9, 0.8),
                  plot.margin=margin(0,1.5,0.5,0, unit="cm"))
    
    pHeatMCWD <- heatmapsMCWD$weak +
            ggtitle("") +
            xlab("Mean annual MCWD (mm)") +
            ylab("Mean hours per year above CAPE threshold") +
            theme(legend.position="inside", legend.position.inside = c(0.9, 0.8),
                  plot.margin=margin(0,1.5,0.5,0, unit="cm"))
    
    pHeatBoth <- heatmapsBoth$vpdmcwd + 
            ggtitle("") +
            xlab("Mean VPD of driest quarter (kPa)") +
            ylab("Mean annual MCWD (mm)") +
            theme(legend.position="inside", legend.position.inside = c(0.15, 0.8),
                  plot.margin=margin(0,1.5,0.5,0, unit="cm"))
    
    pMapCape <- plotList$capeHoursWeakM + 
            ggtitle("Hours above CAPE threshold") +
            xlab("Degrees longitude") + 
            ylab("Degrees latitude")
            
    pMapVPD <- plotList$vpdDryMean_kPa + 
            xlab("Degrees longitude") + 
            ylab("Degrees latitude")
    
    pMapMCWD <- plotList$mcwdAnnualMean +
            xlab("Degrees longitude") + 
            ylab("Degrees latitude")
    
    pHistCape <- plotListHist$capeHoursWeakM
    pHistVPD <- plotListHist$vpdDryMean_kPa
    
    pHistMCWD <- plotListHist$mcwdAnnualMean +
                  scale_x_reverse(expand = expansion(mult = c(0, 0.05)),
                                    breaks=seq(-1200, 0, 200),
                                    labels=seq(-1200, 0, 200))
    
    w <- ggarrange(pHeatVPD, pHeatMCWD, pHeatBoth, 
                    pMapCape, pMapVPD, pMapMCWD,
                    pHistCape, pHistVPD, pHistMCWD,
                    ncol=3, nrow=3, 
                    labels=c("a", "b", "c", "d", "e", "f", "g", "h", "i"),
                    label.x=0.03, label.y=1)
  }

  # add in annotations
  if(figN < 3){
    w <- annotate_figure(pl, 
            left = textGrob("Degrees Latitude", rot = 90, gp = gpar(cex = 1.2)),
            bottom = textGrob("Degrees longitude", gp = gpar(cex = 1.2)))
  }

  if(figType=="supplement"){
    return(w)
  } else if(figType=="standalone"){
    png(fileOut, res=res, units=units, width=widthF, height=heightF)
    print(w)
    dev.off()
  }
}
