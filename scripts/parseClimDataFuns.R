##########################################################
## Purpose: Functions to accompany parseCape.R
##
## Creator: Ian McGregor, March 2024
## Contact: 
## System: R Version 4.1.3
##########################################################

# --------------------------------------------------------------------#
## identifyFocalCells = use Feng raster to crop ERA5 data and return either the
##                      non-NA pixels, the Feng raster itself, or the cropped
##                      ERA5 raster depending on the function call
## prepData = define variable table and other vars for running the functions
## getCapeVals = get all cape values from entire timeframe to determine quantiles
## extractData = extract variable data from era5 netcdf and format
## calcVPD = calculate vpd using Fang et al 2022 equations (link in function)
## getTimes = create vector of timestamps from raster layers of netcdf
## parseNetcdf = main processing function to do first analysis pass of netcdf data,
##               including subsetting cape and vpd data to only what we need
##                at the month-year timesteps (or still everything for vpd). 
##                With this we can save smaller files as an intermediate step.
## calcVPDQuart = calculate driest quarter of each year using vpd
## analyzeVars = second analysis pass of data, takes output from parseNetcdf
##                and aggregates analyses across full timeframe.
## formatTifs = arrange and format variable rasters for plotting
## makePlots = create plots of the different climate variables over the basin
## savePlots = arrange the plots and save them to different pngs
# --------------------------------------------------------------------#
identifyFocalCells <- function(fl, polyClip, type){
  r <- rast(fl[1])[[1]]

  rCrop <- rast(polyClip)

  ## First crop the ERA5 data to the smaller extent
  ## Then resample the Feng data to be same res as ERA5
  ## Finally perform full mask
  r <- crop(r, rCrop)
  p <- resample(rCrop, r)
  r <- mask(r, p)
  
  if(loc=="amazon"){
    writeRaster(r, paste0(loc, "analysisRastTemplate.tif"))
  }

  if(type=="analysis") return(cells(r))
  if(type=="extract") return(r)
  if(type=="plot") return(r[[1]])
}
getBbox <- function(i, r){
  n <- as.vector(ext(r))[i]
  
  if(n < 0){
    if(i %in% c(1,3,4)) n <- floor(n)
    if(i == 2)  n <- ceiling(n)
  } else if(n > 0){
    if(i %in% c(1,3,4)) n <- ceiling(n)
    if(i == 2)  n <- floor(n)
  }
  
  return(n)
}
prepData <- function(loc, crd=NULL){
  climVars = data.table(vars=c("cape", "mean sea level pressure", 
                               "surface pressure", "2m dewpoint", "2m air temp",
                               "10m u-wind", "10m v-wind", "convective rain rate"),
                        varsAbb=c("CAPE", "MSL", "SP", "VAR_2D", "VAR_2T", 
                                  "VAR_10U", "VAR_10V", "CRR"),
                        fileLab=c("cape", "pressureMeanSea", "pressureSurface",
                                  "tempDew", "tempAir", "windU", "windV", 
                                  "convRainRate"))
  
  if(loc=="aus"){
    fBauman <- list.files("aus", pattern=".csv", full.names=TRUE)
    sites <- fread(fBauman[1], skip=2) #Bauman_et_al_2022_Nature_1984_2019_climatology_data
    
    return(list(climVars=climVars, fBauman=fBauman, sites=sites, spatPoly=NULL))
  }
  
  if(loc=="amazon"){
    ## bbox for amazon basin from feng paper
    f <- paste0("amazon/", list.files("amazon/", pattern="Datasets"))
    
    ## this is directly from the Feng data from the paper
    ## Datasets and codes used\current ERA5 CAPE\cape_1990_2019_ERA5_aft_27000m.tif
    f1 <- list.files(list.files(f, full.names=TRUE)[2], full.names=TRUE)
    r <- rast(f1)
    bboxNCAR <- sapply(1:4, getBbox, r)
    
    fl <- list.files("southAmerica/netcdf", full.names=TRUE)[1]

    if(is.null(crd)){
      sites <- identifyFocalCells(fl, polyClip=f1, type="analysis")
    } else {
      sites <- identifyFocalCells(fl, polyClip=f1, type="adriane")
    }

    return(list(climVars=climVars, spatPoly=f1, bboxNCAR=bboxNCAR, sites=sites))
  }

  if(loc=="southAmerica"){
    f <- list.files("southAmerica", pattern="metadata", full.names = TRUE)
    sites <- fread(f)
    
    return(list(climVars=climVars, sites=sites, spatPoly=NULL))
  }
}
getCapeVals <- function(fVar, nCores, sites, nTimeStart, polyClip){
  fData <- fVar[grepl("cape", fVar)]
  fData <- fData[nTimeStart:length(fData)]

  cl <- makeCluster(nCores)
  clusterEvalQ(cl, library(terra))
  clusterEvalQ(cl, library(data.table))
  clusterExport(cl, c("sites", "polyClip", "identifyFocalCells"))

  vals <- parSapply(cl, X=fData, function(d){
    r <- rast(d)
    rCrop <- identifyFocalCells(fl=fData, polyClip, type="extract")
    r <- crop(r, rCrop)
    vals <- data.table(values(r))
    vals <- vals[, ID := 1:nrow(vals)][sites, ]
    g <- melt(vals, id.vars="ID")
    return(g$value)
  })

  stopCluster(cl)
  return(vals)
}
extractData <- function(fileNum, fl, sites, loc, rCrop){
  r <- rast(fl[fileNum])
  
  if(loc=="amazon"){
    r <- crop(r, rCrop)
    vals <- data.table(values(r))
    vals <- vals[, ID := 1:nrow(vals)][sites, ]
  } else if(loc %in% c("aus", "southAmerica")){

    if(loc=="aus"){
      geo <- c("long", "lat")
    } else {
      geo <- c("Longitude.Decimal", "Latitude.Decimal")
    }

    v <- vect(sites, geom=geo, crs="EPSG:4326")

    ## ID is the row number of the original sites table
    vals <- data.table(extract(r, v))

    if(loc=="southAmerica") vals$ID <- sites$Plot.ID
  }

  ## reformat the table
  g <- melt(vals, id.vars="ID")
  return(list(site = g$ID, value=g$value))
}
calcVPD <- function(pressureSurface, tempAir, tempDew){
    ## first convert pressure Surface to hPa (from Pa)
    ## 1 hPa = 100 Pa
    pressureSurfaceHpa <- pressureSurface * 0.01

    ## then convert temps to C (by default are in K)
    # 1 deg C = 273.15 K
    tempAirC <- tempAir - 273.15
    tempDewC <- tempDew - 273.15

    ## the equations below come from here
    # https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022EF003019
    fw <- 1 + (7*(10e-4)) + (3.46 * (10e-6) * pressureSurfaceHpa)
    SVP <- 6.112 * fw * exp((17.67 * tempAirC) / (tempAirC + 243.5))
    AVP <- 6.112 * fw * exp((17.67 * tempDewC) / (tempDewC + 243.5))
    VPD <- SVP - AVP

    ## vpd is output as hPa from the equation. Convert to kPa for analyses
    VPD <- VPD * 0.1
    
    return(VPD)
}
findIndividualTimes <- function(X){
  k <- strsplit(X, "_")[[1]]
  return(k[[length(k)]])
}
getTimes <- function(fVar, climVars){
  ## no wind for now
  nFilesVar <- sapply(climVars$varsAbb[1:5], function(X){
    fVar <- fVar[grepl(X, fVar)]
    return(length(fVar))
  })

  nFilesVar <- nFilesVar[nFilesVar > 0]
  fVar <- fVar[grepl(paste0(names(nFilesVar), collapse="|"), fVar)]

  ## only keep the timestamps that are common to all the variables
  if(all(as.numeric(nFilesVar) == nFilesVar[1])){
    timestamps <- as.character(sapply(fVar, findIndividualTimes))
    timeSt <- unique(timestamps[which(duplicated(timestamps))])
  } else {
    stop(paste0("There are not equal numbers of files for the climate variables.",
              " Please double-check data downloads."))
  }
  
  return(timeSt)
}
parseNetcdf <- function(nTime, fVar, timeSeries, thresholds, sites, loc, locPath){
  fl <- fVar[grepl(paste0(timeSeries[nTime], collapse="|"), fVar)]
  fl <- fl[!grepl("10U|10V", fl)]
  varNames <- sapply(fl, function(q) return(strsplit(q, "\\.")[[1]][2]))
  varNames <- as.vector(varNames)

  if(loc=="amazon"){
    rCrop <- identifyFocalCells(fl, polyClip, type="extract")
  } else {
    rCrop <- NULL
  }
  varData <- lapply(1:length(fl), extractData, fl, sites, loc, rCrop)

  names(varData) <- varNames

  r <- rast(fl[1])
  timestampUTC <- time(r); rm(r)
  timestamp = format(timestampUTC, tz=tz)
  hours <- hour(timestamp)
  timeYear <- year(timestampUTC)[1] #each iteration is 1 month's data
  timeMonth <- month(timestampUTC)[1]

  ## cape feng
  cape <- data.table(site=varData$CAPE$site, cape=varData$CAPE$value)
  cape[, `:=` (timeHours = rep(hours, each=length(unique(site))),
                timeMonth = timeMonth, timeYear=timeYear)]

  capeFeng <- cape[timeHours >= 13 & hours <= 19][, timeHours := NULL]

  ## cape thresh
  capeThresh <- data.table(site=unique(cape$site), 
                      timeYear=timeYear, timeMonth=timeMonth)
  threshType <- c("Weak", "Moderate", "Strong")
  for(th in 1:length(thresholds)){
      thc <- cape[cape > thresholds[th], .N, by=site]
      missSites <- setdiff(capeThresh$site, thc$site)
      thc <- rbind(thc, data.table(site=missSites, N=0))
      capeThresh <- capeThresh[thc, on="site"]
      setnames(capeThresh, old="N", new=paste0("nHours", threshType[th]))
  }

  ## vpd bauman
  vpd <- calcVPD(pressureSurface = varData$SP$value, 
                tempAir = varData$VAR_2T$value,
                tempDew = varData$VAR_2D$value)
  vpdBauman <- cape[, `:=` (vpd=vpd, cape=NULL, timeHours=NULL)
              ][, mean(vpd), by=site
              ][, `:=` (timeYear = timeYear, timeMonth=timeMonth)]
  setnames(vpdBauman, old="V1", new="vpdMean")

  for(i in c("capeFeng", "capeThresh", "vpdBauman")){
    outFile <- paste0(locPath, i, "_", timeYear, "_", timeMonth)
    outFile <- paste0(outFile, ".Rdata")

    ## slower than fwrite, but much more efficient with memory
    if(i=="capeFeng") save(capeFeng, file=outFile)
    if(i=="capeThresh") save(capeThresh, file=outFile)
    if(i=="vpdBauman") save(vpdBauman, file=outFile)
  }

  return(paste0("Done with ", timeYear, "_", timeMonth))
}
calcVPDQuart <- function(s, allData, yrs){
  bl <- allData[site==s, ]

  ## Determine driest quarter of each year by identifying the 3
  ## consecutive months that had the highest mean VPD, then
  ## assign that value to the year
  rollM <- lapply(yrs, function(y){
    bly <- bl[timeYear==y, ]
    rollMean <- frollmean(bly$vpdMean, n=3)
    dryThresh <- max(rollMean, na.rm=TRUE)

    return(data.table(timeYear=y, vpdDryQuart=dryThresh))
  })
  rollM <- rbindlist(rollM)

  ## For an overall spatial analysis, we then take the mean of all
  ## the driest quarter VPDs across time so for each site we only have
  ## one VPD value
  out <- data.table(site=s, vpdDryMean_kPa = mean(rollM$vpdDryQuart))
  return(out)
}
analyzeVars <- function(allVars, locPath, nCores, yearStart, sites){
  ## format / process each variable individually
  for(targetVar in allVars){
    print(paste0("Analyzing ", targetVar))
    f <- list.files(locPath, pattern=targetVar, full.names=TRUE)
    if(targetVar=="capeFeng"){
      ## bring in afternoon hourly values from previous step
      cl <- makeCluster(nCores)
      clusterEvalQ(cl, library(data.table))
      allData <- parLapply(cl, X=f, function(q){
        load(q)
        return(capeFeng)
      })

      stopCluster(cl)
      allData <- rbindlist(allData)

      ## now calc mean for these values for each site for whole timeframe
      capeAft <- allData[timeYear >= yearStart, ][, mean(cape), by=site]
      setnames(capeAft, old="V1", new="capeAftMean")
      rm(allData)
    }
    if(targetVar=="capeThresh"){
      allData <- lapply(f, function(q){
        load(q)
        return(capeThresh)
      })
      allData <- rbindlist(allData)

      meanCapeThreshHours <- allData[timeYear >= yearStart, mean(nHoursWeak), by=site]

      meanCapeThreshHours <- meanCapeThreshHours[
          allData[timeYear >= yearStart, mean(nHoursModerate), by=site], on="site"
      ][allData[timeYear >= yearStart, mean(nHoursStrong), by=site], on="site"]

      colnames(meanCapeThreshHours) <- c("site", "capeHoursWeakM",
                                          "capeHoursModM", "capeHoursStrongM")

      rm(allData)
    }
    if(targetVar=="vpdBauman"){
      allData <- lapply(f, function(q){
        load(q)
        return(vpdBauman)
      })
      allData <- rbindlist(allData)
      allData <- allData[timeYear >= yearStart]
      yrs <- unique(allData$timeYear)

      cl <- makeCluster(nCores)
      clusterEvalQ(cl, library(data.table))
      clusterExport(cl, c("allData", "yrs"), envir=environment())

      if(grepl("southAmerica", locPath)) sites <- sites$Plot.ID

      meanVPDQuartYr <- rbindlist(parLapply(cl, X=sites, calcVPDQuart, allData, yrs))

      stopCluster(cl)
      rm(allData)
    }
    if(targetVar=="mcwd"){
      ## First, bring in mcwd across all years and get mean
      ### pattern= "carryover" or "reset"
      f <- list.files(paste0(loc, "/mcwd"), pattern="carryover", full.names = TRUE)
      r <- rast(f)
      rMean <- mean(r, na.rm=TRUE)
      
      ## Then, resample to match the other variables, extract the data, and
      ## add to table
      template <- rast(paste0(loc, "/analysisRastTemplate.tif"))
      rMeanRes <- resample(rMean, template, "cubic")
      valSite <- extract(rMeanRes, sites)
      meanMCWD <- data.table(site=sites, mcwdAnnualMean = valSite[,1])
    }
  }
  
  ## now bring everything together
  out <- meanCapeThreshHours
  out <- out[meanVPDQuartYr, on="site"][capeAft, on="site"][meanMCWD, on="site"][order(site)]
  
  return(out)
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
  varPlots <- varPlots[!grepl("Carry", varPlots)]
  plotVarList <- lapply(varPlots, formatTifs, r, dt)
  names(plotVarList) <- varPlots
  
  titles <- c(rep("Hours", 3), "VPD (kPa)", "CAPE (J/kg)", "MCWD (mm)")
  main <- c("Low CAPE threshold", "Moderate CAPE threshold", 
            "High CAPE threshold", 
            "Mean VPD in driest quarter", "Mean PM CAPE", "Mean Annual MCWD")

  plotList <- lapply(1:length(plotVarList), function(i){
    if(grepl("capeHours", varPlots[i])){
      cols <- rev(brewer.pal(9, "Blues"))
    } else if(grepl("vpd", varPlots[i])){
      ## lighter = drier, darker = wetter
      cols <- c("#663300", "#CC0000", "orange", "#FFCC00", "#FFFF33", "#FFFFCC")
    } else if(grepl("mcwd", varPlots[i])){
      ## lighter = drier, darker = wetter
      cols <- rev(c("#663300", "#993300", "#CC6633", "#CC9900", "#FFCC99"))
    }

    colPal <- colorRampPalette(cols)

    rf <- plotVarList[[i]]
    df <- as.data.frame(rf, xy=TRUE)
    colnames(df)[3] <- "value"
    g <- ggplot(df, aes(x=x, y=y, fill=value)) + 
      geom_raster() +
      xlab("") +
      ylab("") +
      scale_fill_gradientn(colours=colPal(10), name=titles[i]) +
      ggtitle(main[i]) +
      theme_bw() +
      theme(plot.margin = unit(c(0.1,0.2,0,0.1), 'lines')) #top, right, bottom, left

    return(g)
  })
  
  names(plotList) <- names(plotVarList)
  return(plotList)
}
savePlots <- function(figN, plotList, heatmapsVPD, filePre, figType){
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
      labels <- c("g", "h", "j", "k", "l")
    } else if(figType=="standalone"){
      labels <- c("A", "B", "C", "D")
    }
    
    # pl <- ggarrange(plotList$vpdDryMean_kPa,
    #                 plotList$mcwdAnnualMean,
    #                 ggplot() + theme_void(),
    #                 plotList$capeHoursWeakM,
    #                 plotList$capeHoursModM,
    #                 plotList$capeHoursStrongM,
    #                 ncol=3, nrow=2,
    #                 labels=labels,
    #                 label.x=0.03, label.y=1)
    # pl <- grid.arrange(plotList$vpdDryMean_kPa,
    #                   plotList$mcwdAnnualMean,
    #                   plotList$capeHoursWeakM,
    #                   plotList$capeHoursModM,
    #                   plotList$capeHoursStrongM,
    #                   layout_matrix=matrix(c(6,1,1,2,2,6,
    #                                          3,3,4,4,5,5), nrow=2, byrow=TRUE),
    #                   ncol=3, nrow=2)
                      
    pl <- ggarrange(
      ggarrange(plotList$vpdDryMean_kPa,
                plotList$mcwdAnnualMean,
                ncol=2, labels=labels[1:2],
                label.x=0.03, label.y=1),
      ggarrange(plotList$capeHoursWeakM,
                plotList$capeHoursModM,
                plotList$capeHoursStrongM,
                ncol=3, labels=labels[3:5],
                label.x=0.03, label.y=1),
      nrow=2
    )
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

  # CAPE mean afternoon, 1 panel
  if(figN == 3){
    fileOut <- paste0(filePre, "capeAft.png")
    widthF <- 14
    heightF <- 10
  }

  # CAPE thresholds + mean afternoon CAPE, 4-panel square
  if(figN == 4){
    fileOut <- paste0(filePre, "capeThresh_capeAft.png")
    widthF <- 22
    heightF <- 14

    pl <- ggarrange(plotList$capeHoursWeakM,
                    plotList$capeHoursModM,
                    plotList$capeHoursStrongM,
                    plotList$capeAftMean,
                    ncol=2, nrow=2, align="hv",
                    labels=c("A", "B", "C", "D"))
  }

  if(figN==5){
    fileOut <- paste0(filePre, "stormsPerspective3Pan.png")
    widthF <- 50
    heightF <- 12
    
    p1 <- heatmapsVPD$weak +
            ggtitle("") +
            xlab("Mean VPD of driest quarter (kPa)") +
            ylab("Mean hours per year above CAPE threshold") +
            theme(legend.position="inside", legend.position.inside = c(0.9, 0.8),
                  plot.margin=margin(0,1.5,0,0, unit="cm"))
    
    p2 <- plotList$vpdDryMean_kPa + 
            xlab("Degrees longitude") + 
            ylab("Degrees latitude")
            
    p3 <- plotList$capeHoursWeakM + 
            ggtitle("Hours above CAPE threshold") +
            xlab("Degrees longitude") + 
            ylab("Degrees latitude")
    
    p4 <- plotList$mcwdAnnualMean +
            xlab("Degrees longitude") + 
            ylab("Degrees latitude")
    
    w <- ggarrange(p1, p2, p3, p4,
                    ncol=4, nrow=1, labels=c("a", "b", "c", "d"),
                    label.x=0.03, label.y=1)
  }

  # add in annotations
  if(figN == 3){
    w <- plotList$capeAftMean +
          xlab("Degrees longitude") +
          ylab("Degrees latitude")
  } else if(figN < 5){
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
createStack <- function(varName, r, dt, rastFile){
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
    colors <- rev(brewer.pal(9, "Greys"))
    colLoess <- "#F0F0F0"
  }
  
  p <- ggplot(dtSub) +
    aes(x=dtSub[,get(idVar)], y=value) + 
    geom_bin2d(bins=40) +
    geom_smooth(method = "loess", color=colLoess, linewidth=2) +
    theme_classic() +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab(xlabs[i]) +
    ylab(ylabs[i]) +
    ggtitle(titles[i]) +
    theme(plot.margin = unit(c(0.3,0.2,0,0.3), 'lines'))
  
  if(grepl("vpd", idVar)){
    p <- p +
          scale_fill_gradientn(colors=colors, breaks=breaks, labels=breaksArea, 
                        name="1000s of km2")
  } else {
    p <- p + 
          scale_fill_gradientn(colors=colors, name="1000s of km2")
  }

  return(p)
}
