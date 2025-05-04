# Description of files
Below is a description of the files found in [this google drive folder](https://drive.google.com/drive/u/0/folders/1VGbKbd7jWngFbh7VBVefhp-tQCF0sMF4) used to create Figures XXX for Gora et al., 2025.

## Folders
- amazon
  - climVars
    - `.Rdata` files for each month-year combination from 1971-2019. This is an output from `scripts/parseClimData.R`, which parses the netcdf files and isolates the target variables of CAPE (cape thresholds) and VPD.
  - Datasets and codes used
    - This is directly from the supplemental data for [Feng et al 2023](https://www.nature.com/articles/s41467-022-35570-1).
  - figures
    - figures for the paper
  - mcwd
    -Evapotranspiration (et_stacks_annual) and precipitation (precip_stacks_annual) rasters stacked annually for calculating Maximum Cumulative Water Deficit. Evapotranspiration (pet_penman) and precipitation (pr) data were downloaded from Chelsa (see below). MCWD was reset annually for each raster cell via reset_mask.tif after reaching the month the focal cell experiences the greatest water deficit. Annual calculated values for MCWD are aggregated across our study period in MCWD_precip_reset_2009_2019.tif. 
  - netcdf
    - Climate variables downloaded from Chelsa and ERA5 via the NCAR repository (for more details on the data sources, please see [this markdown file](https://github.com/deadLabCary/Gora-et-al-2025/blob/main/dataSources.md). Please see below for specifics about each variable.
  - `analysisRastTemplate.tif`
    - raster template that we use as the base gridded structure to visualize the data. This is created using the extent of the Feng et al., 2023 Amazon shapefile and the pixel resolution of the ERA5 data.
  - `processedClimVars1990_2019.tif`
    - aggregated raster with each pixel containing one value for each climate variable, averaged across the timeframe of 1990-2019. This is an output from `scripts/parseClimData.R`.
  - `processedClimVars1990.csv`
    - aggregated csv with one value of each climate variable for each pixel in the study area, averaged across the timeframe of 1990-2019. This is an output from `scripts/parseClimData.R`.
 
## Climate Variables
**NB**: Note that ERA5 data is agnostic of an above- / below-canopy designation. Officially, the relevant variables (temp, wind, surface pressure) are all at "Earth's surface" or Xm above "Earth's surface". 
- Convective Available Potential Energy (CAPE)
  - Hourly data, units are in J/kg, see metadata [here](https://codes.ecmwf.int/grib/param-db/59).
- Vapor Pressure Deficit (VPD)
  - VPD (or, how dry is it?) is not natively supplied by ERA5. Instead, we calculated VPD using the equations in [Fang et al (2022)](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022EF003019). Because air pressure (*Pmst*) is available from ERA5 as "surface pressure", we only used Equations 1-5 and therefore only needed additional data for temperature from ERA5. VPD was calculated in units of hPa first, then converted to kPa (hPa \* 0.1) for presentation.
    - Air and dewpoint temperature
      - Hourly data, units are in K and therefore converted to C before calculating VPD
      - See [here](https://codes.ecmwf.int/grib/param-db/167) for air temperature and [here](https://codes.ecmwf.int/grib/param-db/168) for dewpoint metadata.
    - Surface pressure
      - Hourly data, units are in Pa. This was converted to hPa (Pa \* 0.01) for the Fang equations.
      - View the metadata [here](https://codes.ecmwf.int/grib/param-db/134).
- Wind
  - We thought we needed wind at first and downloaded the 10m U and V windspeed components, but we ended up not using them. These are described [here (U)](https://codes.ecmwf.int/grib/param-db/165) and [here (V)](https://codes.ecmwf.int/grib/param-db/166).

### Specific methods of calculation
- CAPE thresholds
  - For our study period of 1990-2019, we calculated the total number of hours above a certain threshold CAPE value per year for every pixel in the Amazon basin, then took the mean of that to yield one value per pixel covering the whole timeframe. Thresholds for weak, moderate, and strong CAPE values were guided by taking quantiles of all pixels' CAPE values over the full timeframe, whereby we used rough quantiles of ~75% (weak, 1023 J/kg), ~94% (moderate, 1900 J/kg), and ~99% (strong, 3000 J/kg).
- CAPE afternoon
  - Following the methods of [Feng et al 2023](https://www.nature.com/articles/s41467-022-35570-1), we calculated the mean CAPE value per pixel over all afternoon hours (defined as 1-7pm) over the entire timeframe of 1990-2019.
- VPD mean driest quarter
  - The calculation of VPD was done as described above from Fang et al 2022. To get the mean VPD of the driest quarter, we did the following:
    - for each pixel, we identified the 3 consecutive months that had the highest mean VPD, then assigned that VPD value for the whole year. We do that for every year in the timeframe, then took the mean across years to have a single VPD value per pixel.
