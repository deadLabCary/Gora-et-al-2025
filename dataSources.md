# Primary data sources
- ERA5
  - Interpolated climate data for the whole world on a grid of 0.25 degrees (27.8 km), based on a mix of observational data and modeled forecasts. For more details, please see the [main website](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=overview).
- NCAR archive
  - This is a Research Data Archive provided by NSF and NCAR, with [this site](https://rda.ucar.edu/datasets/ds633.0/) showing the ERA5 data (listed in the "Data access" tab as "ERA5 atmospheric surface analysis"). For a bounding box of the Amazon basin, we used the shapefile of the basin provided in the supplementary code of [Feng et al (2023)](https://www.nature.com/articles/s41467-022-35570-1#-data-availability-).
- CHELSA
  - Reformatted ERA5 climate data downsampled to 1km resolution, mainly focused on precipitation and adjacent variables. Current version as of our analysis (2.1) is available [here](https://envicloud.wsl.ch/#/).
