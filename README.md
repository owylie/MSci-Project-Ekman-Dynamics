# MSci-Project-Ekman-Dynamics

These are the notebooks and other code files used for my masters project.
They calculate Ekman pumping and Ekman currents on monthly timescales from 2011-2020, and on daily timescales in 2020. 
The daily calculations are done for 2020 in a specific region of the arctic, but the area and year can be changed if one downloads different data and changes the code a bit.
All the files in this repository should be downloaded in the same folder, and as well as this the following data sets need to be downloaded (they're too large to upload here), and put in separate folders.

1) CPOM Geostrophic current data. http://www.cpom.ucl.ac.uk/dynamic_topography/. Request access, and download the file 'Full_DOT_data_Arco.nc'. 
   This file should go in a folder called CPOM_geo.

ncinfo on Full DOT:

    CPOM CryoSat2 Arctic Oceanography data 2020. Written by H Heorton
    
    Geoid: GOCO03s
    
    Surface_elevation: CPOM threshold
    
    Smoothing: 100 km guassian
    
    Time_written: 2021-10-06
    
    publisher_name: UCL_CPOM
    
    publisher_type: institution
    
    publisher_email: a.muir@ucl.ac.uk
    
    publisher_url: http://www.cpom.ucl.ac.uk/dynamic_topography/
    
    Time_dimension: Days since 2000-01-01
    
    dimensions(sizes): time(120), x(334), y(334)
    
    variables(dimensions): float32 Sea_level_anomaly(time, x, y), float32 time(time), float32 DOT_smoothed(time, x, y), float32 DOT_unsmoothed(time, x, y), float32 DOT_uncertainty_estimate(time, x, y), float32 Geo_surf_current_x(time, x, y), float32 Geo_surf_current_y(time, x, y), float32 lons(x, y), float32 lats(x, y), float32 xdist(x, y), float32 ydist(x, y), float32 ang_c(x, y), float32 ang_s(x, y)


2) ERA5 wind data. This can be downloaded as an hourly or monthly data set. To work in accordance with these notebooks, it will have to be dowloaded in the following way.
   
   Monthly: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form. 
     Select 'monthly averaged reanalysis', then '10m u-component of wind' and '10m v-component of wind', then the years 2011 to 2020, and all months, at hour 00:00 (the only option)
     For effiency, I only downloaded from 60 degrees to 90 degrees North, but a larger spatial region shouldn't affect the code. 
     This will download one file with all the data, which should be named 'monthly_2011-21.nc', and put in a folder called 'ERA5'.
   
   Hourly: https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form.
     Select the same 'reanalysis' and 10m wind components. I used 2020, all months, all hours, from 60 to 90 degrees North, but other years can be looked at by changing the code.
     This will download one file with all the data, which should be named '2020_all.nc' (or whatever year you chose), and put in a folder called 'ERA5'.

ncinfo on the file:

    <class 'netCDF4._netCDF4.Dataset'>

    root group (NETCDF3_64BIT_OFFSET data model, file format NETCDF3):

    Conventions: CF-1.6
    
    history: 2021-12-16 15:38:32 GMT 
    
    by grib_to_netcdf-2.23.0: /opt/ecmwf/mars-client/bin/grib_to_netcdf -S param -o /cache/data7/adaptor.mars.internal-1639669111.9061267-31108-1-4e62d4da-6097-40db-8e36-d4c777aeb562.nc /cache/tmp/4e62d4da-6097-40db-8e36-d4c777aeb562-adaptor.mars.internal-1639669111.0741432-31108-1-tmp.grib
    
    dimensions(sizes): longitude(1440), latitude(721), time(2)
    
    variables(dimensions): float32 longitude(longitude), float32 latitude(latitude), int32 time(time), int16 u10(time, latitude, longitude), int16 v10(time, latitude, longitude)
    
    
3) Pathfinder Ice Drift. https://nsidc.org/data/NSIDC-0116/versions/4. You will need to make an account to access this data.
   This should be downloaded in yearly .nc files (dont change the names), and put in a folder called 'Pathfinder'.

ncinfo on one of the yearly ice drift files:

    <class 'netCDF4._netCDF4.Dataset'>
    
    root group (NETCDF4 data model, file format HDF5):

    version: 4.1
    
    release_date: Apr 2021
    
    Conventions: CF-1.4
    
    dataset_doi: 10.5067/INAWUWO7QH7B
    
    dimensions(sizes): x(361), y(361), time(366)
    
    variables(dimensions): float64 x(x), float64 y(y), float64 time(time), int32 crs(), float32 u(time, y, x), float32 v(time, y, x), int16 icemotion_error_estimate(time, y, x), float32 latitude(y, x), float32 longitude(y, x)
    
    

4) NSIDC Ice Concentration. https://nsidc.org/data/NSIDC-0051/versions/1. You will have to manually download the data for each year (it takes a loooong time)
   In the 'filter by date' box, put in the year you want (will need to do at least 2011 to 2020). In 'filter spatially', range from 60 to 90 degrees North.
   Save each year of data files in a folder labelled by the year (eg. '2015'), and put all these folders in a folder called 'NSIDC_nt'.



The files FINAL_MONTHLY.ipynb calculates the monthly averaged Ekman pumping and currents over the entire timescale, while IceGovernor.ipynb compares the different variables. IceGovernor.ipynb focuses only on a specific region of the Arctic (around the Beaufort Gyre area) in 2020, but this can be modified reasonably simply to cover any region or year.
Ice_Drift.py, Ice_Conc.py, Calculating_pumping.py, Ekman_pumping.py and grid_set.py are scripts that contain classes referenced in the two notebooks. 

Ekman_pupming.py and grid_set.py were created by Harry Heorton, and Calculating_pumping.py is heavily influenced by other code written by Harry.

NSIDC_gs.npz, Pathfinder.npz and PS_20km_gs2021.npz are grids created by Harry Heorton, that are also used in the notebooks.

I appologise for the unnecessarily long (and lazy) code in places - please feel free to adapt it and add in a few loops to make it more efficient!!

Finally, a huge thank you to Michel Tsamados and Harry Heorton for their help on this project.
