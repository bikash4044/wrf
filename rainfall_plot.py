import glob
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import cartopy.feature as cfeature
from netCDF4 import Dataset as nc
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
from cartopy.mpl.gridliner import LatitudeFormatter, LongitudeFormatter
import numpy.ma as ma

shapefile_path = r"../west_central.shp"
# shapefile_path = r"../west_central.shp"
gdf = gpd.read_file(shapefile_path)
reader = shpreader.Reader(shapefile_path)
shape_feature = cfeature.ShapelyFeature(reader.geometries(), ccrs.PlateCarree(), edgecolor='white', facecolor='none')
year = os.getcwd().split('/')[-1]
if not os.path.exists('plots'):
    os.makedirs('plots')
files = glob.glob("jun/wrfout_d02*")
lon = nc(files[0]).variables["XLONG"][0]
lat = nc(files[0]).variables["XLAT"][0]
# Create a geometry mask function
def create_mask(lon, lat, shapefile_gdf):
    mask = np.zeros(lon.shape, dtype=bool)
    for i in range(lon.shape[0]):
        for j in range(lon.shape[1]):
            point = Point(lon[i, j], lat[i, j])
            # Check if the point is inside the shapefile geometry
            if shapefile_gdf.contains(point).any():
                mask[i, j] = True
    return mask
def create_mask2(lon, lat, shapefile_gdf):
    # Create an empty mask of the same shape as the lon/lat grid
    mask = np.zeros(lon.shape, dtype=bool)

    # Convert lon/lat into shapely Points and check if they are contained within the shapefile boundary
    points = [Point(x, y) for x, y in zip(lon.flatten(), lat.flatten())]
    mask = np.array([shapefile_gdf.contains(p).any() for p in points]).reshape(lon.shape)
    
    return mask
mask = create_mask(lon, lat, gdf)


files = glob.glob("jun/wrfout_d02*")
ds1 = nc(files[0])
ds2 = nc(files[-1])  
raint1 = (ds2.variables["RAINC"][0] + ds2.variables["RAINNC"][0]) - (ds1.variables["RAINC"][0] + ds1.variables["RAINNC"][0])
        
subset_data = raint1[mask]

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
pcm = ax.pcolormesh(lon, lat,raint1, cmap="jet", transform=ccrs.PlateCarree())
xticks=[70, 76, 82]
yticks=[14, 19, 24]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xticks,fontsize=10)
ax.set_yticklabels(yticks,fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=15)
cbar = plt.colorbar(pcm)
cbar.set_label("Rainfall(mm)")
ax.add_feature(shape_feature, linewidth=1)
pcm.set_clim(vmin=subset_data.min(), vmax=subset_data.max())
ax.set_extent([68,85,14,27])
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
plt.title(f"Rainfall on June {year}", fontsize="20")
plt.savefig(f'plots/rainfall{year}jun.png')



files = glob.glob("jul/wrfout_d02*")
ds1 = nc(files[0])
ds2 = nc(files[-1])  
raint2 = (ds2.variables["RAINC"][0] + ds2.variables["RAINNC"][0]) - (ds1.variables["RAINC"][0] + ds1.variables["RAINNC"][0])
        
subset_data = raint2[mask]

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
pcm = ax.pcolormesh(lon, lat,raint2, cmap="jet", transform=ccrs.PlateCarree())
xticks=[70, 76, 82]
yticks=[14, 19, 24]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xticks,fontsize=10)
ax.set_yticklabels(yticks,fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=15)
cbar = plt.colorbar(pcm)
cbar.set_label("Rainfall(mm)")
ax.add_feature(shape_feature, linewidth=1)
pcm.set_clim(vmin=subset_data.min(), vmax=subset_data.max())
ax.set_extent([68,85,14,27])
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
plt.title(f"Rainfall on July {year}", fontsize="20")
plt.savefig(f'plots/rainfall{year}jul.png')


files = glob.glob("aug/wrfout_d02*")
ds1 = nc(files[0])
ds2 = nc(files[-1])  
raint3 = (ds2.variables["RAINC"][0] + ds2.variables["RAINNC"][0]) - (ds1.variables["RAINC"][0] + ds1.variables["RAINNC"][0])
        
subset_data = raint3[mask]

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
pcm = ax.pcolormesh(lon, lat,raint3, cmap="jet", transform=ccrs.PlateCarree())
xticks=[70, 76, 82]
yticks=[14, 19, 24]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xticks,fontsize=10)
ax.set_yticklabels(yticks,fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=15)
cbar = plt.colorbar(pcm)
cbar.set_label("Rainfall(mm)")
ax.add_feature(shape_feature, linewidth=1)
pcm.set_clim(vmin=subset_data.min(), vmax=subset_data.max())
ax.set_extent([68,85,14,27])
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
plt.title(f"Rainfall on August {year}", fontsize="20")
plt.savefig(f'plots/rainfall{year}aug.png')


files = glob.glob("sep/wrfout_d02*")
ds1 = nc(files[0])
ds2 = nc(files[-1])  
raint4 = (ds2.variables["RAINC"][0] + ds2.variables["RAINNC"][0]) - (ds1.variables["RAINC"][0] + ds1.variables["RAINNC"][0])
        
subset_data = raint4[mask]

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
pcm = ax.pcolormesh(lon, lat,raint4, cmap="jet", transform=ccrs.PlateCarree())
xticks=[70, 76, 82]
yticks=[14, 19, 24]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xticks,fontsize=10)
ax.set_yticklabels(yticks,fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=15)
cbar = plt.colorbar(pcm)
cbar.set_label("Rainfall(mm)")
ax.add_feature(shape_feature, linewidth=1)
pcm.set_clim(vmin=subset_data.min(), vmax=subset_data.max())
ax.set_extent([68,85,14,27])
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
plt.title(f"Rainfall on September {year}", fontsize="20")
plt.savefig(f'plots/rainfall{year}sep.png')


ds1 = nc(files[0])
ds2 = nc(files[-1])  
raint = raint1+raint2+raint3+raint4
        
subset_data = raint[mask]

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.PlateCarree()})
pcm = ax.pcolormesh(lon, lat,raint, cmap="jet", transform=ccrs.PlateCarree())
xticks=[70, 76, 82]
yticks=[14, 19, 24]
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xticks,fontsize=10)
ax.set_yticklabels(yticks,fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=15)
cbar = plt.colorbar(pcm)
cbar.set_label("Rainfall(mm)")
ax.add_feature(shape_feature, linewidth=1)
pcm.set_clim(vmin=subset_data.min(), vmax=subset_data.max())
ax.set_extent([68,85,14,27])
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
plt.title(f"Rainfall on {year}", fontsize="20")
plt.savefig(f'plots/rainfall{year}tot.png')


