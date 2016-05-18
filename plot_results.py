# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:10:52 2016

@author: kang
"""

import sys
sys.path.insert(0, 'Src')
import function_20160515 as k_model
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap

ncfile = 'input/ERA_Reanalysis_NH.nc'
lonname = 'longitude'
latname = 'latitude'

# Read the nc file
    
fh = Dataset(ncfile, mode='r')
    
# Get the lat and lon
#   Set the grid size for lat. and lon. (here is 0.5 degree)
    
Lon_list = fh.variables[lonname][:]; 
Lat_list = fh.variables[latname][:]; 

Lon_list = np.mod((Lon_list + 180.),360.) - 180.;

# make 2-d grid of lons, lats
lons, lats = np.meshgrid(Lon_list,Lat_list)

from matplotlib.collections import LineCollection

import shapefile

r = shapefile.Reader("IPA_Map/Untitled")
shapes = r.shapes()

for iii in range(1979,2016):
    
    year =  iii;
    
    filename = 'Output/'+str(year)+'_Results.csv';
    
    results = np.loadtxt(filename, delimiter=',');
    
    results = np.reshape(results[:,2],[116,720])
    
    results[np.where(results==-999.99)] = np.nan;
    
    fig,ax = plt.subplots(figsize=(6,6))
    
    my_map = Basemap(projection='npstere',boundinglat=28,lon_0=0,
                     resolution='l')
                     
    my_map.drawcoastlines()
    #my_map.drawcountries()
    plt.title('Forced by ERA for '+str(year))
    
    # compute native x,y coordinates of grid.
    x, y = my_map(lons, lats)
    
    # set desired contour levels.
    clevs = np.array([0.0, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0])
    
    CS2 = my_map.contourf(x,y,results, clevs, extend='max')
    
    cbar = my_map.colorbar(CS2)
    cbar.set_label('Active Layer Thickness (m)', rotation=90)
    
    #from matplotlib import path as mpath
    
    for shape in shapes:
        lons1,lats1 = zip(*shape.points)
        data = np.array(my_map(lons1, lats1)).T
     
        if len(shape.parts) == 1:
            segs = [data,]
        else:
            segs = []
            for i in range(1,len(shape.parts)):
                index = shape.parts[i-1]
                index2 = shape.parts[i]
                segs.append(data[index:index2])
            segs.append(data[index2:])
     
        lines = LineCollection(segs,antialiaseds=(1,))
    #    lines.set_facecolors(cm.jet(np.random.rand(1)))
        lines.set_edgecolors('k')
        lines.set_linewidth(1.2)
        ax.add_collection(lines)
        
    fig.savefig('Figure/'+str(year)+'.png', 
                format='png', 
                bbox_inches='tight')