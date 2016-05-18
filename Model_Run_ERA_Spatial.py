# -*- coding: utf-8 -*-
"""

Kudryavtsev Model

     By Kang Wang 
    March 29, 2016
            
Input:
    (1) Location:
        input_lat: Latitude
        input_lon: Longitude
        
    (2) Climate : 
        Ta  : Mean annual air temperature (C)
        Aa  : Amplitude of air temperature (C)
        Hsn : Winter-Averaged Snow Depth (m)
        Rsn : Snow Density (kg/m3)
        vwc : Volumetric Water Content (m3 / m3)
        
    (3) Vegetation:
        Hvgf: Height of vegetation in frozen period (m)
        Hvgt: Height of vegetation in thawed period (m)
        Dvf : Thermal diffusivity of vegetation in frozen period (m2 s)
        Dvt : Thermal diffusivity of vegetation in thawed period (m2 s)

Output:
        1) Mean annual temperature on the top of permafrost (C)
        2) Active Layer Thickness (m)
    
References:
    
    Anisimov, O. A., Shiklomanov, N. I., & Nelson, F. E. (1997). 
        Global warming and active-layer thickness: results from transient general circulation models. 
        Global and Planetary Change, 15(3), 61-77.
    Romanovsky, V. E., & Osterkamp, T. E. (1997). 
        Thawing of the active layer on the coastal plain of the Alaskan Arctic. 
        Permafrost and Periglacial processes, 8(1), 1-22.
    Sazonova, T. S., & Romanovsky, V. E. (2003). 
        A model for regional‐scale estimation of temporal and spatial variability of active layer thickness and mean annual ground temperatures. 
        Permafrost and Periglacial Processes, 14(2), 125-139.
    Sturm, M., Holmgren, J., König, M., & Morris, K. (1997). 
        The thermal conductivity of seasonal snow. Journal of Glaciology, 43(143), 26-41.
    Ling, F., & Zhang, T. (2004). 
        A numerical model for surface energy balance and thermal regime of the active layer and permafrost containing unfrozen water. 
        Cold Regions Science and Technology, 38(1), 1-15.
    Wieder, W.R., J. Boehnert, G.B. Bonan, and M. Langseth. (2014). 
        Regridded Harmonized World Soil Database v1.2. Data set. 
        Available on-line [http://daac.ornl.gov] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, USA.  http://dx.doi.org/10.3334/ORNLDAAC/1247  . 
        
"""

import sys
sys.path.insert(0, 'Src')
import function_20160515 as k_model
import numpy as np
#import matplotlib.pyplot as plt
from time import clock
import csv

lonname    = 'lon';
latname    = 'lat';

input_file = 'Parameters/T_CLAY.nc4'
varname    = 'T_CLAY';
[lat_grid, lon_grid, Clay_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

input_file = 'Parameters/T_SAND.nc4'
varname    = 'T_SAND';
[lat_grid, lon_grid, Sand_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

input_file = 'Parameters/T_SILT.nc4'
varname    = 'T_SILT';
[lat_grid, lon_grid, Silt_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

input_file = 'Parameters/T_OC.nc4'
varname    = 'T_OC';
[lat_grid, lon_grid, Peat_percent] = k_model.import_ncfile(input_file, lonname, latname, varname)

for i in range(1979,2016):
    
    year_label = i;
    
    print'Runing for {0:d}'.format(year_label)
      
    input_file = 'Input/'+str(year_label)+'.csv';
    
    # Input file
    
    start1=clock()
    
    lat_list    = k_model.read_site_location(input_file, col_needed=0)
    lon_list    = k_model.read_site_location(input_file, col_needed=1)
    Ta_list     = k_model.read_site_location(input_file, col_needed=2)
    Aa_list     = k_model.read_site_location(input_file, col_needed=3)
    Hsn_list    = k_model.read_site_location(input_file, col_needed=4)
    Rsn_list    = k_model.read_site_location(input_file, col_needed=5)
    vwc_list    = k_model.read_site_location(input_file, col_needed=6)
    Hvgf_list   = k_model.read_site_location(input_file, col_needed=7)
    Hvgt_list   = k_model.read_site_location(input_file, col_needed=8)
    Dvf_list    = k_model.read_site_location(input_file, col_needed=9)
    Dvt_list    = k_model.read_site_location(input_file, col_needed=10)
    
    finish1=clock()
    
    print'Reading: {0:0.1f} s'.format((finish1-start1))
    
    #obs_file = 'input/for_calibration.csv';
    
    #ALT_Obs = k_model.read_site_location(obs_file, col_needed=2);
    #ALT_Obs = k_model.read_site_location(input_file, col_needed=11)
    
    n_grid = len(lat_list); 
    #n_grid=1
    
    ALT_grid = np.ones(n_grid) * -999.99
    
    cal = 1
    
    if cal ==1:
        
        start2=clock()
    
        for i in range(n_grid):
                
            input_lat   = lat_list[i]
            input_lon   = lon_list[i]
            
            Ta = Ta_list[i]
            Aa = Aa_list[i]
            
            Hsn = Hsn_list[i]
            rho_sn = Rsn_list[i]
        
            vwc = vwc_list[i]
        #    vwc = 0.8;
            
            Hvgf = Hvgf_list[i];
            Hvgt = Hvgt_list[i];
        
            Dvf = Dvf_list[i] * 1.0#E-8
            Dvt = Dvt_list[i] * 1.0#E-8
            
            if ~np.isnan(Ta) & ~np.isnan(Aa) & ~np.isnan(Hsn) & ~np.isnan(rho_sn) & ~np.isnan(vwc):
                
            #=======Calculation Starting=======    
                
                # Extract soil texture from Grid Soil Database (Netcdf files)
                #       According to locations  
                
                [p_clay, p_sand, p_silt, p_peat] = k_model.Extract_Soil_Texture(input_lat, input_lon, 
                         lon_grid, lat_grid, 
                         Clay_percent, Sand_percent, Silt_percent, Peat_percent);
                                         
                p_peat = 0.0
                                         
                #============
                
                input_file2 = 'Parameters/Typical_Thermal_Parameters.csv';
                
                # Estimate_Soil_Heat_Capacity
                [Ct, Cf] = k_model.Estimate_Soil_Heat_Capacity(input_file2, 
                            p_clay, p_sand, p_silt, p_peat, vwc);
                
                # Estimate_Soil_Thermal_Conductivity
                [Kt, Kf] = k_model.Estimate_Soil_Thermal_Conductivity(input_file2, 
                            p_clay, p_sand, p_silt, p_peat, vwc);
                
                # Estimate_Snow_Heat_Capacity_and_Thermal_Conductivity...
                [Csn, Ksn] = k_model.Estimate_Snow_Thermal_Parameters(rho_sn);
                
                # Calculate period of a year, cold season and warm season...
                [tao,tao1,tao2] = k_model.Estimate_Cold_Warm_Season_Length(Ta, Aa);
                
                # Calculate latent heat ...
                L = k_model.Estimate_Latent_Heat(vwc);
                
                #===========
                
                # Calculate snow effect
                [Tvg, Avg] = k_model.Estimate_Snow_Effect(Ta, Aa, Hsn, Csn, Ksn, rho_sn, tao);
                
                # Calculate vegetation effect    
                [Tgs, Ags] = k_model.Estimate_Vegetation_Effect(Tvg, Avg, Hvgf, Hvgt, 
                                            Dvf, Dvt, tao1, tao2, tao);
                
                #   Calculating the temperature at the top of permafrost :        
                Tps_numerator = k_model.Calculate_Tps_numerator(Tgs, Ags, Kf, Kt);
            
                Tps = k_model.Calculate_Temperature_Top_PF(Tps_numerator, Kf, Kt, Cf, Ct);
                
                # Calculating the Active Layer Thickness :
                
                ALT_grid[i] = k_model.Calculate_ALT(Tps_numerator, Tps, Ags, Kf, Kt, Cf, Ct, L, tao);
        
        ##=======Calculation finished=======
        
        finish2=clock()    
        
        print'Calculating: {0:0.1f} s'.format((finish2-start2))
      
    final = np.vstack((np.vstack((lat_list,lon_list)), ALT_grid));
    final = np.transpose(final);
      
    
    
    csvfile = file('Output/'+str(year_label)+'_Results.csv', 'wb')
    writer = csv.writer(csvfile)
    
    writer.writerows(final)
    
    csvfile.close()