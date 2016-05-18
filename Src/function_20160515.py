def read_site_location(input_file, col_needed):

    """ 
    The function is to read latitude and longitude as input;
    INPUTs:
            input_file: a file contains information of some sites.
            col_needed: Column index of needed

    OUTPUTs:
            data_selected: data selected

    DEPENDENTs:
            None 
    """
    
    import numpy as np
    
    data = np.loadtxt(input_file,skiprows=1,delimiter=',')
        
    data_selected = data[:,col_needed]
    
    return data_selected
    
def Estimate_Snow_Thermal_Parameters(rho_sn):
    
    """ 
    The function is to estimate snow thermal conductivity and heat capacity ;
    INPUTs:
            rho_sn: snow density (kg m-3). 

    OUTPUTs:
            Csn: heat capacity of snow  (J kg-1 C-1)
            Ksn: thermal conductivity of snow  (W m-1 C-1)

    DEPENDENTs:
            None 
    """
    
    # Conductivity of snow: 
    #   eq-4, Sturm et al., 1997:    
    Ksn = (rho_sn/1000.)*(rho_sn/1000.)*3.233-1.01*(rho_sn/1000.)+0.138; # Unit: (W m-1 C-1)
    
    # Capacity of snow:
    #   eq-30, Ling et al., 2004; OR Table-1, Goodrich, 1982.
    Csn = 2090.0;                                                # Unit: J m-3 C-1 
    
    return Csn, Ksn

def Estimate_Soil_Heat_Capacity(input_file, p_clay, p_sand, p_silt, p_peat, vwc):
    
    """ 
    The function is to estimate soil heat capacity ;
    INPUTs:
            input_file: file contains typical thermal parameters of clay, sand, silt and peat.
            p_clay: percent of clay (%)
            p_sand: percent of sand (%)
            p_silt: percent of silt (%)
            p_peat: percent of peat (%)
            vwc   : volumetric water content (m3 / m3)

    OUTPUTs:
            Ct: heat capacity of thawed soil (J m-3 C-1) 
            Cf: heat capacity of frozen soil (J m-3 C-1) 

    DEPENDENTs:
            None 
    """    
    
    import numpy as np

    s_data = np.genfromtxt(input_file, names = True, 
                      delimiter=',', dtype=None)
                    
    Bulk_Density_Texture = s_data['Bulk_Density'];
    Heat_Capacity_Texture = s_data['Heat_Capacity'];
    
    # Adjust percent of sand, silt, clay and peat ==
    
    tot_percent = p_sand+p_clay+p_silt+p_peat;
    
    percent_sand = p_sand / tot_percent;
    percent_clay = p_clay / tot_percent;
    percent_silt = p_silt / tot_percent;
    percent_peat = p_peat / tot_percent;
    
    # Calculate heat capacity and bulk density of soil using exponential weighted.
                             
    Heat_Capacity =      Heat_Capacity_Texture[2]*percent_clay + \
                         Heat_Capacity_Texture[1]*percent_sand + \
                         Heat_Capacity_Texture[0]*percent_silt + \
                         Heat_Capacity_Texture[3]*percent_peat       # Unit: J kg-1 C-1 
                       
    Bulk_Density  =      Bulk_Density_Texture[2]*percent_clay + \
                         Bulk_Density_Texture[1]*percent_sand + \
                         Bulk_Density_Texture[0]*percent_silt + \
                         Bulk_Density_Texture[3]*percent_peat        # Unit: kg m-3
    
    # Estimate heat capacity for composed soil
    # based on the empirical approaches suggested by Anisimov et al. (1997)
        
    Ct = Heat_Capacity*Bulk_Density + 4190.*vwc; # eq-15, Anisimov et al. 1997; Unit: J m-3 C-1
    Cf = Heat_Capacity*Bulk_Density + 2025.*vwc; # eq-15, Anisimov et al. 1997; Unit: J m-3 C-1
    
    return Ct, Cf

def Estimate_Soil_Thermal_Conductivity(input_file, p_clay, p_sand, p_silt, p_peat, vwc):

    """ 
    The function is to estimate soil thermal conductivity ;
    INPUTs:
            input_file: file contains typical thermal parameters of clay, sand, silt and peat.
            p_clay: percent of clay (%)
            p_sand: percent of sand (%)
            p_silt: percent of silt (%)
            p_peat: percent of peat (%)
            vwc   : volumetric water content (m3 / m3)

    OUTPUTs:
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 

    DEPENDENTs:
            None 
    """ 

    import numpy as np    
    
    s_data = np.genfromtxt(input_file, names = True, 
                      delimiter=',', dtype=None)
                      
    # Adjust percent of sand, silt, clay and peat ==
        
    tot_percent = p_sand+p_clay+p_silt+p_peat;
    
    percent_sand = p_sand / tot_percent;
    percent_clay = p_clay / tot_percent;
    percent_silt = p_silt / tot_percent;
    percent_peat = p_peat / tot_percent;
                    
    Thermal_Conductivity_Thawed_Texture = s_data['Thermal_Conductivity_Thawed']
    Thermal_Conductivity_Frozen_Texture = s_data['Thermal_Conductivity_Frozen']
    
    # Estimate thermal conductivity for composed soil 
    
    Kt_Soil =  Thermal_Conductivity_Thawed_Texture[0]**percent_silt * \
               Thermal_Conductivity_Thawed_Texture[2]**percent_clay * \
               Thermal_Conductivity_Thawed_Texture[1]**percent_sand * \
               Thermal_Conductivity_Thawed_Texture[3]**percent_peat
              
    Kf_Soil =  Thermal_Conductivity_Frozen_Texture[0]**percent_silt * \
               Thermal_Conductivity_Frozen_Texture[2]**percent_clay * \
               Thermal_Conductivity_Frozen_Texture[1]**percent_sand * \
               Thermal_Conductivity_Frozen_Texture[3]**percent_peat
    
    # Consider the effect of water content on thermal conductivity
                  
    Kt = Kt_Soil**(1.0-vwc)*0.54**vwc #   Unit: (W m-1 C-1)
    Kf = Kf_Soil**(1.0-vwc)*2.35**vwc #   Unit: (W m-1 C-1)
    
    return Kt, Kf
    
    
def Estimate_Cold_Warm_Season_Length(Ta, Aa):
    
    """ 
    The function is to calculate length of cold season and warm season;
    
    INPUTs:
            Ta: mean annual air temperature (C);
            Aa: amplitude of annual air temperature (C);

    OUTPUTs:
            tao : period of a year (second)
            tao1: period of cold season; (second)
            tao2: period of warm season; (second)
    DEPENDENTs:
            None 
    """
    
    import numpy as np
    
    tao = 365.*24.*3600.; # seconds in a year
    
    tao1 = tao*(0.5 - 1./np.pi*np.arcsin(Ta/Aa)); # Cold Season , Page-129, Sazonova, 2003 
    tao2 = tao - tao1;                            # Warm Season , Page-129, Sazonova, 2003
    
    return tao, tao1, tao2

def Estimate_Latent_Heat(vwc):

    """ 
    The function is to calculate latent heat;
    
    INPUTs:
            vwc: volumetric water content (m3 / m3)

    OUTPUTs:
            L : volumetric latent heat of the water  (J/kg)

    DEPENDENTs:
            None 
    """
        
    L = 3.34E8*vwc;                  # Latent Heat, Unit: J kg-1 C-1;
                                     # eq-16, Anisimov et al. 1997
    
    return L
    
def Estimate_Snow_Effect(Ta, Aa, Hsn, Csn, Ksn, rho_sn, tao):

    """ 
    The function is to calculate snow effect;
    
    INPUTs:
            Ta: mean annual air temperature (C)
            Aa: amplitude of air temperature (C)
            Hsn: Snow depth (m)
            Csn: heat capacity of snow  (J m-3 C-1)
            Ksn: thermal conductivity of snow  (W m-1 C-1)
            rho_sn: snow density (kg m-3).             
            tao : period of a year (second)
    OUTPUTs:
            Tvg : Temperature on the top of vegetation (C)
            Avg : Amplitude of temperature on the top of vegetation (C)

    DEPENDENTs:
            None 
    """
    
    import numpy as np
    
    #   Estimating Snow Effects
    
    K_diffusivity = Ksn/(rho_sn*Csn)   
    
#    print K_diffusivity
    
#    print (1.0 - np.exp(-1.0*Hsn*np.sqrt(np.pi/(tao*K_diffusivity))))
    
    deta_Tsn = Aa*(1.0 - np.exp(-1.0*Hsn*np.sqrt(np.pi/(tao*K_diffusivity)))); # eq-7, Anisimov et al. 1997
    deta_Asn = 2.0/np.pi*deta_Tsn; # eq-2, Sazonova et al., 2003
    
    # mean annual temperature and amplitude 
    # bellow snow OR top of vegetation
    Tvg = Ta + deta_Tsn;    # Page-129, Sazonova et al., 2003
    Avg = Aa - deta_Asn;    # Page-129, Sazonova et al., 2003
    
    return Tvg, Avg

def Estimate_Vegetation_Effect(Tvg, Avg, Hvgf, Hvgt, Dvf, Dvt, tao1, tao2, tao):
    
    """ 
    The function is to calculate vegetation effect;
    
    INPUTs:
            Tvg : Temperature on the top of vegetation (C)
            Avg : Amplitude of temperature on the top of vegetation (C)
            Hvgf: Height of vegetation in frozen period (m)
            Hvgt: Height of vegetation in thawed period (m)
            Dvf : Thermal diffusivity of vegetation in frozen period (m2 s)
            Dvt : Thermal diffusivity of vegetation in thawed period (m2 s)
            tao : period of a year (second)
            tao1: period of cold season; (second)
            tao2: period of warm season; (second)
    OUTPUTs:
            Tgs : Temperature on the ground surface (C)
            Ags : Amplitude of temperature on the ground surface (C)

    DEPENDENTs:
            None 
    """
    
    import numpy as np

    # Estimate vegetation effects
    
    # WINTER vegetation thermal effects:  
    #   eq-10, Anisimov et al. 1997
    deta_A1 = (Avg - Tvg) * \
              (1.-np.exp(-1.*Hvgf*np.sqrt(np.pi/(Dvf*2.*tao1))));
              
    # SUMMER vegetation thermal effects:  
    #   eq-11, Anisimov et al. 1997        
    deta_A2 = (Avg  + Tvg) * \
              (1.-np.exp(-1.*Hvgt*np.sqrt(np.pi/(Dvt*2.*tao2))));
              
    # Effects of vegetation on seasonal amplitude of temperature: 
    #   eq-8, Anisimov et al. 1997
    deta_Av = (deta_A1*tao1+deta_A2*tao2) / tao;
    
    # Effects of vegetation on annual mean temperature: 
    #   eq-9, Anisimov et al. 1997    
    deta_Tv = (deta_A1*tao1-deta_A2*tao2) / tao * (2. / np.pi);
    
    # mean annual temperature and amplitude 
    # on the ground surface
    Tgs = Tvg + deta_Tv;    # eq-13, Sazonova et al., 2003
    Ags = Avg - deta_Av;    # eq-14, Sazonova et al., 2003
    
    return Tgs, Ags
    
def Calculate_Tps_numerator(Tgs, Ags, Kf, Kt):
    
    """ 
    The function is to calculate Tps_Numerator;
    
    INPUTs:
            Tgs : Temperature on the ground surface (C)
            Ags : Amplitude of temperature on the ground surface (C)
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 
            
    OUTPUTs:
            Tps_numerator : Numerator in Equation of
                            Temperature on the top of permafrost

    DEPENDENTs:
            None 
    """    
    
    import numpy as np
        
    Tps_numerator = 0.5*Tgs*(Kf+Kt)\
                    +(Ags*(Kt-Kf)/np.pi\
                    *(Tgs/Ags*np.arcsin(Tgs/Ags)\
                    +np.sqrt(1.-(np.pi**2.0/Ags**2.0)))); # eq-14, Anisimov et al. 1997
                    
    return Tps_numerator

def Calculate_Temperature_Top_PF(Tps_numerator, Kf, Kt, Cf, Ct):
    
    """ 
    The function is to calculate temperature on the top of permafrost;
    
    INPUTs:
            Tps_numerator : Numerator in Equation of
                            Temperature on the top of permafrost
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 
            Ct: heat capacity of thawed soil (J m-3 C-1) 
            Cf: heat capacity of frozen soil (J m-3 C-1) 
            
    OUTPUTs:
            Tps : Temperature on the top of permafrost (C)

    DEPENDENTs:
            None 
    """ 
   
#   Calculating the temperature at the top of permafrost :
                    
    if Tps_numerator<=0.0: # It is for PERMAFROST
        K_star = Kf;
        #print "PERMAFROST"
    else:                  # It is for SEASONAL FROZEN GROUND 
        K_star = Kt;
        #print "SEASONAL FROZEN GROUND"
               
    # Temperature at the top of permafrost    
       
    Tps = Tps_numerator/K_star;  # eq-14, Anisimov et al. 1997
    
    return Tps
    
def Calculate_ALT(Tps_numerator, Tps, Ags, Kf, Kt, Cf, Ct, L, tao):
    
    """ 
    The function is to calculate active layer thickness;
    
    INPUTs:
            Tps_numerator : Numerator in Equation of
                            Temperature on the top of permafrost
            Kt: thermal conductivity of thawed soil (W m-1 C-1) 
            Kf: thermal conductivity of frozen soil (W m-1 C-1) 
            Ct: heat capacity of thawed soil (J m-3 C-1) 
            Cf: heat capacity of frozen soil (J m-3 C-1) 
            L : volumetric latent heat of the water  (J/kg)
            tao : period of a year (second)
            
    OUTPUTs:
            Zal : Active Layer thickness (m)

    DEPENDENTs:
            None 
    """ 
    
    import numpy as np
    
    if Tps_numerator<=0.0: # It is for PERMAFROST
        C = Cf;
        K = Kf;
    else:                   # It is for SEASONAL FROZEN GROUND 
        C = Ct;
        K = Kt;
    
    Aps = (Ags - abs(Tps))/np.log((Ags+L/(2.*C)) / \
          (abs(Tps)+L/(2.*C))) - \
          L/(2.*C);                                                           # eq-4, Romanovsky et al. 1997
    
    Zc = (2.*(Ags - abs(Tps))*np.sqrt((K*tao*C)/np.pi)) / (2.*Aps*C + L);     # eq-5, Romanovsky et al. 1997
    
    Zal = (2.*(Ags - abs(Tps))*np.sqrt(K*tao*C/np.pi)\
        +(((2.*Aps*C*Zc+L*Zc)*L*np.sqrt(K*tao/(np.pi*C)))\
        /(2.*Ags*C*Zc + L*Zc +(2.*Aps*C+L)*np.sqrt(K*tao/(np.pi*C)))))\
        /(2.*Aps*C+L);                                                    # Active Layer Thickness, eq-3, Romanovsky et al. 1997
        
    if Tps_numerator>0: # It is seasonal frost depth, not be considered here.
        
        Zal = -999.99
        Tps = -999.99
    
    return Zal
    
def Extract_Soil_Texture(input_lat, input_lon, 
                         lon_grid, lat_grid, 
                         Clay_percent, Sand_percent, Silt_percent, Peat_percent):
    
    """ 
    The function is to extract the soil texture according to input of latitude and longitude;
    INPUTs:
            input_lat: Latitude;
            input_lon: Longitude;
    OUTPUTs:
            p_clay: percent of clay (%)
            p_sand: percent of sand (%)
            p_silt: percent of silt (%)
            p_peat: percent of peat (%)             
    DEPENDENTs:
            function "Extract_Grid_Value"
    """

#    Clay_file = 'Parameters/T_CLAY.nc4';
#    Sand_file = 'Parameters/T_SAND.nc4';
#    Silt_file = 'Parameters/T_SILT.nc4';
#    Peat_file = 'Parameters/T_OC.nc4';
    
    clay_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Clay_percent)
                         
    sand_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Sand_percent)
    
    silt_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Silt_percent)
                                   
    peat_perc = Extract_Grid_Value(input_lat, input_lon, 
                                   lon_grid, lat_grid, Peat_percent)
                         
    return clay_perc, sand_perc, silt_perc, peat_perc
    
def import_ncfile(input_file, lonname,  latname,  varname):
    
    from netCDF4 import Dataset
    #    import numpy as np
        
        # Read the nc file 
        
    fh = Dataset(input_file, mode='r')
        
        # Get the lat and lon
        #   Set the grid size for lat. and lon. (here is 0.5 degree)
        
    lon_grid = fh.variables[lonname][:]; 
    lat_grid = fh.variables[latname][:];
        
    p_data  = fh.variables[varname][:];
        
    return lat_grid,lon_grid,p_data
    
def Extract_Grid_Value(input_lat, input_lon, lon_grid, lat_grid, p_data): 

    """ 
    The function is to extract the grid value from NetCDF file,
    according to input of latitude and longitude;
    
    INPUTs:
            input_lat: Latitude;
            input_lon: Longitude;
            input_file: grid data file (NetCDF file)
            lonname: name of longitude in "input_file"
            latname: name of latitude in "input_file"
            lon_grid_scale: grid size of longitude
            lat_grid_scale: grid size of latitude
            varname: name of variable should be extracted.
            
    OUTPUTs:
            p_data: grid value   
                    
    DEPENDENTs:
            None 
    """
    
#    from netCDF4 import Dataset
    import numpy as np
    
    lon_grid_scale = 0.05;
    lat_grid_scale = 0.05;
    
    lon_grid_top = lon_grid + lon_grid_scale / 2.0;
    lat_grid_top = lat_grid + lat_grid_scale / 2.0;
    
    lon_grid_bot = lon_grid - lon_grid_scale / 2.0;
    lat_grid_bot = lat_grid - lat_grid_scale / 2.0;
    
    # Get the index of input location acccording to lat and lon inputed
    
    idx_lon = np.where((input_lon <= lon_grid_top) & (input_lon >= lon_grid_bot))          
    idx_lat = np.where((input_lat <= lat_grid_top) & (input_lat >= lat_grid_bot))
    
    idx_lon = np.array(idx_lon)
    idx_lat = np.array(idx_lat)
    
    if np.size(idx_lon) >= 1 and np.size(idx_lat) >= 1:
        q_data  = p_data[idx_lat[0,0], idx_lon[0,0]]
    else:
        q_data  = np.nan;
    
    return q_data
    
    