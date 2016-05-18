clc; clear; close all;

Forcing_dataset_filename = 'ERA_Reanalysis_NH.nc';
Land_Mask_filename = 'Land_Mask_NH.nc';

Tair_filename = 'ERA_Reanalysis_NH.nc';
Snod_filename = 'ERA_Reanalysis_NH.nc';
Vwc_filename  = 'ERA_Reanalysis_NH.nc';

% Make Mask ===
Land_Sea_Mask = ncread(Land_Mask_filename, 'lsm');
Land_Sea_Mask(Land_Sea_Mask==0) = NaN;

% Get the time variable ===

Time = ncread(Forcing_dataset_filename, 'time');
Time_Units = 'hours since 1900-01-01 00:00:0.0';
Time = Time/24 + datenum(1900,1,1);
Time = datevec(double(Time));
Time = Time(:,1:2);

Year_all = unique(Time(:,1));

% Get Latitudes and Longitudes ===
Lat_list = ncread(Forcing_dataset_filename, 'latitude');
Lon_list = ncread(Forcing_dataset_filename, 'longitude');

Lon_list = mod((Lon_list + 180),360) - 180;

Lat = repmat(Lat_list,1,720); Lat = Lat';
Lon = repmat(Lon_list,1,116);

% For Each Year ===
for iiii = 1:37
    
    year = Year_all(iiii); % Taking 1979 as an example
    
    idx  = find(Time(:,1)==year);
    
    % Get Monthly Air Temperature ===
    
    Tair = ncread(Tair_filename, 't2m', [1,1,idx(1)], [720, 116, numel(idx)]);
    Tair = Tair - 273.15;
    Tair = Tair .* repmat(Land_Sea_Mask,[1,1, numel(idx)]);
    Tair_mean = nanmean(Tair, 3);
    Tair_amplitude = (nanmax(Tair,[],3) - nanmin(Tair,[],3))/2.0;
    
    % Get Monthly Snow Depth ===
    
    Snod = ncread(Snod_filename, 'sd', [1,1,idx(1)], [720, 116, numel(idx)]);
    Snod = Snod .* repmat(Land_Sea_Mask,[1,1, numel(idx)]);
    Rsno = ncread(Snod_filename, 'rsn', [1,1,idx(1)], [720, 116, numel(idx)]);
    Rsno(Rsno<100) = NaN;
    Snod = Snod ./ (Rsno * 0.001); Snod(Snod==0) = NaN; Rsno(Snod==0) = NaN;
    Snod = nanmean(Snod, 3); % Average snow depth
    Rsno = nanmean(Rsno,3);  % Average snow density
    
    % Get VWC ===
    
    Vwc = ncread(Snod_filename, 'swvl1', [1,1,idx(1)], [720, 116, numel(idx)]);
    Vwc = Vwc .* repmat(Land_Sea_Mask,[1,1, numel(idx)]);
    Vwc = nanmax(Vwc,[], 3); % Maximum VWC
    
    % Write out===
    n_grid  = numel(Land_Sea_Mask);
    lat_out = reshape(Lat,n_grid,1);
    lon_out = reshape(Lon,n_grid,1);
    Ta_out  = reshape(Tair_mean,n_grid,1);
    Aa_out  = reshape(Tair_amplitude,n_grid,1);
    Sd_out  = reshape(Snod,n_grid,1);
    rho_out = reshape(Rsno,n_grid,1);
    VWC_out = reshape(Vwc,n_grid,1);
    Hvgf_out = zeros(n_grid,1);
    Hvgt_out = zeros(n_grid,1);
    Dvf_out  = ones(n_grid,1)*1.39E-6;
    Dvt_out  = ones(n_grid,1)*5.56E-8;
    
    fid = fopen([num2str(year,'%d'),'.csv'],'wt');
    
    fprintf(fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',...
        'Lat','Lon', 'Ta', 'Aa', 'Sd', 'rho', 'VWC',...
        'Hvgf','Hvgt','Dvf','Dvt');
    
    fprintf(fid, '%f,%f,%f,%f,%f,%d,%f,%f,%f,%0.12f,%0.12f\n', ...
        [lat_out lon_out Ta_out Aa_out Sd_out rho_out VWC_out ...
        Hvgf_out Hvgt_out Dvf_out Dvt_out]');
    
    fclose(fid);
    
end