clc; clear; close all;

Forcing_dataset_filename = 'Input/ERA_Reanalysis_NH.nc';
Land_Mask_filename = 'Input/Land_Mask_NH.nc';

% Make Mask ===
Land_Sea_Mask = ncread(Land_Mask_filename, 'lsm');
Land_Sea_Mask(Land_Sea_Mask==0) = NaN;

% Get Latitudes and Longitudes ===
Lat_list = ncread(Forcing_dataset_filename, 'latitude');
Lon_list = ncread(Forcing_dataset_filename, 'longitude');

Lon_list = mod((Lon_list + 180),360) - 180;

Lat = repmat(Lat_list,1,720); Lat = Lat';
Lon = repmat(Lon_list,1,116);

year = 2015;

results = dlmread(['Output/',num2str(year,'%d'),'_results.csv'],',',0,0);
results(results==-999.99) = NaN;
results(results<=0) = NaN;

nanmean(nanmean(results(:,3)))
nanmax(nanmax(results(:,3)))
nanmin(nanmin(results(:,3)))

results = reshape(results(:,3), 720, 116);

h = figure;

clev = [0:0.1:1.0 1.25:0.25:2.0 2.0:0.5:4];

m_proj('stereographic','lat',90,'long',180,'radius',60);
m_contourf(Lon_list, Lat_list, results',clev,'linestyle','none');
m_coast();
% m_plotbndry();
m_grid('box','fancy','tickdir','in');
colormap((jet));
cb = colorbar;
cb.Label.String = 'Active Layer Thickness (m)';

set_axis

%% Print Out the figure:
position = [0 0 15 15];
papersize= [15 15];
set(h,'paperunits','centimeters');
set(h,'PaperSize',papersize)
set(h,'paperposition',position);
print(h, '-dpng', ['ALT_',num2str(year,'%d'),'.png']);