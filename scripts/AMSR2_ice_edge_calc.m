clear all
close all

%% Read in the file 
addpath '/home/....'

% In a folder list all the AMSR2 files you want to loop over
MyFolderInfo = dir('*.nc');
 for j = 1:length(MyFolderInfo)
     grr = string(MyFolderInfo(j).name);  % List all files in the directory and loop over them
     file(:,j) = grr;

%% Read in the major variables

    lon = ncread(file(j),'longitude');
    lat = ncread(file(j),'latitude');
    ice = ncread(file(j),'sea_ice_concentration');
    land = ncread(file(j),'land');
    
%% Create a lat limit

    % this limits the region where the 0 % (or 15 %) are looked for
     latlim = [-58 -56];       % (just change this)
     
     inds = not(abs(sign(sign(latlim(1) - lat) + sign(latlim(2) - lat))));
     ind_1 = find(inds == 1);
     lon = lon(ind_1);
     lat = lat(ind_1);
     ice = ice(ind_1);
     land = land(ind_1);    
    
%% Read in the lonlimit
    
    % this limits the region where the 0 % (or 15 %) are looked for
    lonlim = [-5 5];           % (just change this)
    
    inds = not(abs(sign(sign(lonlim(1) - lon) + sign(lonlim(2) - lon))));
    ind_1 = find(inds == 1);
    lon = lon(ind_1);
    lat = lat(ind_1);
    ice = ice(ind_1);
    land = land(ind_1);
    
%% Create array - ice edge can be taken as 0 % SIC or 15 % SIC (just change)

    ind_land = find(land == 1);         % Find the points that are sea
    lon_new = lon(ind_land);            % Only have sea
    lat_new = lat(ind_land);            % Only have sea
    ice_new = ice(ind_land);            % Only have sea
    ind_ice = find(ice ==0);            % Find where the ice concentration is zero
    %ind_ice = find(ice >=15 & ice<=16); % Find where the ice concentration is 15 %
    lon_new = lon(ind_ice);             % 0 ice
    lat_new = lat(ind_ice);             % 0 ice
    ice_new = ice(ind_ice);             % 0 ice

    lon_unique = unique(lon_new);  % Find unique values of lon

    for i =1:length(lon_unique);
        ind_point = find(lon_unique(i) == lon_new);   % Find the indice location of chosen longitude
        lat_points = lat_new(ind_point);              % Index lat
        lat_min = min(lat_points);                    % find the value of the southern most latitude
        coord(j,:,i) = [lon_unique(i) lat_min];       % Cordinate of ice edge
    end
end
