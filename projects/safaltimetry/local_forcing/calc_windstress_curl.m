%% Load ECMWF ERA-Interim
ecmwfpath = '/ltraid4/ecmwf/era-interim/monthly/taux.tauy.1979.2013.nc';
variablenames = {'latitude','longitude','time','inss','iews'};
for var = 1:length(variablenames)
    varname = variablenames{var};
    fprintf('Extracting %s\n',varname)
    ecmwf.(varname) = double(nc_varget(ecmwfpath,varname));    
end
%%
maskpath = '/ltraid4/ecmwf/era-interim/landmask.nc';
landmask = nc_varget(maskpath,'lsm');
[ntime nlat nlon] = size(ecmwf.inss);

%% Calculate windstress curl and store
ecmwf.taucurl = zeros(size(ecmwf.inss));
[longrid latgrid] = meshgrid(ecmwf.longitude,ecmwf.latitude);
for t = 1:ntime
    print_progress(t,ntime,20)
        [null null null tauydx] = gradientm(latgrid, ...
        longrid, squeeze(ecmwf.inss(t,:,:)));
    [null null tauxdy null] = gradientm(latgrid, ...
        longrid, squeeze(ecmwf.iews(t,:,:)));
    ecmwf.taucurl(t,:,:) = tauydx - tauxdy;
    
end
%%
ecmwf.taucurl(:,landmask==1)=NaN;
[topo.lat topo.lon topo.z] = satbath(10,[-90 -30],[0 359.9]);
topo.z(topo.z>0) = 0;
%%
clf
worldmap([-90 -30],[0 360])

pcolorm(topo.lat,topo.lon,topo.z);
% contourm(ecmwf.latitude,ecmwf.longitude,squeeze(ecmwf.taucurl(1,:,:))*1e7,[-1 1],'LineColor','Black')
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])