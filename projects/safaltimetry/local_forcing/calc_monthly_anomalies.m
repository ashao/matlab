taufile = '/ltraid4/ecmwf/era-interim/monthly/taux.tauy.1979.2013.nc';
fields = {'latitude','longitude','time','iews','inss'};
nfields = length(fields);
for t = 1:nfields
    field = fields{t};
    forcing.(field) = double(nc_varget(taufile,field));
end
forcing.time = forcing.time./24 + datenum(1900,1,1);
%% Calculate wind stress first
forcing.taucurl = zeros(size(forcing.inss));
[longrid latgrid] = meshgrid(forcing.longitude,forcing.latitude);
for t = 1:ntime
    print_progress(t,ntime,20)
        [null null null tauydx] = gradientm(latgrid, ...
        longrid, squeeze(forcing.inss(t,:,:)));
    [null null tauxdy null] = gradientm(latgrid, ...
        longrid, squeeze(forcing.iews(t,:,:)));
    forcing.taucurl(t,:,:) = tauydx - tauxdy;
    
end
%%
[ntime nlat nlon] = size(forcing.iews);
forcing.iews_anomaly = zeros(size(forcing.iews));
forcing.inss_anomaly = zeros(size(forcing.iews));
forcing.curl_anomaly = zeros(size(forcing.iews));
for mon=1:12
    nreps = length(mon:12:ntime);
    forcing.iews_anomaly(mon:12:end,:,:) = forcing.iews(mon:12:end,:,:) - ...
        repmat(mean(forcing.iews(mon:12:end,:,:)),[nreps 1 1]);
    forcing.inss_anomaly(mon:12:end,:,:) = forcing.inss(mon:12:end,:,:) - ...
        repmat(mean(forcing.inss(mon:12:end,:,:)),[nreps 1 1]);
    forcing.curl_anomaly(mon:12:end,:,:) = forcing.taucurl(mon:12:end,:,:) - ...
        repmat(mean(forcing.taucurl(mon:12:end,:,:)),[nreps 1 1]);
end
%% Truncate it to the right time period and then save
timeidx = forcing.time >= datenum(1993,1,1) & forcing.time < datenum(2013,1,1);
ecmwf.time = forcing.time(timeidx);
ecmwf.lat = forcing.latitude;
ecmwf.lon = forcing.longitude;
ecmwf.iews_anomaly = forcing.iews_anomaly(timeidx,:,:);
ecmwf.inss_anomaly = forcing.inss_anomaly(timeidx,:,:);
ecmwf.curl_anomaly = forcing.curl_anomaly(timeidx,:,:);
%%
save ~/uw-apl/projects/saf_altimetry/ecmwf.tauanom.mat ecmwf
