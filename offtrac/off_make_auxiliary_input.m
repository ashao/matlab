infile.ocean = '/ltraid1/ashao/HIM/hyak_store/COMBINE/month/ocean_month.nc';
infile.ice = '/ltraid1/ashao/HIM/hyak_store/COMBINE/month/ice_month.nc';
infile.p = '/ltraid3/ashao/uw-apl/data/COREv2/clim/slp.15JUNE2009.nc';
infile.u10 = '/ltraid3/ashao/uw-apl/data/COREv2/clim/u_10.15JUNE2009.nc';
infile.v10 = '/ltraid3/ashao/uw-apl/data/COREv2/clim/v_10.15JUNE2009.nc';

fillarray4d = zeros(12,49,210,360);
fillarray3d = zeros(12,210,360);
temp = fillarray4d;
salt = fillarray4d;
ice = fillarray3d;
wind = fillarray3d;
p = fillarray3d;
%%
% Temperature, salt, ice, first since they don't need to be
% regrid
for mon = 1:12;
    fprintf('Month %d/%d: Temperature...',mon,12);
    start = [mon-1 0 0 0];
    count = inf*[1 1 1 1];
    stride = [12 1 1 1];
    temp(mon,:,:,:)=mean(nc_varget(infile.ocean,'temp', ... 
        start,count,stride));    
    fprintf('Salt...');
    salt(mon,:,:,:) = mean(nc_varget(infile.ocean,'salt', ... 
        start,count,stride));
    fprintf('Ice...')
    ice(mon,:,:) = mean(sum(nc_varget(infile.ice,'CN', ... 
        start,count,stride),2));     
    fprintf('Done!\n')
end
    
%% Regrid Wind and pressure
load metrics.mat
templon = nc_varget(infile.p,'LON');
templat = nc_varget(infile.p,'LAT');
temp_p = nc_varget(infile.p,'SLP');
temp_wind = sqrt(nc_varget(infile.u10,'U_10_MOD').^2 + ...
    nc_varget(infile.v10,'V_10_MOD').^2);
time = nc_varget(infile.p,'TIME');

% Fudge the northernmost latitude so that regridding works;
templat(templat==max(templat(:)))=max(metrics.geolat.data(:));
[templon templat] = meshgrid(templon,templat);

templon = [fliplr(templon-360) templon];
templat = [fliplr(templat) templat];

sday = [0 cumsum(eomday(2003,1:11))+1];
eday = cumsum(eomday(2003,1:12))+1;

for mon=1:12
    fprintf('Regridding wind: %d/%d\n',mon,12);
    idx = find(time>sday(mon) & time<eday(mon));
    wind_temp = squeeze(mean(temp_wind(idx,:,:)));
    wind_temp = [fliplr(wind_temp) wind_temp];
    wind(mon,:,:) = griddata(templon,templat,wind_temp, ...
        metrics.geolon.data, metrics.geolat.data);
    
    fprintf('Regridding pressure: %d/%d\n',mon,12);
    p_temp = squeeze(mean(temp_p(idx,:,:)));
    p_temp = [fliplr(p_temp) p_temp];
    p(mon,:,:) = griddata(templon,templat,p_temp, ...
        metrics.geolon.data, metrics.geolat.data);
end

%% Output all variables
outfile = '/ltraid3/ashao/uw-apl/data/offtrac/input/auxiliary_fields.nc';
delete(outfile)
nc_create_empty(outfile);

nc_adddim(outfile,'Time',0);
nc_adddim(outfile,'xh',360);
nc_adddim(outfile,'yh',210);
nc_adddim(outfile,'zl',49);

time = eomday(2002,1:12)/2+cumsum(eomday(2002,1:12))-31;
nc_addvar(outfile,nc_getvarinfo(infile.ocean,'Time'));
nc_varput(outfile,'Time',time);
nc_addvar(outfile,nc_getvarinfo(infile.ocean,'xh'));
nc_varput(outfile,'xh',nc_varget(infile.ocean,'xh'));
nc_addvar(outfile,nc_getvarinfo(infile.ocean,'yh'));
nc_varput(outfile,'yh',nc_varget(infile.ocean,'yh'));
nc_addvar(outfile,nc_getvarinfo(infile.ocean,'zl'));
nc_varput(outfile,'zl',nc_varget(infile.ocean,'zl'));
%
nc_addvar(outfile,nc_getvarinfo(infile.ocean,'temp'));
nc_varput(outfile,'temp',temp);
nc_addvar(outfile,nc_getvarinfo(infile.ocean,'salt'));
nc_varput(outfile,'salt',salt);
%
disp('SLP')
pinf = nc_getvarinfo(infile.p,'SLP'); % Store this as a prototype for the windspeed
pinf.Dimension={'Time' 'yh' 'xh'};
pinf.Size = [12 210 360];
pinf=rmfield(pinf,'Chunking');
pinf.Name = 'p_surf';
nc_addvar(outfile,pinf);
nc_varput(outfile,'p_surf',p);
%
disp('CN')
tempinf = nc_getvarinfo(infile.ice,'CN');
tempinf.Dimension={'Time' 'yh' 'xh'};
tempinf.Size = pinf.Size;
tempinf.Dimension = pinf.Dimension;
nc_addvar(outfile,tempinf);
nc_varput(outfile,'CN',ice);
%
disp('V_10')
tempinf = nc_getvarinfo(infile.v10,'V_10_MOD');
tempinf.Dimension = {'Time' 'yh' 'xh'};
tempinf.Size = pinf.Size;
tempinf=rmfield(tempinf,'Chunking');
tempinf.Name = 'wind';
nc_addvar(outfile,tempinf);
nc_varput(outfile,'wind',wind);


