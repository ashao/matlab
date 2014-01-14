%%
pfile = '/ltraid2/ashao/uw-apl/data/offtrac/input/normalyear/p_an.monthly.nc';
pdata = nc_varget(pfile,'p_an');
plat = nc_varget(pfile,'lat');
plon = nc_varget(pfile,'lon');
[nmonths nz nlon nlat] = size(pdata);

load metrics
outdata = zeros(nmonths,nz,210,360);
wet = logical(metrics.wet.data);

for mon = 1:12
    for zi = 1:nz
        fprintf('Month %d Layer%d\n',mon,zi);
        outdata(mon,zi,:,:) = interp2(plon,plat,squeeze(pdata(mon,zi,:,:)), ...
            mod(metrics.geolon.data,360),metrics.geolat.data);
    end
end

outdata(isnan(outdata))=0.0;
%%
outfile = '/ltraid2/ashao/uw-apl/data/offtrac/input/normalyear/p_an.himgrid.monthly.nc';
delete(outfile)
nc_create_empty(outfile);
nc_adddim(outfile,'time',12);
nc_adddim(outfile,'depth',nz);
nc_adddim(outfile,'lat',210);
nc_adddim(outfile,'lon',360);
time = cumsum(eomday(2000,1:12))-eomday(2000,1:12)/2;
v1.Name = 'p_an';
v1.Datatype = 'double';
v1.Dimension = {'time','depth','lat','lon'};
nc_addvar(outfile,v1);
nc_varput(outfile,'p_an',outdata)