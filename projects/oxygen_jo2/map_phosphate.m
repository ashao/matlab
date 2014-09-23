phos_in.file = '/home/ashao/uw-apl/data/offtrac/input/normalyear/woalevpo4.nc'
phos_out.file = '/home/ashao/uw-apl/data/offtrac/input/normalyear/p_an.himgrid.monthly.nc';

phos_in.lat = nc_varget(phos_in.file,'lat');
phos_in.lon = nc_varget(phos_in.file,'lon');
phos_in.conc = nc_varget(phos_in.file,'p_an');

load metrics

[ntime nlay nlat nlon] = size(phos_in.conc);

%%

phos_out.conc = zeros([ntime,nlay,size(metrics.geolon.data)]);
templon = mod(phos_in.lon+180,360)-180;
interplon = [phos_in.lon-360 ; phos_in.lon];
[longrid latgrid] = meshgrid(interplon,phos_in.lat);
%%

for mon = 1:ntime
    fprintf('Month %d Layer: ',mon);
    for lay = 1:nlay
        fprintf('%d...',lay);
        data = squeeze(phos_in.conc(mon,lay,:,:));
        data = [data data];
        data(isnan(data))=0.0;
        phos_out.conc(mon,lay,:,:) = interp2( ...
            double(interplon),double(phos_in.lat),double(data), ...
            metrics.geolon.data,metrics.geolat.data);        
    end
    fprintf('DONE!\n');
end
%%
phos_out.conc(:,:,~logical(metrics.wet.data))=0;
phos_out.conc(isnan(phos_out.conc))=0;
nc_varput(phos_out.file,'p_an',phos_out.conc)