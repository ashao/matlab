infile = '/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/o_an.monthly.nc';

oxy_in.lat = double(nc_varget(infile,'lat'));
oxy_in.lon = double(nc_varget(infile,'lon'));
oxy_in.data = double(nc_varget(infile,'o_an'));
load metrics
wetmask = logical(metrics.wet.data);

%%
clf
axesm('MapProjection','Robinson','MapLatLimit',[-90 90],'MapLonLimit',[-280 80]);
contourfm(oxy_in.lat,oxy_in.lon,squeeze(oxy_in.data(1,1,:,:)))

%% Fill in NaNs due to the bathymetry being too shallow

[nmonths nlay nlat nlon] = size(oxy_in.data);

for mon = 1:nmonths
    for lay = 2:nlay
        
        nanidx = isnan(squeeze(oxy_in.data(mon,lay,:,:)));
        oxy_in.data(mon,lay,nanidx)=oxy_in.data(mon,lay-1,nanidx);
        
    end
end

oxy_in.data = mean(oxy_in.data);
[nmonths nlay nlat nlon] = size(oxy_in.data);
%% Perform bilinear interpolation for most of the ocean
[longrid latgrid] = meshgrid(oxy_in.lon,oxy_in.lat);
longrid = [longrid-360 longrid ];
latgrid = [latgrid latgrid];

oxy_out.data = zeros(nmonths,nlay,210,360);
geolonpt = metrics.geolon.data(wetmask);
geolatpt = metrics.geolat.data(wetmask);

for mon = 1:nmonths
    for lay = 1:nlay
        
        indata = squeeze(oxy_in.data(mon,lay,:,:));
        indata = [indata indata];
        oxy_out.data(mon,lay,wetmask) = interp2(longrid,latgrid,indata,geolonpt,geolatpt);
        
    end
end

%% Interpolate by the weighted average of a 2-degree box

for imon = 1:nmonths
    for ilay = 1:nlay
        data = squeeze(oxy_out.data(imon,ilay,:,:));
        intpdata = squeeze(oxy_in.data(imon,ilay,:,:));
        intpdata = [intpdata intpdata];
        interpidx = find( (isnan(data) | data==0) & wetmask);
        npts = length(interpidx);
        fprintf('Month: %d Layer %d NPoints: %d ... ',imon,ilay,npts)
        tic;
        geolonpt = metrics.geolon.data(interpidx);
        geolatpt = metrics.geolat.data(interpidx);
        temparray = zeros(npts,1);
        
        for ipt = 1:npts
            interp_scale = 1;
            indata = [];
            
            while isempty(indata) | temparray(ipt)==0
                
                dist = sqrt( (longrid-geolonpt(ipt)).^2 + ...
                    (latgrid - geolatpt(ipt)).^2);
                idx = find(dist <= interp_scale);
                
                indata = intpdata(idx);
                notnan = ~isnan(indata);
                indata = indata(notnan);
                interp_scale = interp_scale+1;                
            
            
                dist = dist(idx);
                wts = exp( -dist.^2 / (2*interp_scale) );
                wts = wts(notnan);
                wts = wts/sum(wts);
                
                temparray(ipt) = sum(wts.*indata);
            end
                        
        end        
        oxy_out.data(imon,ilay,interpidx)=temparray;
        fprintf('Finished in %e seconds\n',toc);
    end
end
%%

save /ltraid3/ashao/uw-apl/projects/oxygen_jo2/oxy_init.mat oxy_out
