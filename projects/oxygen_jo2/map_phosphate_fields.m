infile = '/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/p_an.monthly.nc';
phos_in.lat = double(nc_varget(infile,'lat'));
phos_in.lon = double(nc_varget(infile,'lon'));
phos_in.data = double(nc_varget(infile,'p_an'));

% %%
% clf
% axesm('MapProjection','Robinson','MapLatLimit',[-90 90],'MapLonLimit',[-280 80]);
% contourfm(phos_in.lat,phos_in.lon,squeeze(phos_in.data(1,1,:,:)))
% 
%% Fill in NaNs due to the bathymetry being too shallow

[nmonths nlay nlat nlon] = size(phos_in.data);

for mon = 1:nmonths
    for lay = 2:nlay
        
        nanidx = isnan(squeeze(phos_in.data(mon,lay,:,:)));
        phos_in.data(mon,lay,nanidx)=phos_in.data(mon,lay-1,nanidx);
        
    end
end

%% Interpolate by the weighted average of a 2-degree box
[longrid latgrid] = meshgrid(phos_in.lon,phos_in.lat);
longrid = [longrid-360 longrid ];
latgrid = [latgrid latgrid];

testdata = squeeze(phos_in.data(1,1,:,:));
testdata =[testdata testdata];
load metrics

wetidx = find(logical(metrics.wet.data));
ptgeolon = metrics.geolon.data(wetidx);
ptgeolat = metrics.geolat.data(wetidx);

npts = length(wetidx);
phos_out.data = zeros([nmonths,nlay,size(metrics.wet.data)]);


for imon = 1:nmonths
    for ilay = 1:nlay
        fprintf('Month: %d Layer: %d\n',imon,ilay);
        data = squeeze(phos_in.data(imon,ilay,:,:));
        data = [data data];
        
        temparray = zeros(npts,1);
        tic;
        for ipt = 1:npts
            
            doneflag = 1;
            interp_scale = 1;
            indata = [];
            
            while isempty(indata)
                
                dist = sqrt( (longrid-ptgeolon(ipt)).^2 + ...
                    (latgrid - ptgeolat(ipt)).^2);
                idx = find(dist <= interp_scale);
                
                indata = data(idx);
                notnan = ~isnan(indata);
                indata = indata(notnan);
                interp_scale = interp_scale+1;                
            end
            
                dist = dist(idx);
                wts = exp( -dist.^2 / (2*interp_scale) );
                wts = wts(notnan);
                wts = wts/sum(wts);
                
                temparray(ipt) = sum(wts.*indata);
            
                        
        end        
        phos_out.data(imon,ilay,wetidx)=temparray;
        fprintf('Finished in %e seconds\n',toc);
    end
end
save /ltraid3/ashao/uw-apl/projects/oxygen_jo2/ phos_out