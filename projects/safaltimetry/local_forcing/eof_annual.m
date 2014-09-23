function [ pc ] = eof_annual( data, time, lat, lon, latrange, lonrange)

ntime = length(time);
latidx = lat >= min(latrange) & lat < max(latrange);
lonidx = lon >= min(lonrange) & lon < max(lonrange);

lat = lat(latidx);
lon = lon(lonidx);

data = data(:,latidx,lonidx);
[longrid latgrid] = meshgrid(lon,lat);
wts = cosd(latgrid);
wts = wts / sum(wts(:));

wts_grid(1,:,:) = wts;
wts_grid = repmat(wts_grid,[ntime 1 1]);

% whos
data = data.*wts_grid;
data = data(:,:);

data = detrend(data')';

[U S V] = svds(data,1);
pcs = U*S;

sidx = 1;
eidx = 12;
counter = 0;
while eidx <= ntime
    
    counter = counter +1;
    pc(counter) = mean(pcs(sidx:eidx));
    sidx = sidx + 12;
    eidx = sidx + 11;
    
end

end