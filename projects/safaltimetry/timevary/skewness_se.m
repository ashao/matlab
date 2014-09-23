stime = datenum(1993,1,1)-1;
offset = datenum(1950,1,1);
etime = datenum(2013,1,1);
inpath = '/ltraid4/aviso/alongtrack/sla/vxxc_matlab/';
files = dir([inpath 't*.mat']);
nfiles = length(files)

annual_skewness = zeros(254,20,3129);
latgrid = zeros(254,3129);

for fidx = 1:nfiles
    fprintf('Track %03d\n',fidx)
    load([inpath files(fidx).name])
    counter = 0;
    [ntime npt] = size(track.sla);
    avgtime = nanmean(track.time + datenum(1950,1,1),2);
    latgrid(fidx,1:npt) = track.lat(1:npt);
    for year = 1993:2012
        counter = counter + 1;       
        tidx = avgtime > datenum(year,1,1) & avgtime < datenum(year+1,1,1);        
        annual_skewness(fidx,counter,1:npt) = skewness(track.sla(tidx,1:npt));
    end
        
    
end
%%
clf
skewness_se = squeeze(nanstd(annual_skewness,0,2)./sqrt(1));
skewness_se(2:2:end,:)  = fliplr(skewness_se(1:2:end,:));
pcolor(skewness_se);shading flat
colorbar