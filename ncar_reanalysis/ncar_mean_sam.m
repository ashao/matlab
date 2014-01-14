function array = ncar_mean_sam(ufile,vfile,samfile)

if nargin<2
    ufile='/scratch/data/ncep_reanalysis/daily/combined/uwnd.sig995.all.nc';
    vfile='/scratch/data/ncep_reanalysis/daily/combined/vwnd.sig995.all.nc';    
    samfile='/scratch/data/climate_indices/aao.mat';
end

array.lat=nc_varget(ufile,'lat');
array.lon=nc_varget(ufile,'lon');
array.time=datenum(1,1,1)+nc_varget(ufile,'time')/24;
array.uwnd=nc_varget(ufile,'uwnd');
array.vwnd=nc_varget(vfile,'vwnd');
[array.longrid array.latgrid]=meshgrid(array.lon,array.lat);
load(samfile)

ncar_sam_idx=interp1(aao.dates,aao.index,array.time);
posidx=find(ncar_sam_idx>=2);
negidx=find(ncar_sam_idx<=-2);

% length(posidx)
% length(negidx)
% pause

array.sam_pos_mean_u=squeeze(nanmean(array.uwnd(posidx,:,:),1));
array.sam_pos_mean_v=squeeze(nanmean(array.vwnd(posidx,:,:),1));
array.sam_neg_mean_u=squeeze(nanmean(array.uwnd(negidx,:,:),1));
array.sam_neg_mean_v=squeeze(nanmean(array.vwnd(negidx,:,:),1));

array.sam_pos_std_u=squeeze(nanstd(array.uwnd(posidx,:,:),1));
array.sam_pos_std_v=squeeze(nanstd(array.vwnd(posidx,:,:),1));
array.sam_neg_std_u=squeeze(nanstd(array.uwnd(negidx,:,:),1));
array.sam_neg_std_v=squeeze(nanstd(array.vwnd(negidx,:,:),1));

end