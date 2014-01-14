ncfile='/scratch/models/offtrac/runs/fukushima_inert/FUKUSHIMA_INERT.0722.nc';
datapath='/scratch/data/offtrac/input/';
geolat=nc_varget([datapath 'metrics.nc'],'geolat');
geolon=nc_varget([datapath 'metrics.nc'],'geolon');
time=nc_varget(ncfile,'Time');
ntime=length(time);

maxinv_loc=zeros(ntime,2);
maxconc_loc=zeros(ntime,3);
cs137_inv_time=zeros([ntime 210 360]);
for i=1:length(time)
    disp(sprintf('Month %d of %d',i,ntime))
    mn_cs137=nc_varget(ncfile,'mn_cs137',[i-1 0 0 0],[1 -1 -1 -1]);
    mn_h=nc_varget(ncfile,'mn_h',[i-1 0 0 0],[1 -1 -1 -1]);
    depth=cumsum(mn_h,2);
    [maxval maxidx]=max(reshape(mn_cs137(:),[prod(size(mn_cs137)) 1]));
    [depthidx latidx lonidx]=ind2sub(size(mn_cs137),maxidx);
    maxconc_loc(i,:)=[geolon(latidx,lonidx) geolat(latidx,lonidx) depth(depthidx,latidx,lonidx)];
    inv_cs137=squeeze(sum(mn_cs137.*mn_h,1));
    cs137_inv_time(i,:,:)=inv_cs137;
    [maxval maxidx]=max(inv_cs137(:));
    [maxi maxj]=ind2sub(size(geolat),maxidx);
    maxinv_loc(i,:)=[geolon(maxi,maxj) geolat(maxi,maxj)];    
end
    
