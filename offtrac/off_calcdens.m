function [ sigmat ] = off_calcdens(datapath)
if nargin<1
    datapath='/scratch/data/offtrac/input/';
end

%% Calculate sigma-t for mixed layer
start4d=zeros(1,4);
end4d=[12 -1 -1 -1];

T=nc_varget([datapath 'ts.nc'],'temp',start4d,end4d);
S=nc_varget([datapath 'ts.nc'],'salt',start4d,end4d);
H=nc_varget([datapath 'H-clim.nc'],'h',[1 0 0 0],end4d);
lath=nc_varget([datapath 'metrics.nc'],'lath');
lonh=nc_varget([datapath 'metrics.nc'],'lonh');
geolat=nc_varget([datapath 'metrics.nc'],'geolat');
geolon=nc_varget([datapath 'metrics.nc'],'geolon');
wet=nc_varget([datapath 'metrics.nc'],'wet');
depth=cumsum(H,2);


[ntime nlay nlat nlon]=size(T);
geolat4d(1,1,:,:)=geolat;
geolat4d=repmat(geolat4d,[ntime nlay 1 1]);

pres=sw_pres( depth, geolat4d );

sigmat=sw_pden(S,T,pres,zeros(size(pres)));


end