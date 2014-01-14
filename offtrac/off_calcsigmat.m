function sigmat = calc_sigmat(datapath)

if nargin<1
    datapath='/scratch/data/offtrac/input/';
end

%% Calculate sigma-t for mixed layer
start4d=zeros(1,4);
end4d=[12 4 -1 -1];

T=nc_varget([datapath 'ts.nc'],'TEMPCLIM',start4d,end4d);
S=nc_varget([datapath 'ts.nc'],'SALTCLIM',start4d,end4d);
H=nc_varget([datapath 'H-clim.nc'],'HCLIM',[1 0 0 0],end4d);
lath=nc_varget([datapath 'metrics.nc'],'lath');
lonh=nc_varget([datapath 'metrics.nc'],'lonh');
geolat=nc_varget([datapath 'metrics.nc'],'geolat');
geolon=nc_varget([datapath 'metrics.nc'],'geolon');
wet=nc_varget([datapath 'metrics.nc'],'wet');

[wetlat wetlon]=find(wet==0);
% T(:,:,wetlat,wetlon)=NaN;

depth=cumsum(H,2);
depth(:,1,:,:)=0;

[ntime nlay nlat nlon]=size(T);
ndat=nlay*nlat*nlon;

latvec(1,:,:)=geolat;
latvec=reshape(repmat(latvec,[nlay 1 1]),[ndat 1]);

sigmat=zeros(size(T))*NaN;

for i=1:ntime
    disp(sprintf('Month %d/%d',i,ntime))
    depthvec=reshape(depth(i,:,:,:),[ndat 1]);
    pres=reshape(sw_pres(depthvec,latvec),[nlay nlat nlon]);
    sigmat(i,:,:,:)=sw_pden(squeeze(S(i,:,:,:)),squeeze(T(i,:,:,:)),pres,2000*ones(size(pres)))-1000;
end
%% All layers below mixed layer have fixed densities, but are referenced to 2000dbar
sigma2=nc_varget([datapath 'Global_HIM_IC.nc'],'Layer')-1000;
sigma2mat(1,:,1,1)=sigma2(5:end);
sigma2mat=repmat(sigma2mat,[ntime 1 nlat nlon]);
nlay=49;
sigmat(:,5:nlay,:,:)=sigma2mat;

