function [ array ] = off_procsec(datapath, ncfile, months, latrange, lonrange, layers)
year=1936;
s_mon=min(months)-1;
n_mon=length(months);

% time=floor(s_mon/12)*

lath=nc_varget(ncfile,'lath');
lonh=nc_varget(ncfile,'lonh');

latidx=find(lath>=min(latrange) & lath<=max(latrange));
lonidx=find(lonh>=min(lonrange) & lonh<=max(lonrange));

s_latidx=min(latidx)-1;
n_lat=length(latidx);
s_lonidx=min(lonidx)-1;
n_lon=length(lonidx);


s_layer=min(layers)-1;
n_layer=length(min(layers):max(layers));


disp('Initializing Netcdf Arrays')

array.lath=lath(latidx);
array.lonh=lonh(lonidx);
array.months=months;
mldepth=nc_varget(ncfile,'mn_h',[s_mon 0 s_latidx s_lonidx],[n_mon 2 n_lat n_lon]);
array.time=nc_varget(ncfile,'Time');
array.geolat=nc_varget([datapath 'metrics.nc'],'geolat',[s_latidx s_lonidx],[n_lat n_lon]);
array.geolon=nc_varget([datapath 'metrics.nc'],'geolon',[s_latidx s_lonidx],[n_lat n_lon]);

disp('Initializing Fill Arrays')
fillmat_3D=single(zeros(length(months), n_lat, n_lon));
fillmat_4D=single(zeros(length(months),n_layer,n_lat,n_lon));
% size(fillmat_3D)

array.depth=fillmat_4D;

array.mn_cfc11=fillmat_4D;
array.mn_cfc12=fillmat_4D;
array.mn_sf6=fillmat_4D;
array.mn_h=fillmat_4D;
array.cfc11_sat_tracer=fillmat_4D;
array.cfc12_sat_tracer=fillmat_4D;

array.cfc11_relsat=fillmat_3D;
array.cfc12_relsat=fillmat_3D;
array.sf6_relsat=fillmat_3D;

array.cfc11_flux=fillmat_3D;
array.cfc12_flux=fillmat_3D;
array.sf6_flux=fillmat_3D;

array.mldepth=fillmat_3D;

time=year+floor(min(months)/12);

for i=1:length(months)
     disp(sprintf('%d / %d \n',i,length(months)))
    midx=mod(array.months(i),12);
    start_3D=[array.months(i) s_latidx s_lonidx];
    end_3D=[1 n_lat n_lon];
    start_4D=[array.months(i) s_layer s_latidx s_lonidx];
    end_4D=[1 n_layer n_lat n_lon];
    
    if midx==0
        midx=12;
    end
    
    array.mldepth(i,:,:)=single(squeeze(sum(mldepth(i,:,:,:),2)));
    time=time+eomday(1,midx)/365;

    array.depth(i,:,:,:)=cumsum(nc_varget(ncfile,'mn_h',start_4D,end_4D));
    
    %   Extract tracer concentrations and surface saturation value
    array.mn_cfc11(i,:,:,:)=single(nc_varget(ncfile,'mn_cfc11',start_4D,end_4D));
    array.mn_cfc12(i,:,:,:)=single(nc_varget(ncfile,'mn_cfc12',start_4D,end_4D));
    array.mn_sf6(i,:,:,:)=single(nc_varget(ncfile,'mn_sf6',start_4D,end_4D));
    array.mn_h(i,:,:,:)=single(nc_varget(ncfile,'mn_h',start_4D,end_4D));
    
    array.cfc11_flux(i,:,:)=single(nc_varget(ncfile,'mn_c11flux',start_3D,end_3D));
    array.cfc12_flux(i,:,:)=single(nc_varget(ncfile,'mn_c12flux',start_3D,end_3D));
    array.sf6_flux(i,:,:)=single(nc_varget(ncfile,'mn_sf6flux',start_3D,end_3D));
    
    cfc11sat=single(nc_varget(ncfile,'cfc11_sat',start_3D,end_3D));
    cfc12sat=single(nc_varget(ncfile,'cfc12_sat',start_3D,end_3D));
    sf6sat=single(nc_varget(ncfile,'sf6_sat',start_3D,end_3D));
    
    array.cfc11_relsat(i,:,:)=single(squeeze(array.mn_cfc11(i,1,:,:))./cfc11sat);
    array.cfc12_relsat(i,:,:)=single(squeeze(array.mn_cfc12(i,1,:,:))./cfc12sat);
    array.sf6_relsat(i,:,:)=single(squeeze(array.mn_sf6(i,1,:,:))./sf6sat);
    
    array.cfc11_sat_tracer(i,:,:,:)=single(nc_varget(ncfile,'cfc11relsat',start_4D,end_4D));
    array.cfc12_sat_tracer(i,:,:,:)=single(nc_varget(ncfile,'cfc12relsat',start_4D,end_4D));
    array.sf6_sat_tracer(i,:,:,:)=single(nc_varget(ncfile,'sf6relsat',start_4D,end_4D));
end

array.depth(:,1,:,:)=0;
