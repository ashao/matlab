function array = off_comp_inventory(ncfile,metricfile)

array.time=nc_varget(ncfile,'Time');
Ah(1,:,:)=nc_varget(metricfile,'Ah');
wet=nc_varget(metricfile,'wet');
Ah(1,~wet)=NaN;
ntime=length(array.time)
Ah_3dgrid=repmat(Ah,[49 1 1]);
Ah_2dgrid=squeeze(Ah);

% h=waitbar(0,'Extracting data from timestep');

array.cfc11_mon_global_flux=zeros(ntime,1);
array.cfc12_mon_global_flux=zeros(ntime,1);
array.sf6_mon_global_flux=zeros(ntime,1);
array.cfc11_global_inv=zeros(ntime,1);
array.cfc12_global_inv=zeros(ntime,1);
array.sf6_global_inv=zeros(ntime,1);

cfc11_mon_global_flux=zeros(ntime,1);
cfc12_mon_global_flux=zeros(ntime,1);
sf6_mon_global_flux=zeros(ntime,1);
cfc11_global_inv=zeros(ntime,1);
cfc12_global_inv=zeros(ntime,1);
sf6_global_inv=zeros(ntime,1);

seconds_in_day=60*60*24;

for t=1:ntime
    tic;
    disp(sprintf('%d/%d',t,ntime))
    mon=mod(t,12);
    if mon==0
        mon=12;
    end
    
    start4d=[t-1 0 0 0];
    count4d=[1 -1 -1 -1];
    start3d=[t-1 0 0];
    count3d=[1 -1 -1];
    
%     waitbar(t/ntime,h,sprintf('Month %d/%d',t,ntime));
    hsnap=nc_varget(ncfile,'mn_h',start4d,count4d);
    mldepth=squeeze(sum(hsnap(1:2,:,:),1));
    vol_grid=hsnap.*Ah_3dgrid;
    
    %% Calculate global flux
    c11flux=nc_varget(ncfile,'mn_c11flux',start3d,count3d).*Ah_2dgrid.*mldepth;
    c12flux=nc_varget(ncfile,'mn_c12flux',start3d,count3d).*Ah_2dgrid.*mldepth;
    sf6flux=nc_varget(ncfile,'mn_sf6flux',start3d,count3d).*Ah_2dgrid.*mldepth;
    cfc11_mon_global_flux(t)=nansum(nansum(c11flux));
    cfc12_mon_global_flux(t)=nansum(nansum(c12flux));
    sf6_mon_global_flux(t)=nansum(nansum(sf6flux));
    
    %% Calculate current global inventory
    
    mn_cfc11=nc_varget(ncfile,'mn_cfc11',start4d,count4d).*vol_grid;
    mn_cfc12=nc_varget(ncfile,'mn_cfc12',start4d,count4d).*vol_grid;
    mn_sf6=nc_varget(ncfile,'mn_sf6',start4d,count4d).*vol_grid;
    
    cfc11_global_inv(t)=nansum(nansum(nansum(mn_cfc11)));
    cfc12_global_inv(t)=nansum(nansum(nansum(mn_cfc12)));
    sf6_global_inv(t)=nansum(nansum(nansum(mn_sf6)));        
    toc;
    
end
array.cfc11_mon_global_flux=(cfc11_mon_global_flux);
array.cfc12_mon_global_flux=(cfc12_mon_global_flux);
array.sf6_mon_global_flux=(sf6_mon_global_flux);
array.cfc11_global_inv=cfc11_global_inv;
array.cfc12_global_inv=cfc12_global_inv;
array.sf6_global_inv=sf6_global_inv;

end