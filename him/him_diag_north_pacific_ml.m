function [ array,year_ssh ] = him_diag_temp(ncfile,metricfile,hfield,tfield,sfield)

lonidx=[135 80  135 80  108];
latidx=[158 158 128 128 143];
    array.time=nc_varget(ncfile,'Time'); 
    ntime=length(array.time);

%     Get metrics.nc information
    geolon=nc_varget(metricfile,'geolon');
    geolat=nc_varget(metricfile,'geolat');
    
    array.lon=geolon(latidx,lonidx);
    array.lat=geolat(latidx,lonidx);
    
%     Allocate array for global mean ssh
    array.temp=zeros(ntime,length(lonidx));
    array.salt=zeros(ntime,length(lonidx));
    array.mldepth=zeros(ntime,5);
    h=waitbar(0);
    
    for t=1:ntime
       waitbar(t/ntime,h,sprintf('Time Index %d/%d',t,ntime))
%        disp(sprintf('Month %d/%d',t,ntime))
       hnow=nc_varget(ncfile,hfield,[t-1 0 0 0],[1 2 -1 -1]);
       temp_grid=nc_varget(ncfile,tfield,[t-1 0 0 0],[1 2 -1 -1]);       
       salt_grid=nc_varget(ncfile,sfield,[t-1 0 0 0],[1 2 -1 -1]);             
       
       mldepth=squeeze(sum(hnow,1));
       for pt=1:length(latidx)
           pt_latidx=latidx(pt);
           pt_lonidx=lonidx(pt);
           mldepth_pt=mldepth(pt_latidx,pt_lonidx);
           temp_pt=temp_grid(:,pt_latidx,pt_lonidx).*hnow(:,pt_latidx,pt_lonidx);
           salt_pt=salt_grid(:,pt_latidx,pt_lonidx).*hnow(:,pt_latidx,pt_lonidx);
           array.temp(t,pt)=sum(temp_pt)./mldepth_pt;
           array.salt(t,pt)=sum(salt_pt)./mldepth_pt;
           array.mldepth(t,pt)=mldepth_pt;
       end
       
    end
    close(h)

end
