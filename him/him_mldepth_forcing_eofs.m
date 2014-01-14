load metrics
hindcast.geolon=metrics.geolon.data;
hindcast.geolat=metrics.geolat.data;

%%
hindpath='/ltraid1/ashao/HIM/hyak_store/HINDCAST/';
begyear=1948;
endyear=2007;

nyears=length(begyear:endyear);
hindcast.mldepth=zeros(nyears*12,210,360);
hindcast.temp=zeros(nyears*12,210,360);
hindcast.salt=zeros(nyears*12,210,360);
hindcast.temp=zeros(nyears*12,210,360);
hindcast.ustar=zeros(nyears*12,210,360);
hindcast.heatflux=zeros(nyears*12,210,360);
hindcast.time=zeros(nyears*12,1);
n=0;
for year=begyear:endyear
    n=n+1;
    fprintf('Extracting year %d\n',year);
    
    oceanfile=strcat(hindpath,sprintf('ocean_month.%d.nc',year));
    begidx=(n-1)*12+1;
    endidx=n*12;
    idx=begidx:endidx;
        fprintf('Time...')
    hindcast.time(idx)=nc_varget( ...
        oceanfile,'Time');    
    fprintf('Mixed Layer Depth...')
    hindcast.mldepth(idx,:,:)=(sum(nc_varget( ...
        oceanfile,'h',[0 0 0 0],[-1 2 -1 -1]),2));    
    fprintf('SST...')
    hindcast.temp(idx,:,:)=(nc_varget(...
        oceanfile,'temp',[0 0 0 0],[-1 1 -1 -1]));
    fprintf('SSS...')
    hindcast.salt(idx,:,:)=(nc_varget(...
        oceanfile,'salt',[0 0 0 0],[-1 1 -1 -1]));
    fprintf('Ustar...');
    hindcast.ustar(idx,:,:)=(nc_varget(...
        oceanfile,'ustar'));
    fprintf('Heat flux\n');
    hindcast.heatflux(idx,:,:)=(nc_varget(...
        oceanfile,'LwLatSens'));
    
end