function array = him_diag_pme(ncfile,metricfile,pmefield)

Ah_grid=nc_varget(metricfile,'Ah');
Ah_vec=makevec(Ah_grid);
array.geolat=nc_varget(metricfile,'geolat');
array.geolon=nc_varget(metricfile,'geolon');
array.time=nc_varget(ncfile,'Time');
nmon=length(array.time);
wet_grid=nc_varget(metricfile,'wet');
wet_vec=makevec(wet_grid);
dry_vec=~wet_vec;
dry_grid=~wet_grid;
Ah_vec(dry_vec)=NaN;
Ah_grid(dry_grid)=NaN;
nyears=floor(nmon/12);
array.global_mean_net_pme=zeros(nyears,1);
array.net_annual_pme=zeros([nyears size(array.geolon)]);
array.global_net_pme=zeros(nyears,1);
seconds_in_day=24*60*60;
seconds_in_month=eomday(1900,1:12)*seconds_in_day;
seconds_in_month_grid(:,1,1)=seconds_in_month;
seconds_in_month_grid=repmat(seconds_in_month_grid,[1 size(array.geolon)]);

for t=0:nyears-1
    start_3d=[12*t 0 0];
    count_3d=[12 -1 -1];
    pme=nc_varget(ncfile,pmefield,start_3d,count_3d).*seconds_in_month_grid;
    pme(:,dry_grid)=NaN;
%     array.global_net_pme(t+1)=sum(sum(pme,2),3);
    array.net_annual_pme(t+1,:,:)=sum(pme,1);
    array.global_mean_net_pme(t+1,:)=weight_mean(makevec(sum(pme,1)),Ah_vec);
    array.global_net_pme(t+1)=nansum(nansum(nansum(pme)));
end
    
    array.net_annual_pme(:,dry_grid)=NaN;    
    
    
end

