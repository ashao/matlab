function [ array ] = off_calcsla(datapath,startidx,endidx)
%     407,563
    ncfile=[filesep 'H.nc'];
    ntime=endidx-startidx;
    start=[startidx 0 0 0];
    count=[ntime -1 -1 -1];
    fillmatrix3d=zeros(ntime,210,360);
    SSH=squeeze(sum(nc_varget([datapath ncfile],'H',start,count),2));
    [ntime nlat nlon]=size(SSH);
    meanSSH=mean(SSH);
    SLA=SSH-repmat(meanSSH,[ntime,1,1]);
    SLA_detrend=fillmatrix3d;
    
    for i=1:nlat
        SLA_detrend(:,i,:)=detrend(squeeze(SLA(:,i,:)));
    end
    
    array.lath=nc_varget([datapath 'metrics.nc'],'lath');
    array.lonh=nc_varget([datapath 'metrics.nc'],'lonh');
    array.geolat=nc_varget([datapath 'metrics.nc'],'geolat');
    array.geolon=nc_varget([datapath 'metrics.nc'],'geolon');
    array.SLA=SLA_detrend;
    
end