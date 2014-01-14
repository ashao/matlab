function [ array ] = fu_extract( ncfile, datapath, latrange, lonrange, daterange, layers, fields )

    geolat=nc_varget([datapath 'metrics.nc'],'geolat');
    geolon=nc_varget([datapath 'metrics.nc'],'geolon');
    lath=nc_varget([datapath 'metrics.nc'],'lath');
    lonh=nc_varget([datapath 'metrics.nc'],'lonh');
    
    latidx=findrange(lath,min(latrange),max(latrange));
    lonidx=findrange(lonh,min(lonrange),max(lonrange));
    
    array.lath=single(lath(latidx));
    array.lonh=single(lonh(lonidx));
    array.geolat=single(geolat(latidx,lonidx));
    array.geolon=single(geolon(latidx,lonidx));
    
    ntime=length(min(daterange):max(daterange));
    nlayers=length(min(layers):max(layers));
    nlat=length(latidx);
    nlon=length(lonidx);
    start4d=[min(daterange) min(layers) min(latidx) min(lonidx)]-1;
    end4d=[ntime nlayers nlat nlon];

    for i=1:length(fields)
                
        fieldname=char(fields(i));
        disp(sprintf('Extracting field %s',fieldname))
        array.(fieldname)=single(nc_varget(ncfile,fieldname,start4d,end4d));
        
    end
    array.depth=cumsum(array.mn_h,2);
    array.depth(:,1,:,:)=0;
    
end

