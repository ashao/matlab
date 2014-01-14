function array = off_corr_mldepth_relsat(array)

[ntime nlat nlon]=size(array.cfc11_relsat);
corr_array=zeros([nlat,nlon]);
minsat=zeros([ntime/12 nlat,nlon]);
maxml=zeros([ntime/12 nlat,nlon]);

for lat=1:nlat
    disp(sprintf('%d/%d',lat,nlat))
    for lon=1:nlon
        
        minsat(:,lat,lon)=binfunc(array,lat,lon,'cfc12_relsat',@min);
        maxml(:,lat,lon)=binfunc(array,lat,lon,'mldepth',@max);
        
        corrvals=corrcoef(squeeze(minsat(:,lat,lon)),squeeze(maxml(:,lat,lon)));
        corr_array(lat,lon)=corrvals(2,1);
    end
end

array.corr_array=corr_array;
array.minsat=minsat;
array.maxml=maxml;

end

    function corr_val = calc_corr(array,lat,lon,field1,field2)
        
        data1=array.(field1);
        data1=squeeze(data1(:,lat,lon));
        
        data2=array.(field2);
        data2=squeeze(data2(:,lat,lon));
        
        corr_val=corrcoef(data1,data2);
        corr_val=corr_val(2,1);
        
    end
    
    function binval = binfunc( array, lat, lon, field, func )
        data=array.(field);
        data=squeeze(data(:,lat,lon));        
        data=reshape( data, [numel(data)/12 12 ]);
        binval=func(data,[],2);                
    
    end