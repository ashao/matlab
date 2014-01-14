
   
    datapath='/scratch/data/offtrac/input/';
    uhclim=nc_varget([datapath 'UH-clim.nc'],'UHCLIM',[0 0 0 0],[12 1 -1 -1]);
    vhclim=nc_varget([datapath 'UH-clim.nc'],'UHCLIM',[0 0 0 0],[12 1 -1 -1]);
    geolat=nc_varget([datapath 'metrics.nc'],'geolat');
    geolon=nc_varget([datapath 'metrics.nc'],'geolon');
    
    mtrans=double(sqrt((uhclim.^2+vhclim.^2)));
    
    
    