function array = calcskew_longterm_sla(datapath)
% Calcuate skewness from aviso gridded files
% WARNING: Need ~18GB RAM

%%  Initialize Matrix
    files=dir([datapath '*.nc']);
    array.lat=nc_varget([datapath filesep files(1).name],'NbLatitudes');
    array.lon=nc_varget([datapath filesep files(1).name],'NbLongitudes');
    ntime=length(files);
    nlat=length(array.lat);
    nlon=length(array.lon);
    tempSLA=zeros([ntime nlon nlat],'single');    
    h=waitbar(0,sprintf('File %d/%d',0,ntime));
    %% Extract 
    for i=1:ntime
        waitbar(i/ntime, h,sprintf('File %d/%d',i,ntime))
        tempSLA(i,:,:)=single(nc_varget([datapath filesep files(i).name],'Grid_0001'));
    end       
    
    array.skewness=squeeze(skewness(tempSLA,1));
end
    
    
    