% function [ ] = calc_skew_all(datapath,outpath)

    idxdir = '/ltraid2/ashao/uw-apl/data/skewjetmodel/indices/';

    aaofile = 'ftp://ftp.cpc.ncep.noaa.gov/cwlinks/norm.daily.aao.index.b790101.current.ascii';
    onifile= 'http://www.esrl.noaa.gov/psd/data/correlation/oni.data';                
    
    urlwrite(aaofile,[idxdir 'aao.data']);
    urlwrite(onifile,[idxdir 'oni.data']);
    
    aaotemp=importdata([idxdir 'aao.data']);
    aao.time=datenum(aaotemp(:,1),aaotemp(:,2),aaotemp(:,3));
    aao.val=aaotemp(:,4);
    %%
    idoni=fopen([idxdir 'oni.data']);
    range=fscanf(idoni,'%d %d',2);
    nyears = length(min(range):max(range));
    [monthgrid yeargrid]=meshgrid(1:12,min(range):max(range));
    cyear = min(range);
    
    idx = 0;
               
    onitemp=fscanf(idoni,'%f',[nyears 13]);        
%     tempstr=fgets(idoni)
    fclose(idoni);
% end