function [ array,year_ssh ] = him_diag_temp(ncfile,metricfile,hfield,tfield,sfield)
% Calculates mean SSH for each timestep in ncfile weighted by area of the
% grid cell
% ncfile: path to datafile
% metricfile: path to metrics.nc
% hfield: name of thickness field usually H or TCLIM
% tfield: name of temperature field usually T or TEMPCLIM

    array.time=nc_varget(ncfile,'Time'); 
    ntime=length(array.time);

%     Get metrics.nc information
    Ah=nc_varget(metricfile,'Ah');
    Ah=reshape(Ah,numel(Ah),1);
    array.geolon=nc_varget(metricfile,'geolon');
    array.geolat=nc_varget(metricfile,'geolat');
    wet=nc_varget(metricfile,'wet');
    wet=reshape(wet,numel(wet),1);
    wetidx=find(wet==1);
    
%     Allocate array for global mean ssh
    array.global_mean_temp=zeros(ntime,1);
    array.global_mean_salt=zeros(ntime,1);
    array.global_inst_ssh=zeros(ntime,1);
%     array.global_year_ssh=zeros(ntime,1);
    
    for t=1:ntime
        
       disp(sprintf('Month %d/%d',t,ntime))
       hnow=nc_varget(ncfile,hfield,[t-1 0 0 0],[1 -1 -1 -1]);
       temp=nc_varget(ncfile,tfield,[t-1 0 0 0],[1 -1 -1 -1]);       
       salt=nc_varget(ncfile,sfield,[t-1 0 0 0],[1 -1 -1 -1]);
%        year_ssh=nc_varget(ncfile,year_sshfield,[t-1 0 0],[1 -1 -1]);       
       
       temp=temp.*hnow;
       temp=sum(temp,1);
       depth=sum(hnow,1);
       temp=temp./depth;
       temp=reshape(temp,numel(temp),1);
       array.global_mean_temp(t)=weight_mean(temp(wetidx),Ah(wetidx));
       
       salt=salt.*hnow;
       salt=sum(salt,1);       
       salt=salt./depth;
       salt=reshape(salt,numel(salt),1);
       array.global_mean_salt(t)=weight_mean(salt(wetidx),Ah(wetidx));
       
%        depth=reshape(depth,numel(depth),1);       
%        array.global_inst_ssh(t)=weight_mean(depth(wetidx),Ah(wetidx));
       
%        year_ssh=reshape(year_ssh,numel(year_ssh),1);
%        array.global_year_ssh(t)=weight_mean(year_ssh(wetidx),Ah(wetidx));
    end
    

end
