function [ atmcon ] = trac_atmospheric(tracer, time, lat)
% TRAC_ATMOSPHERIC: Calculate atmospheric concentration for a given time
% and latitude
% INPUT:
%   tracer (struct) fields:
%      Nval: northern hemisphere values
%      Sval: southern hemisphere values
%      year: year corresponding to the values
%   time: year (or year + fraction of year) to find atmospheric
%       concentration
%   latitude: latitude of interest

% Boundaries of Northern/Southern hemisphere
% Nlat=10.2;
% Slat=-10.2;
% 
% year=floor(time);
% yearfrac=time-year;
% 
% [null curr_idx]=min(abs(tracer.year-time));
% % curr_idx
% prev_idx=curr_idx-1;
% next_idx=curr_idx+1;
% 
% if lat>Nlat
%     prev_val=tracer.Nval(prev_idx);
%     curr_val=tracer.Nval(curr_idx);
%     next_val=tracer.Nval(next_idx);
% elseif lat<Slat
%     prev_val=tracer.Sval(prev_idx);
%     curr_val=tracer.Sval(curr_idx);
%     next_val=tracer.Sval(next_idx);
% else
%     prev_val=interp1([Slat, Nlat],[tracer.Sval(prev_idx) tracer.Nval(prev_idx)],lat);
%     curr_val=interp1([Slat, Nlat],[tracer.Sval(curr_idx) tracer.Nval(curr_idx)],lat);
%     next_val=interp1([Slat, Nlat],[tracer.Sval(next_idx) tracer.Nval(next_idx)],lat);
% end
% 
% atmcons=[prev_val curr_val next_val];
% atmyear=[tracer.year([prev_idx curr_idx next_idx])];
% tracatm=interp1(atmyear,atmcons,time);

    nlat=10.2;
    slat=-10.2;
    nidx=find(lat>nlat);
    sidx=find(lat<slat);
    eidx=find(lat<=nlat & lat>=slat);
    
    if time<tracer.year(length(tracer.year))
        
        nval=interp1(tracer.year,tracer.Nval,time);
        sval=interp1(tracer.year,tracer.Sval,time);
        
        atmcon=zeros(size(lat));
    else
        nval=tracer.Nval(length(tracer.year));
        sval=tracer.Sval(length(tracer.year));
    end
    atmcon(nidx)=nval;
    atmcon(sidx)=sval;
    atmcon(eidx)=interp1([nlat slat],[nval sval],lat(eidx));
end