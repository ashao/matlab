function [ trac_sat ] = trac_calcsat( tracer, time, T, S, P, lat, lon )
%TRAC_CALCSAT Calculates saturation of a tracer gas in the global mixed layer
% INPUT:
%   tracer (struct) fields:
%      Nval: northern hemisphere values
%      Sval: southern hemisphere values
%      year: year corresponding to the values
%   time: year (or year + fraction of year) to find atmospheric
%       concentration
%   T (nlayers x lat x lon): Temperature in Celsius 
%   S (nlayers x lat x lon): Salinity in PSU
%   P (lat x lon): in atmospheres
%   lat: latitudes
%   lon: longitudes

if length(size(T))<3
    nlay=1;
    [nlat nlon]=size(T);
else
    [nlay nlat nlon]=size(T);
end
trac_sol=trac_calcsol(T+273.15,S,tracer.sol_coeffs);
trac_atmcon=zeros(nlay,nlat,nlon);

for i=1:nlay
    Pout(i,:,:)=P;
    for j=1:nlat
%         disp(lat(j))
        lat_atmcon=trac_atmospheric(tracer, time, lat(j));
        trac_atmcon(i,j,:)=ones(1,nlon)*lat_atmcon;
    end
end

Pout=squeeze(Pout);
trac_atmcon=squeeze(trac_atmcon);
trac_sat=trac_atmcon.*Pout.*trac_sol;
end

