function [ out ] = off_extrtimeseries(tspath, ncfile, tracpath, smon, emon, syear, fieldnames, lonrange, latrange)
% OFF_TIMESERIES Creates a number of given plot type for a given time range
% and field
% INPUT:
%   tsfile: path to tsfile directory
%   ncfile: path/fname of input netcdf
%   smon: start month
%   emon: end month
%   fieldnames.name: names of fields to use
%   nw_coord: lon, lat of northwest corner of geographic box;
%   se_coord: lon, lat of southeast corner of geographic box;

nfields=length(fieldnames);
out=struct;

for i=smon:emon
    disp([num2str(i) '/' num2str(length(smon:emon))])
    array=off_proc(tspath,ncfile,i,syear,tracpath);
    latidx=find(array.lath>=min(latrange) & array.lath<=max(latrange));
    lonidx=find(array.lonh>=min(lonrange) & array.lonh<=max(lonrange));
    
    out(i).lat=array.geolat(latidx,lonidx);
    out(i).lon=array.geolon(latidx,lonidx);
    out(i).mixedlyr=array.mn_h(4,latidx,lonidx);
    out(i).T=array.T(1,latidx,lonidx);
    out(i).S=array.S(1,latidx,lonidx);
    out(i).cfc11_kg=array.cfc11_kg(1,latidx,lonidx);
    out(i).cfc12_kg=array.cfc12_kg(1,latidx,lonidx);
    out(i).sf6_kg=array.sf6_kg(1,latidx,lonidx);
    out(i).cfc11_relsat=array.cfc11_relsat(1,latidx,lonidx);
    out(i).cfc12_relsat=array.cfc12_relsat(1,latidx,lonidx);
    out(i).sf6_relsat=array.sf6_relsat(1,latidx,lonidx);
    
end
    
    