function [ ] = off_timeseries(tspath, ncfile, smon, emon, field, plotfun, options, nplotx, nploty)
% OFF_TIMESERIES Creates a number of given plot type for a given time range
% and field
% INPUT:
%   tsfile: path to tsfile directory
%   ncfile: path/fname of input netcdf
%   smon: start month
%   emon: end month
%   field: name of field to use
%   plotfun: handle to function control processing
%   options: structure handling configuration of plotfun
%   nplotx: number of subplots (horizontal)
%   nploty: number of subplots (vertical)

figure
nplots=(emon-smon)+1;
counter=1;
for i=smon:emon
    array=off_proc(tspath,ncfile,i);
    subplot(nplotx,nploty,counter)
    title([ field ' in month ' num2str(i)])
    plotfun(array,field,options)
    counter=counter+1;
end
    