function [ ] = off_depthcontour(array, field, options)
% OFF_DEPTHCONTOUR: Contour plot of field by depth along geographic lat/lon
% INPUT:
%   array: structure from off_proc
%   field: name of field to plot
%   options: structure containing fields
%       type (string): 'lat' or 'lon' for lat/lon slice
%       geoslice: value of lat or lon to make section
%       contourlevel (opt): Levels to contour the data
%       axis (opt): Range of plot
%       caxis (opt): Colorbar axis

data=getfield(array,field);

[nlay nlat nlon]=size(data);

if ~isfield(options,'contourlevel')
    interval=(max(data(:))-min(data(:)))/25;
    options.contourlevel=min(data(:)):interval:max(data(:));
end

if ~isfield(options,'axis')
    options.axis='auto';
end

if ~isfield(options,'caxis')
    options.caxis='auto';
end

if strcmpi(options.type,'lat')
    [null latidx]=min(abs(array.lath-options.geoslice));   
    depthtab=squeeze(array.depth(:,latidx,:));
    lontab=ones(nlay,1)*array.lonh';
    datatab=squeeze(data(:,latidx,:));
    [c h]=contour(lontab,-depthtab,datatab,options.contourlevel)
    caxis(options.caxis);
    axis(options.axis);
    clabel(c,h);
    colorbar
    title([field ' at latitude ' num2str(options.geoslice) ' for month ' num2str(array.month)])    
end

if strcmpi(options.type,'lon')
    [null lonidx]=min(abs(array.lonh-options.geoslice));
%     lonidx
    depthtab=squeeze(array.depth(:,:,lonidx));
    lattab=ones(nlay,1)*array.lath';
    datatab=squeeze(data(:,:,lonidx));
    colormap('jet')
    [c h]=contourf(lattab,-depthtab,datatab,options.contourlevel);
    caxis(options.caxis);
    axis(options.axis);
    clabel(c,h,'Color','w');
    colorbar
    title([field ' at latitude ' num2str(options.geoslice) ' for month ' num2str(array.month)])    
end

end