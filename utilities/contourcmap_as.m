function hndlout = contourcmap(varargin)
%CONTOURCMAP Contour colormap and colorbar for current axes
% 
%   CONTOURCMAP(CMAPSTR) updates the figure's colormap for the current axes
%   with the colormap specified by the string, CMAPSTR. Valid entries for
%   CMAPSTR include  'pink', 'hsv', 'jet', or any similar MATLAB colormap
%   function. If the axes contains Mapping Toolbox contour objects, the
%   resultant colormap contains the same number of colors as the original
%   colormap. Otherwise, the resultant colormap contains 10 colors.
% 
%   CONTOURCMAP(CMAPSTR, CDELTA) updates the figure's colormap with colors
%   varying according to CDELTA. If CDELTA is a scalar, it represents a
%   step size, and colors are generated at multiples of CDELTA. If CDELTA
%   is a vector of evenly spaced values, colors are generated at those
%   values; otherwise an error is issued. If the axes contains Mapping
%   Toolbox contour objects, the value of CDELTA is ignored.
%
%   CONTOURCMAP(..., PARAM1, VAL1, PARAM2, VAL2,...) allows you to add a
%   colorbar and control the colorbar's properties. Parameter names can be
%   abbreviated and are case-insensitive. See the table below for a list of
%   available parameters.
%
%   Optional Parameters
%   -------------------
%  
%   'Colorbar'          String with values 'on' or 'off' specifying 
%                       whether a colorbar is present, 'on', or absent
%                       from the axes, 'off'.
%
%   'Location'          String specifying the location of the colorbar.
%                       Permissible values are 'vertical' (the default), 
%                       'horizontal', or 'none'.
%
%   'ColorAlignment'    String specifying the alignment of the labels in
%                       the colorbar. Permissible values are 'center',
%                       where the labels are centered on the color bands or
%                       'ends' where labels are centered on the color
%                       breaks. If the axes contains Mapping Toolbox
%                       contour objects, the ColorAlignment will be set
%                       automatically to 'center' for contour lines and
%                       'ends' for filled contours, and cannot be modified.
%
%   'SourceObject'      Handle of the graphics object which is used to 
%                       determine the color limits for the colormap. The
%                       SourceObject value is the handle of a currently
%                       displayed object. If omitted, gca is used.
%
%   'TitleString'       String specifying the title of the colorbar axes.
%
%   'XLabelString'      String specifying the X label of the colorbar axes.
%
%   'YLabelString'      String specifying the Y label of the colorbar axes.
%
%   'ZLabelString'      String specifying the Z label of the colorbar axes.
%                        
%   In addition, properties and values that can be applied to the title and
%   labels of the colorbar axes are valid.
% 
%   h = CONTOURCMAP(...) returns a handle to the colorbar axes.
% 
%   Examples
%   --------
%   load topo
%   R = georasterref('RasterSize', size(topo), ...
%       'Latlim', [-90 90], 'Lonlim', [0 360]);
%   figure('Color','white')
%   worldmap(topo, R)
%   contourfm(topo, R);
%   contourcmap('jet', 'Colorbar', 'on', ...
%      'Location', 'horizontal', ...
%      'TitleString', 'Contour Intervals in Meters');
%
%   load topo
%   load coast
%   R = georasterref('RasterSize', size(topo), ...
%       'Latlim', [-90 90], 'Lonlim', [0 360]);
%   figure('Color','white')
%   worldmap(topo, R)
%   geoshow(topo, R, 'DisplayType', 'texturemap');
%   contourcmap('summer', 2000, 'Colorbar', 'on', ...
%      'Location', 'horizontal', ...
%      'TitleString', 'Contour Intervals in Meters');
%   geoshow(lat, long, 'Color', 'black')
%
%   See also CLABELM, CLEGENDM, COLORMAP, CONTOURCBAR, CONTOURFM, CONTOURM.

% Copyright 1996-2011 The MathWorks, Inc.
% $Revision: 1.6.4.12 $  $Date: 2011/05/17 02:19:09 $

% Parse required parameters.
error(nargchk(1, Inf, nargin, 'struct'))
if isnumeric(varargin{1})
    error(nargchk(2, Inf, nargin, 'struct'))
end
[cdelta, cmapfunc, params] = parseParams(varargin{:});

% Check for parameter/value pairs.
[h, isContourHandle, cbarProps] = parseOptionalParams(params);

% Construct (or delete) a contourcmap colorbar.
if isContourHandle
    hCbarAxes = contourcmapContour(h, cdelta, cmapfunc, cbarProps);
else
    hCbarAxes = contourcmapSurface(h, cdelta, cmapfunc, cbarProps);
end

% Return handle output if requested.
if nargout > 0;
    hndlout = hCbarAxes;
end

end

%================== General Parsing and Validation ========================

function [cdelta, cmapfunc, params] = parseParams(varargin)

% CMAPSTR
if isnumeric(varargin{1})
    cmapstr = varargin{2};
    pos = [2,1];
    varargin(2) = [];
else
    cmapstr = varargin{1};
    varargin(1) = [];
    pos = [1,2];
end

% Validate CMAPSTR
validateattributes(cmapstr, {'char'}, {'nonempty'}, ...
    mfilename, 'CMAPSTR', pos(1));
assert(exist(cmapstr,'file') == 2, message(...
    'map:contourcmap:invalidParam','CMAPSTR'))
cmapfunc = str2func(cmapstr);

% CDELTA
if ~isempty(varargin) && isnumeric(varargin{1})
    cdelta = varargin{1};
    validateattributes(cdelta, {'numeric'}, {'real','nonempty'}, ...
        mfilename, 'CDELTA', pos(2));
    params = varargin(2:end);
else
    cdelta = [];
    params = varargin;
end

end

%--------------------------------------------------------------------------

function [h, isContourHandle, cbarProps] = parseOptionalParams(params)
% Parse the parameter/value pairs.

% Validate the optional parameter/value pairs.
cbarProps = validateOptionalParams(params);

% Set default colorbar location.
if isequal(cbarProps.colorbar, 'off')
    cbarProps.location = 'none';
end

% Find the source handle and determine if it is a GeoContourGroupHandle.
[h, isContourHandle] = findSourceObject(cbarProps.sourceobject);

end

%--------------------------------------------------------------------------

function cbarProps = validateOptionalParams(params)
% Validate the optional parameter/value pairs.

assert(mod(numel(params),2) == 0, message('map:contourcmap:invalidParams'))

% Assign default values.
cbarProps = struct( ...
    'titlestring', '', ...
    'xlabelstring', '', ...
    'ylabelstring', '', ...
    'zlabelstring', '', ...
    'titleParams', '', ...
    'location', 'vertical', ...
    'coloralignment', '', ...
    'colorbar', 'off', ...
    'sourceobject', []);

% Assign values.
if ~isempty(params)
    todelete = false(size(params));
    for i=1:2:length(params)
        params{i} = canonicalProps(params{i});
        switch params{i}
            case 'location'
                cbarProps.(params{i}) = validatestring(params{i+1},  ...
                    {'none','vertical','horizontal'}, ...
                    mfilename, '''Location''');
                todelete([i i+1]) = true;
                
            case 'coloralignment'
                cbarProps.(params{i}) = validatestring(params{i+1},...
                    {'center','ends'}, mfilename, '''ColorAlignment''');
                todelete([i i+1]) = true;
                
            case 'colorbar'
                cbarProps.(params{i}) = validatestring(params{i+1}, ...
                    {'on','off'}, mfilename, '''Colorbar''');
                todelete([i i+1]) = true;
                
            case {'xlabelstring', 'ylabelstring', ...
                  'zlabelstring', 'titlestring', ...
                  'sourceobject'}
                cbarProps.(params{i}) = params{i+1};
                todelete([i i+1]) = true;
        end
    end
    params(todelete) = [];
end
cbarProps.titleParams = params;

if isempty(cbarProps.sourceobject)
    cbarProps.sourceobject = gca;
end

end

%--------------------------------------------------------------------------

function out = canonicalProps(in)
% Expand property names to canonical names

try
    out = lower(validatestring(in, ...
        {'Location', 'ColorAlignment', 'SourceObject', 'TitleString',...
        'XLabelString', 'YLabelString', 'ZLabelString', 'Colorbar'}, ...
        mfilename));
catch e
    if strcmp(e.identifier,'MATLAB:contourcmap:unrecognizedStringChoice') ...
            && ~isempty(e.cause)
        rethrow(e)
    else
        out = in;
    end
end

end

%--------------------------------------------------------------------------

function [h, isContourHandle] = findSourceObject(h)
% Find the source handle from H and determine if the source handle is a
% GeoContourGroup handle.

assert(ishghandle(h), message('map:contourcmap:invalidHandle','SourceObject'))

isContourHandle = false;

g = findobj(h, 'Type', 'hggroup');
if ~isempty(g)
    k = 1;
    while ~isContourHandle && (k <= numel(g))
        if isappdata(g(k), 'mapgraph')
            mapgraph = getappdata(g(k), 'mapgraph');
            if isa(mapgraph, 'internal.mapgraph.ContourGroup')
                isContourHandle = true;
                h = g(k);
            end
        end
        k = k + 1;
    end
end


end

%====================== Contour-Specific Functions ========================

function hCbarAxes = contourcmapContour(h, cdelta, cmapfunc, cbarProps)
% Construct (or delete) a contourcmap colorbar for the case in which h
% contains a contour object constructed by contourm.

% Determine the axes containing the data. (If h is a handle to an axes,
% then ax will simply equal h.)
ax = ancestor(h,'axes');

% Overwrite h with a handle to the mapgraph.ContourGroup object.
h = getappdata(h,'mapgraph');

% Validate the special contour parameters, assigning a value to the
% coloralignment parameter, if not set.
cbarProps.coloralignment ...
    = validateContourParams(h.Fill, cdelta, cbarProps.coloralignment);

% Update the figure's colormap while preserving the number of colors.
% This change should trigger listeners set up by contourm that will in
% turn update both h.LineColormap and h.FillColormap.
f = ancestor(h.HGGroup, 'figure');
cmap = cmapfunc(size(get(f,'Colormap'),1));
set(f,'Colormap',cmap)

% Return levels that fit in-range and a P-by-3 array of RGB colors.
[clevels, colors] = colorbarContourLevelsAndColors(h, strcmp(h.Fill,'on'));

% Create or delete the colorbar.
if isequal(cbarProps.colorbar, 'on')
    % Create the colorbar.
    hCbarAxes = createColorbar(ax, clevels, cbarProps, @(hAxes, location) ...
        contourImage(hAxes, strcmp(location,'vertical'), colors));
    listenToContourColormap(h, hCbarAxes)
else
    % Delete the colorbar, if found.
    deleteColorbar(ax);
    hCbarAxes = [];
end

end

%--------------------------------------------------------------

function contourImage(hAxes, colorbarIsVertical, colors)
% Create an image with the correct location, shape, and (RGB) cdata in the
% contour case.

n = size(colors,1);
if colorbarIsVertical
    xLimits = [0 1];
    yLimits = [1 n];
    % n-by-1 RGB image
    cdata = reshape(colors, [n 1 3]);
else
    xLimits = [1 n];
    yLimits = [0 1];
    % 1-by-n RGB image
    cdata = reshape(colors, [1 n 3]);
end
contourf(xLimits, yLimits, cdata, ...
    'Parent',hAxes,'Tag','CONTOURCMAP','DeleteFcn',@deleteim);
end

%--------------------------------------------------------------------------

function coloralignment ...
    = validateContourParams(fill, cdelta, coloralignment)

% Display a warning if CDELTA has been set.
if ~isempty(cdelta)
    warning(message('map:contourcmap:ignoringCDELTA','CDELTA'))
end

% Assign a value to the coloralignment parameter, if not set.
if isempty(coloralignment)
    if isequal(fill,'on')
        coloralignment = 'ends';
    else
        coloralignment = 'center';
    end
end

% Display a warning if coloralignment is not set properly.
if isequal(fill, 'on') && ~isequal(coloralignment, 'ends')
    warning(message('map:contourcmap:ignoringCenterColorAlignment', ...
        'ColorAlignment', 'ends'))
    coloralignment = 'ends';
    
elseif isequal(fill, 'off') && ~isequal(coloralignment, 'center')
    warning(message('map:contourcmap:ignoringEndsColorAlignment', ...
        'ColorAlignment', 'center'))
    coloralignment = 'center';
end

end

%--------------------------------------------------------------------------

function listenToContourColormap(h, hndl)
% Given the handle to an internal.mapgraph.ContourGroup object H and the
% handle to a CONTOURCMAP axes HNDL, add a listener to replace the image in
% the HNDL axes whenever the LineColormap or FillColor map property of H is
% reset (and update/restore axes properties as needed). Also set the
% DeleteFcn of HNDL to clean up by deleting the listener.
%
% Without the DeleteFcn, the sequence of (1) deleting the CONTOURCMAP axes
% (without deleting H) and then (2) changing a colormap property of H could
% result in an "Invalid handle" error.

fill = strcmp(h.Fill,'on');
if fill
    hListener = addlistener(h, 'FillColormapUpdate', @updateColors);
else
    hListener = addlistener(h, 'LineColormapUpdate', @updateColors);
end

set(hndl,'DeleteFcn',@deleteColormapListener)

    %----------------- Nested callback functions --------------------------
    
    function updateColors(~,~)
        % Construct a new colorbar image with updated colors, 
        % and update colorbar axes properties as required.
        
        % Return levels that fit in-range and a P-by-3 array of RGB colors.
        fill = strcmp(h.Fill,'on');
        [clevels, colors] = colorbarContourLevelsAndColors(h, fill);
        
        % Find the current colorbar image object, determine its
        % orientation, and delete it after disabling its DeleteFcn.
        hImage = findobj(hndl,'Type','image');
        colorbarIsVertical = (size(get(hImage,'CData'),2) == 1);
        set(hImage,'DeleteFcn',[])
        delete(hImage)
        
        % Replace the colorbar image.
        contourImage(hndl, colorbarIsVertical, colors)
        
        % Reset/restore colobar axes properties as needed.
        if fill
            tickloc = getTickLocation('ends', numel(clevels));
        else
            tickloc = getTickLocation('center', numel(clevels));
        end
        
        if colorbarIsVertical
            set(hndl, 'YTick', tickloc, 'YTickLabel', clevels, ...
                'XTick', [], 'YDir', 'normal', 'YAxisLocation', 'right')
        else
            set(hndl, 'XTick', tickloc, 'XTickLabel', clevels, ...
                'YTick', [], 'YDir', 'normal')
        end
    end

    function deleteColormapListener(~,~)
        if isvalid(hListener)
            delete(hListener)
        end
    end
end

%====================== Surface-Specific Functions ========================

function hCbarAxes = contourcmapSurface(h, cdelta, cmapfunc, cbarProps)
% Construct (or delete) a contourcmap colorbar for the case in which h
% contains a surface object.

% Create a contour colormap for the object in H and return the colors
% and levels required for the colorbar.

% Determine the axes containing the data. (If h is a handle to an axes,
% then ax will simply equal h.)
ax = ancestor(h,'axes');

% Assign a value to the coloralignment parameter, if not set.
if isempty(cbarProps.coloralignment)
    cbarProps.coloralignment = 'ends';
end

% Validate the color limits.
climits = validateColorLimits(h);

% Compute the color levels for labels with nice increments.
clevels = computeColorLevels(cdelta, climits, ax);

% Update the colormap based on the number of colors.
f = ancestor(ax, 'figure');
updateColormap(f, cmapfunc, cbarProps.coloralignment, clevels);

% Obtain the color index and the number of colors.
cindex = 1:size(colormap,1);
numcolors = size(cindex,1);

% Create or delete the colorbar.
if isequal(cbarProps.colorbar, 'on')
    % Create the colorbar.
    hCbarAxes = createColorbar(ax, clevels, cbarProps, @createImage);
else
    % Delete the colorbar, if found.
    deleteColorbar(ax);
    hCbarAxes = [];
end

    %--------------------------------------------------------------
    % Nested function to create an image with the correct location,
    % shape, and cdata in the non-contour case.
    function createImage(hAxes, location)
        if strcmp(location,'vertical')
            xLimits = [0 1];
            yLimits = [1 numcolors];
            cindex = cindex';
        else
            xLimits = [1,numcolors];
            yLimits = [0 1];
        end
        image(xLimits, yLimits, cindex, ...
            'Parent',hAxes,'Tag','CONTOURCMAP','DeleteFcn',@deleteim);
    end

end

%--------------------------------------------------------------------------

function climits = validateColorLimits(h)
% Check that h is a handle and get color limits

if ishghandle(h,'axes')
    set(h,'CLimMode','auto')
    climits = get(h, 'CLim');
else
    try
        cdata = get(h,'CData');
    catch e %#ok<NASGU>
        error(message('map:contourcmap:invalidHandleObject'))
    end
    
    if isempty(cdata) || all(isnan(cdata(:)))
        error(message('map:contourcmap:invalidCDATA','CData'))
    end
    
    climits = [min(cdata(:)) max(cdata(:))];
    if diff(climits) == 0;
        error(message('map:contourcmap:invalidCDATAProperty','CData'))
    end
end

end

%--------------------------------------------------------------------------

function clevels = computeColorLevels(cdelta, climits, ax)
% Compute the color levels and reset CLim for the data axes ax.

if isscalar(cdelta)
	cmin = floor(climits(1));
	multfactor = floor(cmin/cdelta);
	cmin = cdelta*multfactor;
	cmax = ceil(climits(2));
	clevels = cmin:cdelta:cmax;
    if max(clevels) < cmax
        clevels = [clevels max(clevels)+cdelta];
    end
elseif isempty(cdelta)
    cmin = floor(climits(1));
    cmax = ceil(climits(2));
    cdelta = abs(diff(climits))/10;
    clevels = cmin:cdelta:cmax;
else
	% check to see that spacing is the same
	tolerance = eps;
    if abs(max(diff(diff(cdelta)))) > tolerance
        error(message('map:contourcmap:invalidCDATASize','CDELTA'))
    end
	clevels = cdelta;
end	

% round numbers less than epsilon to zero
clevels(abs(clevels) < eps) = 0;

% Reset the climits based on new increments
set(ax, 'CLim', [min(clevels) max(clevels)])

end

%--------------------------------------------------------------------------

function updateColormap(f, cmapfunc, coloralignment, clevels)
% Update the colormap based on the color alignment.

if isequal(coloralignment, 'ends')
    % 'ends'
    numberofcolors = length(clevels)-1;
else
    % 'center'
    numberofcolors = length(clevels);
end
cmap = cmapfunc(numberofcolors);
set(f,'Colormap',cmap);

end

%======================= Colorbar Construction ============================

function hAxes = createColorbar(hax, clevels, cbarProps, createImage)
% Create the colorbar.

hfig = ancestor(hax,'figure');

% Get the axes units and change them to normalized.
axesunits = get(hax,'Units');
set(hax, 'Units', 'normalized')

% Create the colorbar, if present.
deleteColorbar(hax);

% Create the structure for callbacks.
axesinfo.h = hax;
axesinfo.units = axesunits;
axesinfo.origpos = get(hax, 'Position');

switch cbarProps.location 
    case 'none'
        hAxes = [];
        
    case 'vertical'
        hAxes = createVerticalColorbar(axesinfo, hax, ...
            cbarProps.coloralignment, clevels, createImage);
        
    case 'horizontal'
        hAxes = createHorizontalColorbar(axesinfo, hax, ...
            cbarProps.coloralignment, clevels, createImage);
end

% Set delete function.
setAxesDeleteFcn(hAxes, axesinfo)

% Reset the axes units.
set(hax, 'Units', axesunits)

% Activate the initial axes.
set(hfig,'CurrentAxes',hax)

% Set text properties of colorbar axes and title if provided.
if ~isempty(hAxes)
    setColorbarProps(hAxes, cbarProps);
end

end

%--------------------------------------------------------------------------

function  hAxes = createVerticalColorbar( ...
    axesinfo, hax, coloralignment, clevels, createImage)

% Shrink length by 10 percent.
pos = axesinfo.origpos;
pos(3) = axesinfo.origpos(3)*0.90;
set(hax, 'Position', pos)

% Calculate the position of the colorbar axes.
len   = axesinfo.origpos(3)*0.05;
width = axesinfo.origpos(4);
axesPos = [axesinfo.origpos(1)+axesinfo.origpos(3)*0.95 axesinfo.origpos(2) len width];

% Compute the tick locations.
ytickloc = getTickLocation(coloralignment, numel(clevels));

% Create the colorbar axes.
hAxes = axes('Position', axesPos);

% Create the image.
createImage(hAxes,'vertical')

% Create the axes properties.
set(hAxes,...
    'YTick', ytickloc, 'YTickLabel', clevels, ...
    'XTick', [], 'YDir', 'normal', 'YAxisLocation', 'right');

end

%--------------------------------------------------------------------------

function hAxes = createHorizontalColorbar( ...
    axesinfo, hax, coloralignment, clevels, createImage)

% Shrink width by 10 percent.
pos = axesinfo.origpos;
pos(4) = axesinfo.origpos(4)*0.90;
pos(2) = axesinfo.origpos(2) + axesinfo.origpos(4)*0.10;
set(hax, 'Position', pos)

% Calculate the position of the colorbar axes.
width = axesinfo.origpos(4)*0.05;
len = axesinfo.origpos(3);
axesPos = [axesinfo.origpos(1) axesinfo.origpos(2) len width];

% Compute the tick locations.
xtickloc = getTickLocation(coloralignment, numel(clevels));

% Create the colorbar axes.
hAxes = axes('Position', axesPos);

% Create the image.
createImage(hAxes,'horizontal')

% Set the axes properties.
set(hAxes, ...
    'XTick', xtickloc, 'XTickLabel', clevels, ...
    'YTick', [], 'XDir', 'normal', 'XAxisLocation', 'bottom');

% Use appdata to keep track of the association between the colorbar axes
% with the data axes (its "peer").
setappdata(hAxes, 'peer', hax)

end

%--------------------------------------------------------------------------

function tickloc = getTickLocation(coloralignment, numticks)
% Calculate the locations for the tick marks.

if isequal(coloralignment, 'ends')
    numticks = numticks - 1;
end

switch coloralignment
    case 'ends'
        lowerlim = 1-0.5;
        upperlim = numticks+0.5;
    case 'center'
        lowerlim = 1;
        upperlim = numticks;
end
delta = 1;
tickloc = lowerlim:delta:upperlim;

end

%--------------------------------------------------------------------------

function setColorbarProps(hndl, cbarProps)
% Set the special properties of the colorbar axes.

set(get(hndl,'Title'), 'String',cbarProps.titlestring);
set(get(hndl,'Xlabel'),'String',cbarProps.xlabelstring);
set(get(hndl,'Ylabel'),'String',cbarProps.ylabelstring);
set(get(hndl,'Zlabel'),'String',cbarProps.zlabelstring);

if ~isempty(cbarProps.titleParams)
    set(get(hndl,'Title'),cbarProps.titleParams{:});
    set(hndl,cbarProps.titleParams{:});
end

end

%--------------------------------------------------------------------------

function deleteColorbar(ax)
% Delete the colorbar corresponding to data axes ax.
    
% Find the colorbar image.
f = ancestor(ax,'figure');
hImage = findobj(f,'tag','CONTOURCMAP','type','image');

for k = 1:numel(hImage)
    h = hImage(k);
    if ishghandle(hImage)
        % See if data axes ax is the "peer' of the images axes ancestor.
        axColorbar = ancestor(h, 'axes');
        peer = getappdata(axColorbar, 'peer');
        if isequal(peer, ax)
            % Delete the colorbar axes and image.
            delete(axColorbar);
        end
    end
end 

end

%-------------------------------------------------------------------------

function setAxesDeleteFcn(ax, axesinfo)

set(ax, 'DeleteFcn', @deleteax)

    function deleteax(~, ~)
        % DeleteFcn callback for the axes holding the colorbar.
        
        if ~isempty(axesinfo) && isfield(axesinfo,'h')
            hDataAxes = axesinfo.h;
            if ishghandle(hDataAxes,'axes');
                % Restore the data axes to its original position.
                set(hDataAxes,'Units','normalized')
                set(hDataAxes,'Position',axesinfo.origpos)
                set(hDataAxes,'Units',axesinfo.units)
            end
        end
    end
end

%-------------------------------------------------------------------------

function deleteim(h, ~)
% DeleteFcn callback for the image object used in the colorbar.

% Double check that h is a 'contourcmap' image, then delete parent axes.
% Let the axes' delete function do the actual work.
if ishghandle(h,'image') && strcmpi(get(h,'Tag'),'contourcmap')
    ax = get(h,'Parent');
    delete(ax)
end

end
