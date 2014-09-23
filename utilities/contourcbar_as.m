function cbar = contourcbar(varargin)
%contourcbar Color bar for filled contour map display
%
%   H = contourcbar(...) creates a colorbar associated with a filled
%   contour display created with CONTOURFM, CONTOURM, CONTOUR3M, or
%   GEOSHOW. It supports the same syntax and usage options as the MATLAB
%   function COLORBAR. If a 'peer' axes is specified when calling
%   contourcbar, it should be a map axes containing an object created using
%   one of the Mapping Toolbox functions listed previously. Otherwise the
%   current axes should contain such an object.
%
%   If a Mapping Toolbox contour object is present, then the colorbar is
%   filled with solid blocks of color which bound each other at the
%   contour levels used in the plot. Thus, the contour levels bounding a
%   fill polygon of a given color can be inferred graphically by inspecting
%   the upper and lower limits of the corresponding block in the colorbar.
%   In the absence of a Mapping Toolbox contour object an ordinary color
%   bar is created.
%
%   If multiple Mapping Toolbox contour objects are present in the same
%   axes, then the levels used to divide the colorbar into blocks will
%   correspond to the first contour object that is found. This situation
%   could occur when a larger data set is broken up into multiple grid
%   tiles, for example, but as long the tiles all use the same contour
%   level list, the colorbar will correctly represent them all.
%
%   Examples
%   --------
%   % Topography of North America
%   figure('Color','white')
%   worldmap('north america')
%   load topo
%   R = georasterref('RasterSize',[180 360], ...
%         'LatitudeLimits',[-90 90],'LongitudeLimits',[0 360]);
%   contourfm(topo, R, -7000:1000:3000)
%   caxis([-8000 4000])
%   contourcbar
%
%   % Geoid with non-uniform levels
%   figure('Color','white')
%   ax = worldmap('world');
%   setm(gca,'MLabelParallel',-90)
%   setm(gca,'MLabelLocation',90)
%   load geoid
%   R = georasterref('RasterSize',[180 360], ...
%         'LatitudeLimits',[-90 90],'LongitudeLimits', [0 360]);
%   levels = [-70 -40 -20 -10 0 10 20 40 70];
%   geoshow(geoid, R, 'DisplayType', 'contour',...
%       'LevelList',levels,'Fill','on','LineColor','black')
%   coast = load('coast.mat');
%   geoshow(coast.lat, coast.long, 'Color', 'white', 'LineWidth', 1.5)
%   cb = contourcbar('peer',ax,'Location','southoutside');
%   caxis([-110 90])
%   colormap(hsv)
%   set(get(cb,'XLabel'),'String','Geoid Undulation in Meters')
%
%   See also CLEGENDM, COLORBAR, CONTOURFM

% Copyright 2011-2013 The MathWorks, Inc.

if ((nargin == 1) || ((nargin == 2) ...
    && isscalar(varargin{1}) && isColorbar(varargin{1}))) ...
    && offHideOrDelete(varargin{end})
    % colorbar('off')
    % colorbar('hide')
    % colorbar('delete')
    % colorbar(cbar_handle,'off')
    % colorbar(cbar_handle,'hide')
    % colorbar(cbar_handle,'delete')
    
    % In this case, colorbar will error if called with an output argument,
    % and this function should also.
    nargoutchk(0,0)
    colorbar(varargin{:});
else
    % Construct a standard colorbar, leveraging the
    % input validation provide by the colorbar function.
    cb = colorbar(varargin{:});
    
    h = findMapgraphContourGroup(varargin);
    if ~isempty(h)
        % There's a geographic contour group in the map axes, so modify the
        % colorbar to include only the discrete colors that are used for
        % fill polygons in the contour display.
        if matlab.graphics.internal.isGraphicsVersion1
            [xLimits,yLimits] = updateColorbarWithFillColors1(cb,h);
            setUpListeners(cb,h,xLimits,yLimits)
        else
            updateColorbarWithFillColors(cb,h)
            addlistener(h,'FillColormapUpdate', ...
                @(~,~) updateColorbarWithFillColors(cb,h));
        end
    end
    
    % Leave output cbar assigned unless the function is called with an
    % output argument.
    if nargout > 0
        cbar = cb;
    end
end

end

%------------------------------------------------------------------------

function tf = isColorbar(h)
% True if scalar graphics object with handle h is a colorbar

if matlab.graphics.internal.isGraphicsVersion1
    tf = ishghandle(h,'axes') && strcmpi(get(h,'Tag'),'colorbar');
else
    tf = ishghandle(h,'colorbar');
end

end

%------------------------------------------------------------------------

function tf = offHideOrDelete(str)
% Return true if STR matches 'off', 'hide', or 'delete'.

try
    validatestring(str, {'off','hide','delete'}, 'offHideOrDelete');
    tf = true;
catch e
    if strcmp(e.identifier,'MATLAB:offHideOrDelete:unrecognizedStringChoice')
        tf = false;
    else
        rethrow(e) % Re-throw any unexpected error
    end
end

end

%------------------------------------------------------------------------

function h = findMapgraphContourGroup(args)
% Starting from the "peer" axes if designated in the input argument cell
% array ARGS, or gca otherwise, locate an internal.mapgraph.ContourGroup
% object H, if present. Otherwise return [].

% Locate the "peer" axes object. (The colorbar function itself requires
% that 'peer' be spelled out in full, but the case is arbitrary.)
p = find(strcmpi('Peer',args));
if any(p)
    % At least one peer axes is designated. Use the last one 'Peer'
    % appears twice in the argument list.
    p = p(end);
    peer = args{p + 1};
else
    peer = gca;
end

% Try to find an internal.mapgraph.ContourGroup in the peer axes.
hg = findobj(peer,'Type','hggroup');
isContourGroup = false(size(hg));
for k = 1:numel(hg)
    if isappdata(hg(k),'mapgraph')
        g = getappdata(hg(k),'mapgraph');
        isContourGroup(k) = isa(g,'internal.mapgraph.ContourGroup');
    else
        isContourGroup(k) = false;
    end
end

% If any objects are found to have associated appdata containing an
% internal.mapgraph.ContourGroup, re-assign h to the first one of these.
if any(isContourGroup)
    hg = hg(isContourGroup);
    h = getappdata(hg(1),'mapgraph');
else
    h = [];
end

end

%------------------------------------------------------------------------

function [xLimits, yLimits] = updateColorbarWithFillColors1(cb,h)
% Update the XLim (horizontal bar) or YLim (vertical bar) property of the
% colorbar axes, and the XData (horizontal bar) or YData (vertical bar) of
% the colorbar image, and the CData of the colorbar region of solid color
% corresponding to each contour interval. These regions will not
% necessarily be all the same size, but will instead reflect the relative
% sizes of the contour intervals themselves.

[bigmap, limits] = deriveColorbarColormap(h);
ncolors = size(bigmap,1);

% Find the image that the colorbar function created -- it should be the one
% and only child of the colorbar axis, so the findobj line should not be
% needed.
c = get(cb,'Children');
if isscalar(c)
    hColorbarImage = c;
else
    hColorbarImage = findobj(c,'Type','image','Tag','TMW_COLORBAR');
end   

% Make the image invisible (rather than deleting it, to avoid
% triggering listeners set up by colorbar).
set(hColorbarImage,'Visible','off')

if isappdata(cb,'FillColorsImage')
    hFillColorsImage = getappdata(cb,'FillColorsImage');
else
    hFillColorsImage = image('Parent',cb, ...
        'XData',get(hColorbarImage,'XData'), ...
        'YData',get(hColorbarImage,'YData'), ...
        'CData',get(hColorbarImage,'CData'), ...
        'Tag','FillColors','HandleVisibility','off');
    setappdata(cb,'FillColorsImage',hFillColorsImage)
end

% Check the colorbar orientation.
colorbarIsHorizontal = any(strcmpi(get(cb,'Location'), ...
    {'north','south','northoutside','southoutside'}));

% Choose xdata or ydata values such that the fill colors image will fit
% exactly within the horizontal or vertical limits of the colorbar axes.
xydata = limits + [0.5 -0.5] * diff(limits) / ncolors;

if colorbarIsHorizontal
    set(cb,'Xlim',limits)
    set(cb,'YLim',[-1 2])
    set(hFillColorsImage,'XData',xydata,'YData',[0 1], ...
        'CData',reshape(bigmap,[1 ncolors 3]))
else
    set(cb,'XLim',[-1 2])
    set(cb,'Ylim',limits)
    set(hFillColorsImage,'Xdata',[0 1],'YData',xydata, ...
        'CData',reshape(bigmap,[ncolors 1 3]))
end

xLimits = get(cb,'XLim');
yLimits = get(cb,'YLim');

end

%------------------------------------------------------------------------

function setUpListeners(cb,h,xLimits,yLimits)

hListener = addlistener(h, ...
    'FillColormapUpdate', @(~,~) updateFillColors(cb,h));

xLimListener = addlistener(cb, 'XLim', 'PostSet', @(~,~) restore(cb));
yLimListener = addlistener(cb, 'YLim', 'PostSet', @(~,~) restore(cb));
addlistener(cb, 'Location', 'PostSet', @(~,~) updateFillColors(cb,h));

set(cb,'DeleteFcn',@deleteFillColormapListener)

    %----------------- Nested callback functions ------------------------
    
    function updateFillColors(cb,h)
        set(xLimListener,'Enable','off')
        set(yLimListener,'Enable','off')
        [xLimits, yLimits] = updateColorbarWithFillColors1(cb,h);
        set(xLimListener,'Enable','on')
        set(yLimListener,'Enable','on')
    end

    function restore(cb)
        % Restore image invisibility and colorbar axes limits.
        hColorbarImage = findobj(cb,'Type','image','Tag','TMW_COLORBAR');
        if ishghandle(hColorbarImage,'image')
            set(hColorbarImage,'Visible','off')
        end
        set(cb,'XLim',xLimits)
        set(cb,'YLim',yLimits)
    end

    function deleteFillColormapListener(~,~)
        if isvalid(hListener)
            delete(hListener)
        end
    end
end
 