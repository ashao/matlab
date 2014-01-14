function theAxis = subplot(nRows, nCols, plotId, varargin)
    %SUBPLOT Create axes in tiled positions.
    %   H = SUBPLOT(m,n,p), or SUBPLOT(mnp), breaks the Figure window
    %   into an m-by-n matrix of small axes, selects the p-th axes for
    %   the current plot, and returns the axes handle.  The axes are
    %   counted along the top row of the Figure window, then the second
    %   row, etc.  For example,
    %
    %       SUBPLOT(2,1,1), PLOT(income)
    %       SUBPLOT(2,1,2), PLOT(outgo)
    %
    %   plots income on the top half of the window and outgo on the
    %   bottom half. If the CurrentAxes is nested in a uipanel the
    %   panel is used as the parent for the subplot instead of the
    %   current figure.
    %
    %   SUBPLOT(m,n,p), if the axes already exists, makes it current.
    %   SUBPLOT(m,n,p,'replace'), if the axes already exists, deletes it and
    %   creates a new axes.
    %   SUBPLOT(m,n,p,'align') places the axes so that the plot boxes
    %   are aligned, but does not prevent the labels and ticks from
    %   overlapping.
    %   SUBPLOT(m,n,P), where P is a vector, specifies an axes position
    %   that covers all the subplot positions listed in P.
    %   SUBPLOT(H), where H is an axes handle, is another way of making
    %   an axes current for subsequent plotting commands.
    %
    %   SUBPLOT('position',[left bottom width height]) creates an
    %   axes at the specified position in normalized coordinates (in
    %   in the range from 0.0 to 1.0).
    %
    %   SUBPLOT(..., PROP1, VALUE1, PROP2, VALUE2, ...) sets the
    %   specified property-value pairs on the subplot axes. To add the
    %   subplot to a specific figure pass the figure handle as the
    %   value for the 'Parent' property.
    %
    %   If a SUBPLOT specification causes a new axes to overlap an
    %   existing axes, the existing axes is deleted - unless the position
    %   of the new and existing axes are identical.  For example,
    %   the statement SUBPLOT(1,2,1) deletes all existing axes overlapping
    %   the left side of the Figure window and creates a new axes on that
    %   side - unless there is an axes there with a position that exactly
    %   matches the position of the new axes (and 'replace' was not specified),
    %   in which case all other overlapping axes will be deleted and the
    %   matching axes will become the current axes.
    %
    %   SUBPLOT(111) is an exception to the rules above, and is not
    %   identical in behavior to SUBPLOT(1,1,1).  For reasons of backwards
    %   compatibility, it is a special case of subplot which does not
    %   immediately create an axes, but instead sets up the figure so that
    %   the next graphics command executes CLF RESET in the figure
    %   (deleting all children of the figure), and creates a new axes in
    %   the default position.  This syntax does not return a handle, so it
    %   is an error to specify a return argument.  The delayed CLF RESET
    %   is accomplished by setting the figure's NextPlot to 'replace'.
    %
    %   Be aware when creating subplots from scripts that the Position
    %   property of subplots is not finalized until either a drawnow
    %   command is issued, or MATLAB returns to await a user command.
    %   That is, the value obtained for subplot i by the command
    %   get(h(i),'Position') will not be correct until the script
    %   refreshes the plot or exits.
    %
    %   See also  GCA, GCF, AXES, FIGURE, UIPANEL
    
    %   SUBPLOT(m,n,p,H) when H is an axes will move H to the specified
    %   position.
    %   SUBPLOT(m,n,p,H,PROP1,VALUE1,...) will move H and apply the
    %   specified property-value pairs
    %
    %   SUBPLOT(m,n,p) for non-integer p places the subplot at the
    %   fraction p-floor(p) between the positions floor(p) and ceil(p)
    
    %   Copyright 1984-2010 The MathWorks, Inc.
    
    narg = nargin;
    
    % First we check whether Handle Graphics uses MATLAB classes
    isHGUsingMATLABClasses = feature('HGUsingMATLABClasses');
    
    if isHGUsingMATLABClasses
        if nargout == 0
            if narg == 0
                subplotHGUsingMATLABClasses();
            elseif narg == 1
                subplotHGUsingMATLABClasses(nRows);
            elseif narg == 2
                subplotHGUsingMATLABClasses(nRows, nCols);
            elseif narg == 3
                subplotHGUsingMATLABClasses(nRows, nCols, plotId);
            else
                subplotHGUsingMATLABClasses(nRows, nCols, plotId, varargin{:});
            end
        else
            if narg == 0
                theAxis = subplotHGUsingMATLABClasses();
            elseif narg == 1
                theAxis = subplotHGUsingMATLABClasses(nRows);
            elseif narg == 2
                theAxis = subplotHGUsingMATLABClasses(nRows, nCols);
            elseif narg == 3
                theAxis = subplotHGUsingMATLABClasses(nRows, nCols, plotId);
            else
                theAxis = subplotHGUsingMATLABClasses(nRows, nCols, plotId, varargin{:});
            end
        end
    else
        % we will kill all overlapping axes siblings if we encounter the mnp
        % or m,n,p specifier (excluding '111').
        % But if we get the 'position' or H specifier, we won't check for and
        % delete overlapping siblings:
        killSiblings = 0;
        createAxis = true;
        moveAxis = false;
        delayDestroy = false;
        useAutoLayout = true;
        tol = sqrt(eps);
        parent = get(0, 'CurrentFigure');
        ancestorFigure = parent;
        if ~isempty(parent) && ~isempty(get(parent, 'CurrentAxes'))
            parent = get(get(parent, 'CurrentAxes'), 'Parent');
            ancestorFigure = parent;
            if ~strcmp(get(ancestorFigure, 'Type'), 'figure')
                ancestorFigure = ancestor(parent, 'figure');
            end
        end
        pvpairs = {};
        preventMove = false;
        % This is the percent offset from the subplot grid of the plotbox.
        inset = [.2, .18, .04, .1]; % [left bottom right top]
        inset = [.02 .05 .02 .04]; % ashao: found online
        %check for encoded format
        h = [];
        position = [];
        explicitParent = false;
        explicitPosition = false;
        
        if narg == 0 % make compatible with 3.5, i.e. subplot == subplot(111)
            nRows = 111;
            narg = 1;
        end
        
        if narg == 1
            % The argument could be one of 3 things:
            % 1) a 3-digit number 100 < num < 1000, of the format mnp
            % 2) a 3-character string containing a number as above
            % 3) an axes handle
            arg = nRows;
            
            % turn string into a number:
            if(ischar(arg))
                arg = str2double(arg);
            end
            
            % Check for NaN and Inf.
            if (~isfinite(arg))
                error(id('SubplotIndexNonFinite'), 'Index must be a finite 3-digit number of the format mnp.')
            end
            
            % number with a fractional part can only be an identifier:
            if(rem(arg, 1) > 0)
                h = arg;
                if ~ishghandle(h, 'axes')
                    error(id('InvalidAxesHandle'), 'Requires valid axes handle for input.')
                end
                createAxis = false;
                % all other numbers will be converted to mnp format:
            else
                % Check for input out of range
                if (arg <= 100 || arg >= 1000)
                    error(id('SubplotIndexOutOfRange'), 'Index must be a 3-digit number of the format mnp.')
                end
                
                plotId = rem(arg, 10);
                nCols = rem(fix(arg - plotId) / 10, 10);
                nRows = fix(arg / 100);
                if nRows * nCols < plotId
                    error(id('SubplotIndexTooLarge'), 'Index exceeds number of subplots.');
                end
                killSiblings = 1;
                if (arg == 111)
                    createAxis = false;
                    delayDestroy = true;
                else
                    createAxis = true;
                    delayDestroy = false;
                end
            end
            
        elseif narg == 2
            % The arguments MUST be the string 'position' and a 4-element vector:
            if (strcmpi(nRows, 'position'))
                pos_size = size(nCols);
                if (pos_size(1) * pos_size(2) == 4)
                    position = nCols;
                    explicitPosition = true;
                else
                    error(id('InvalidPositionParameter'), ...
                        'Position must be of the form [left bottom width height].')
                end
            else
                error(id('UnknownOption'), 'Unknown command option.')
            end
            killSiblings = 1; % Kill overlaps here also.
            useAutoLayout = false;
            
        elseif narg == 3
            % passed in subplot(m,n,p) -- we should kill overlaps
            % here too:
            killSiblings = 1;
            
        elseif narg >= 4
            if ~ischar(nRows)
                arg = varargin{1};
                if ~ischar(arg)
                    % passed in subplot(m,n,p,H,...)
                    h = arg;
                    if ~ishghandle(h, 'axes') || ...
                            isa(handle(h), 'scribe.colorbar') || ...
                            isa(handle(h), 'scribe.legend')
                        error(id('InvalidAxesHandle'), 'Requires valid axes handle for input.')
                    end
                    parent = get(h, 'Parent');
                    ancestorFigure = ancestor(h, 'figure');
                    % If the parent is passed in explicitly, don't create a new figure
                    % when the "NextPlot" property is set to "new" in the figure.
                    explicitParent = true;
                    set(ancestorFigure, 'CurrentAxes', h);
                    moveAxis = true;
                    createAxis = false;
                    if narg >= 5 && strcmpi(varargin{2}, 'PreventMove')
                        preventMove = true;
                        pvpairs = varargin(3 : end);
                    else
                        pvpairs = varargin(2 : end);
                    end
                elseif strncmpi(arg, 'replace', 1)
                    % passed in subplot(m,n,p,'replace')
                    killSiblings = 2; % kill nomatter what
                elseif strcmpi(arg, 'align')
                    % passed in subplot(m,n,p,'align')
                    % since obeying position will remove the axes from the grid just set
                    % useAutoLayout to false to skip adding it to the grid to start with
                    useAutoLayout = false;
                    killSiblings = 1; % kill if it overlaps stuff
                elseif strcmpi(arg, 'v6')
                    % passed in subplot(m,n,p,'v6')
                    % since obeying position will remove the axes from the grid just set
                    % useAutoLayout to false to skip adding it to the grid to start with
                    warning(['MATLAB:', mfilename, ':DeprecatedV6Argument'],...
                        ['The ''v6'' argument to %s is deprecated,',...
                        ' and will no longer be supported in a future release.'], upper(mfilename));
                    useAutoLayout = false;
                    killSiblings = 1; % kill if it overlaps stuff
                else
                    % passed in prop-value pairs
                    killSiblings = 1;
                    pvpairs = varargin;
                    par = find(strncmpi('Parent', pvpairs(1 : 2 : end), 6));
                    if any(par)
                        % If the parent is passed in explicitly, don't create a new figure
                        % when the "NextPlot" property is set to "new" in the figure.
                        explicitParent = true;
                        parent = varargin{2 * par(1)};
                        ancestorFigure = ancestor(parent, 'figure');
                    end
                end
            else
                % Passed in "Position" syntax with P/V pairs
                % The arguments MUST be the string 'position' and a 4-element vector:
                if (strcmpi(nRows, 'position'))
                    pos_size = size(nCols);
                    if (pos_size(1) * pos_size(2) == 4)
                        position = nCols;
                        explicitPosition = true;
                    else
                        error(id('InvalidPositionParameter'), ...
                            'Position must be of the form [left bottom width height].')
                    end
                else
                    error(id('UnknownOption'), 'Unknown command option.')
                end
                killSiblings = 1; % Kill overlaps here also.
                useAutoLayout = false;
                pvpairs = [{plotId}, varargin];
                par = find(strncmpi('Parent', pvpairs(1 : 2 : end), 6));
                if any(par)
                    % If the parent is passed in explicitly, don't create a new figure
                    % when the "NextPlot" property is set to "new" in the figure.
                    explicitParent = true;
                    parent = pvpairs{2 * par(1)};
                    ancestorFigure = ancestor(parent, 'figure');
                end
            end
        end
        
        % if we recovered an identifier earlier, use it:
        if ~isempty(h) && ~moveAxis
            parent = get(h, 'Parent');
            ancestorFigure = ancestor(h, 'figure');
            set(ancestorFigure, 'CurrentAxes', h);
        else  % if we haven't recovered position yet, generate it from mnp info:
            if isempty(parent)
                parent = gcf;
                ancestorFigure = parent;
            end
            if isempty(position)
                if min(plotId) < 1
                    error(id('SubplotIndexTooSmall'), 'Illegal plot number.')
                elseif max(plotId) > nCols * nRows
                    error(id('SubplotIndexTooLarge'), 'Index exceeds number of subplots.');
                else
                    
                    row = (nRows - 1) - fix((plotId - 1) / nCols);
                    col = rem(plotId - 1, nCols);
                    
                    % get default axes position in normalized units
                    % If we have checked this quanitity once, cache it.
                    if ~isappdata(ancestorFigure, 'SubplotDefaultAxesLocation')
                        if ~strcmp(get(ancestorFigure, 'DefaultAxesUnits'), 'normalized')
                            tmp = axes;
                            set(tmp, 'Units', 'normalized')
                            def_pos = get(tmp, 'Position');
                            delete(tmp)
                        else
                            def_pos = get(ancestorFigure, 'DefaultAxesPosition');
                        end
                        setappdata(ancestorFigure, 'SubplotDefaultAxesLocation', def_pos);
                    else
                        def_pos = getappdata(ancestorFigure, 'SubplotDefaultAxesLocation');
                    end
                    
                    % compute outerposition and insets relative to figure bounds
                    rw = max(row) - min(row) + 1;
                    cw = max(col) - min(col) + 1;
                    width = def_pos(3) / (nCols - inset(1) - inset(3));
                    height = def_pos(4) / (nRows - inset(2) - inset(4));
                    inset = inset .* [width, height, width, height];
                    outerpos = [def_pos(1) + min(col) * width - inset(1), ...
                        def_pos(2) + min(row) * height - inset(2), ...
                        width * cw, height * rw];
                    
                    % compute inner position
                    position = [outerpos(1 : 2) + inset(1 : 2), ...
                        outerpos(3 : 4) - inset(1 : 2) - inset(3 : 4)];
                    
                end
            end
        end
        
        % kill overlapping siblings if mnp specifier was used:
        nextstate = get(ancestorFigure, 'NextPlot');
        
        if strncmp(nextstate, 'replace', 7)
            nextstate = 'add';
        elseif strncmp(nextstate, 'new', 3)
            killSiblings = 0;
        end
        
        if killSiblings
            if delayDestroy
                if nargout
                    error(id('TooManyOutputs'), ...
                        'Function called with too many output arguments.')
                else
                    set(ancestorFigure, 'NextPlot', 'replace');
                    return
                end
            end
            sibs = datachildren(parent);
            newcurrent = [];
            for i = 1 : length(sibs)
                % Be aware that handles in this list might be destroyed before
                % we get to them, because of other objects' DeleteFcn callbacks...
                if ishghandle(sibs(i), 'axes')
                    units = get(sibs(i), 'Units');
                    sibpos = get(sibs(i), 'Position');
                    % If a legend or colorbar has resized the axes, use the original axes
                    % position as the "Position" property:
                    if ~explicitPosition
                        if isappdata(sibs(i), 'LegendColorbarExpectedPosition') && ...
                                isequal(getappdata(sibs(i), 'LegendColorbarExpectedPosition'), get(sibs(i), 'Position'))
                            sibinset = getappdata(sibs(i), 'LegendColorbarOriginalInset');
                            if isempty(sibinset)
                                % during load the appdata might not be present
                                sibinset = get(get(sibs(i), 'Parent'), 'DefaultAxesLooseInset');
                            end
                            sibinset = offsetsInUnits(sibs(i), sibinset, 'normalized', get(sibs(i), 'Units'));
                            if strcmpi(get(sibs(i), 'ActivePositionProperty'), 'position')
                                pos = get(sibs(i), 'Position');
                                loose = get(sibs(i), 'LooseInset');
                                opos = getOuterFromPosAndLoose(pos, loose, get(sibs(i), 'Units'));
                                if strcmp(get(sibs(i), 'Units'), 'normalized')
                                    sibinset = [opos(3 : 4), opos(3 : 4)] .* sibinset;
                                end
                                sibpos = [opos(1 : 2) + sibinset(1 : 2), opos(3 : 4) - sibinset(1 : 2) - sibinset(3 : 4)];
                            end
                        end
                    end
                    if ~strcmp(units, 'normalized')
                        sibpos = hgconvertunits(ancestorFigure, sibpos, units, 'normalized', parent);
                    end
                    intersect = 1;
                    if ((position(1) >= sibpos(1) + sibpos(3) - tol) || ...
                            (sibpos(1) >= position(1) + position(3) - tol) || ...
                            (position(2) >= sibpos(2) + sibpos(4) - tol) || ...
                            (sibpos(2) >= position(2) + position(4) - tol))
                        intersect = 0;
                    end
                    if intersect
                        % position is the proposed position of an axes, and
                        % sibpos is the current position of an existing axes.
                        % Since the bounding boxes of position and sibpos overlap,
                        % we must determine whether to delete the sibling sibs(i)
                        % whose normalized position is sibpos.
                        
                        % First of all, we check whether we must kill the sibling
                        % "no matter what."
                        if (killSiblings == 2)
                            delete(sibs(i));
                            
                            % If the proposed and existing axes overlap exactly, we do
                            % not kill the sibling.  Rather we shall ensure later that
                            % this sibling axes is set as the 'CurrentAxes' of its
                            % ancestorFigure.
                            
                            % Next we check for a partial overlap.
                        elseif (any(abs(sibpos - position) > tol))
                            % The proposed and existing axes partially overlap.
                            % Since the proposed and existing axes could each be
                            % "grid-generated" or "explicitly-specified", we must
                            % consider four possibilities for the overlap of
                            % "proposed" vs. "existing", i.e.
                            % (1) "grid-generated" vs. "grid-generated"
                            % (2) "grid-generated" vs. "explicitly-specified"
                            % (3) "explicitly-specified" vs. "grid-generated"
                            % (4) "explicitly-specified" vs. "explicitly-specified"
                            
                            % If the position of the proposed axes is
                            % "explicitly-specified", then the only condition that
                            % avoids killing the sibling is an exact overlap.
                            % However, we know that the overlap is partial.
                            if (explicitPosition)
                                delete(sibs(i));
                            else
                                % We know that the position of the proposed axes is
                                % "grid-generated".
                                
                                grid = getappdata(parent, 'SubplotGrid');
                                % The SubplotGrid maintains an array of axes
                                % handles, one per grid location.  Axes that span
                                % multiple grid locations do not store handles in
                                % the SubplotGrid.
                                
                                if isempty(grid) || ~any(grid(:) == sibs(i)) || ...
                                        size(grid, 1) ~= nRows || size(grid, 2) ~= nCols || ...
                                        ~isscalar(row) || ~isscalar(col)
                                    % If the sibling cannot be found in the grid, we
                                    % kill the sibling.  Otherwise, the proposed and
                                    % existing axes are "grid-generated".  If we
                                    % are changing the size of the grid, we kill
                                    % the sibling.  Otherwise, "plotId" may be a
                                    % vector of multiple grid locations, which
                                    % causes a partial overlap between the proposed
                                    % and existing axes, so we kill the sibling.
                                    
                                    % This check recognizes that there may be
                                    % labels, colorbars, legends, etc. attached to
                                    % the existing axes that have affected its
                                    % position.  In such a case, we do not kill the
                                    % sibling.
                                    delete(sibs(i));
                                end
                            end
                        end
                        if ishghandle(sibs(i))
                            % if this axes overlaps the other one exactly then
                            if ~isempty(newcurrent) && ishghandle(newcurrent)
                                delete(newcurrent);
                            end
                            newcurrent = sibs(i);
                        end
                    end
                end
            end
            if ~isempty(newcurrent) && ishghandle(newcurrent)
                set(ancestorFigure, 'CurrentAxes', newcurrent);
                createAxis = false;
            end
            set(ancestorFigure, 'NextPlot', nextstate);
        end
        
        % create the axes:
        if createAxis
            if strcmp(nextstate, 'new') && ~explicitParent
                parent = figure;
                ancestorFigure = parent;
            end
            ax = axes('Units', 'normalized', 'Position', position, ...
                'LooseInset', inset, 'Parent', parent);
            % TODO: Get axes to accept position args on command line
            set(ax, 'Units', get(double(ancestorFigure), 'DefaultAxesUnits'))
            if useAutoLayout
                addAxesToGrid(ax, nRows, nCols, row, col, position);
            end
            if ~isempty(pvpairs)
                set(ax, pvpairs{:});
            end
        elseif moveAxis && ~preventMove
            ax = h;
            units = get(h, 'Units');
            set(h, 'Units', 'normalized', 'Position', position, ...
                'LooseInset', inset, 'Parent', parent);
            set(h, 'Units', units);
            if useAutoLayout
                addAxesToGrid(ax, nRows, nCols, row, col, position);
            end
            if ~isempty(pvpairs)
                set(h, pvpairs{:});
            end
        else
            % this should only happen with subplot(H)
            ax = get(ancestorFigure, 'CurrentAxes');
        end
        
        % return identifier, if requested:
        if(nargout > 0)
            theAxis = ax;
        end
    end
end


% Create subplot listeners to align plot boxes automatically
function createListeners(p, axlisth)
    setappdata(p, 'SubplotListeners', [])
    fig = p;
    if ~isequal(get(fig, 'Type'), 'figure')
        fig = ancestor(fig, 'figure');
    end
    
    list = [ ...
        handle.listener(axlisth, findprop(axlisth(1), 'Units'), ... % must be first
        'PropertyPostSet', @axesUnitsPostSet);
        handle.listener(axlisth, findprop(axlisth(1), 'Units'), ...
        'PropertyPreSet', @axesUnitsPreSet);
        handle.listener(axlisth, findprop(axlisth(1), 'Position'), ...
        'PropertyPostSet', @axesMoved);
        handle.listener(axlisth, findprop(axlisth(1), 'ActivePositionProperty'), ...
        'PropertyPreSet', @axesMoved);
        handle.listener(axlisth, findprop(axlisth(1), 'Parent'), ...
        'PropertyPreSet', @axesMoved);
        handle.listener(axlisth, 'AxisInvalidEvent', ...
        {@subplotlayoutInvalid, p});
        handle.listener(handle(fig), 'FigureUpdateEvent', ...
        {@subplotlayout, p})];
    
    for k = 1 : length(axlisth)
        ax = axlisth(k);
        if ~isappdata(double(ax), 'SubplotDeleteListener')
            setappdata(double(ax), 'SubplotDeleteListener', ...
                handle.listener(ax, 'ObjectBeingDestroyed', ...
                @axesDestroyed));
        end
    end
    setappdata(p, 'SubplotListeners', list)
end


% Add ax to a matrix of handles in the specified location.
% The grid is stored on the parent appdata.
% Also store the insets in ax appdata.
% Only stores the axes if it is in a 1-by-1 cell and
% the grid size matches any existing grid.
function addAxesToGrid(ax, nRows, nCols, row, col, position)
    p = get(ax, 'Parent');
    grid = getappdata(p, 'SubplotGrid');
    if isempty(grid)
        grid = zeros(nRows, nCols);
    end
    if any(size(grid) ~= [nRows, nCols])
        return
    end
    if length(row) ~= 1 || length(col) ~= 1
        return
    end
    if round(row) ~= row || round(col) ~= col
        return
    end
    if grid(row + 1, col + 1) == ax
        return
    end
    grid(row + 1, col + 1) = ax;
    list = grid(:);
    list(list == 0) = []; % remove root
    list(~ishghandle(list)) = []; % remove invalid handles
    createListeners(p, handle(list));
    setappdata(p, 'SubplotGrid', grid)
    setappdata(ax, 'SubplotPosition', position); % normalized
    subplotlayoutInvalid(handle(ax), [], p);
end

% Remove ax from grid of subplots in p
function removeAxesFromGrid(p, ax)
    grid = getappdata(p, 'SubplotGrid');
    if ~isempty(grid)
        n = grid == ax;
        if any(n(:))
            grid(n) = 0;
            list = grid(:);
            list(list == 0) = []; % remove root
            list(~ishghandle(list)) = [];
            if isempty(list)
                rmappdata(p, 'SubplotListeners');
                rmappdata(p, 'SubplotGrid');
            else
                setappdata(p, 'SubplotGrid', grid);
            end
        end
    end
end

% Callback when axes moves to remove it from subplot layout grid
function axesMoved(~, evdata)
    ax = double(evdata.AffectedObject);
    % If the legend or colorbar is causing the move, do not remove the axes
    % from the subplot grid. Do, however, update it's cached position:
    if (isappdata(ax, 'inLayout') && ~isempty(getappdata(ax, 'inLayout'))) || ...
            isappdata(ax, 'LegendColorbarReclaimSpace')
        setappdata(ax, 'SubplotPosition', get(ax, 'Position'));
    else
        removeAxesFromGrid(get(ax, 'Parent'), ax);
    end
end

% Callback when axes changes units
function axesUnitsPreSet(hSrc,evdata) %#ok
    ax = double(evdata.AffectedObject);
    p = get(ax,'Parent');
    list = getappdata(p,'SubplotListeners');
    if ~isempty(list)
        if isa(list,'handle.listener')
            set(list(2:end),'enable','off');
        else
            for i=2:length(list)
                list(i).Enabled = false;
            end
        end
    end
end

% Callback when axes done changing units
function axesUnitsPostSet(hSrc,evdata) %#ok
    ax = double(evdata.AffectedObject);
    p = get(ax,'Parent');
    list = getappdata(p,'SubplotListeners');
    if ~isempty(list)
        if isa(list,'handle.listener')
            set(list(2:end),'enable','on');      
        else
            for i=2:length(list)
                list(i).Enabled = true;
            end
        end
    end
end

% Callback when axes destroyed
function axesDestroyed(hSrc, ~)
    ax = double(hSrc);
    p = get(ax, 'Parent');
    
    if strcmp(get(p, 'BeingDeleted'), 'off')
        removeAxesFromGrid(p, ax);
    elseif isappdata(p, 'SubplotListeners')
        rmappdata(p, 'SubplotListeners');
        rmappdata(p, 'SubplotGrid');
    end
end

function str = id(str)
    str = ['MATLAB:subplot:' str];
end

%----------------------------------------------------------------%
% Convert units of offsets like LooseInset or TightInset
% Note: Copied from legendcolorbarlayout.m
function out = offsetsInUnits(ax, in, from, to)
    fig = ancestor(ax, 'figure');
    par = get(ax, 'Parent');
    p1 = hgconvertunits(fig, [0, 0, in(1 : 2)], from, to, par);
    p2 = hgconvertunits(fig, [0, 0, in(3 : 4)], from, to, par);
    out = [p1(3 : 4), p2(3 : 4)];
end

%----------------------------------------------------------------%
% Compute reference OuterPos from pos and loose. Note that
% loose insets are relative to outerposition
% Note: Copied from legendcolorbarlayout.m
function outer = getOuterFromPosAndLoose(pos, loose, units)
    if strcmp(units, 'normalized')
        % compute outer width and height and normalize loose to them
        w = pos(3) / (1 - loose(1) - loose(3));
        h = pos(4) / (1 - loose(2) - loose(4));
        loose = [w, h, w, h] .* loose;
    end
    outer = [pos(1 : 2) - loose(1 : 2), pos(3 : 4) + loose(1 : 2) + loose(3 : 4)];
end
