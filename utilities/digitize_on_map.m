function [ outlat, outlon ] = digitize_on_map( h, varargin )
% Parse output from inputm to create vectors of the digitized points
% Mouse button 1 enters a point
% Mouse button 2 deletes the previous point
% Mouse button 3 ends the digitizing
% Inputs:
%   h: figure handle (optional)
%   varargin: Options to plot the digitized poitns
if nargin < 1
    h = gcf;
end
button = 1;

% Preallocate a large array, not really necessary
outlat = nan(5000,1);
outlon = nan(5000,1);

counter = 0;

figure(h),hold on; % Set active figure to acquire from
if ishold
    rehold = true;
end

while button ~=3
    
    
    [inlat inlon button] = inputm;
    switch button
        
        case 1, % Store a point
            counter = counter+1;
            outlat(counter) = inlat;
            outlon(counter) = inlon;
            if counter > 1
                h = plotm([outlat(counter-1) outlat(counter)],[outlon(counter-1) outlon],varargin{:})
            end
        case 2,
            counter = counter-1;
            outlat(counter) = [];
            outlon(counter) = [];
            delete(h)
        case 3,
            fprintf('Digitized %d points\n',counter);
            break;
        otherwise,
            disp('Invalid button, no point stored')
            
            
    end
    
end

% Delete any extra points
outlat(counter:end) = [];
outlon(counter:end) = [];

if ~rehold
    hold off;
end