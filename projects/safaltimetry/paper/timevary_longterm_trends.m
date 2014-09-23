inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/timevary/annual/processed/';
files = dir([inpath 't*.mat']);
nfiles = length(files);
trackfile = '/ltraid4/aviso/alongtrack/sla/vxxc_matlab/groundpath.mat';
outpath = '/ltraid3/ashao/uw-apl/matlab/projects/safaltimetry/paper/figs/';
load(trackfile)
nwindows = 20;
%% Extract by track
extents = {'north','mean','south'};
nextents = length(extents);
R2threshold = 0.4;
years =1993:2012;
% Initialize data arrays
for extidx = 1:nextents
    
    extent = extents{extidx};
    alongtrack.(extent).lat = zeros(nfiles,nwindows);
    alongtrack.(extent).lon = zeros(nfiles,nwindows);
    alongtrack.(extent).R2 = zeros(nfiles,nwindows);
    
end

% Load the output and store it into the alongtrack structure
for tidx = 1:nfiles
    fprintf('Working on Track %d\n',tidx)
    % Load the track
    load([inpath files(tidx).name]);
    
    for extidx = 1:nextents
        extent = extents{extidx};
        % Loop over the windows
        for winidx = 1:nwindows
            % Extract the data for the current time window
            optpar = cell2mat(opt_track(winidx).optpar');
            
            % Apply two QC criteria
            minlat = min(groundpath(tidx).lat);
            R2 = cell2mat(opt_track(winidx).R2);
            delidx = optpar(:,2) < minlat | R2' < R2threshold';
            optpar(delidx,:) = [];
            
            % Check to see that we do actually have valid points
            if ~isempty(optpar)
                
                switch extent
                    case 'north'
                        alongtrack.(extent).lat(tidx,winidx) = ...
                            max(optpar(:,2));
                    case 'mean'
                        alongtrack.(extent).lat(tidx,winidx) = ...
                            mean(optpar(:,2));
                    case 'south'
                        alongtrack.(extent).lat(tidx,winidx) = ...
                            min(optpar(:,2));
                end
                
                % Interpolate to find the alongtrack longitude
                alongtrack.(extent).lon(tidx,winidx) = ...
                    intrplon(groundpath(tidx).lat,groundpath(tidx).lon,...
                    alongtrack.(extent).lat(tidx,winidx));
            else % Fill in values with NaN if optpar is empty
                
                alongtrack.(extent).lon(tidx,winidx) = NaN;
                alongtrack.(extent).lat(tidx,winidx) = NaN;
            end
        end
        % Interpolate to fill in the time gaps
        intrpidx = isnan(alongtrack.(extent).lat(tidx,:));
        if any(intrpidx)
            alongtrack.(extent).lat(tidx,intrpidx) = ...
                interp1(years(~intrpidx), ...
                alongtrack.(extent).lat(tidx,~intrpidx), ...
                years(intrpidx));
            alongtrack.(extent).lon(tidx,intrpidx) = ...
                intrplon(groundpath(tidx).lat,groundpath(tidx).lon,...
                alongtrack.(extent).lat(tidx,intrpidx));
        end
        
        
        % If the time gap was at the end, fill in the NaN values with the
        % mean
        fillidx = isnan(alongtrack.(extent).lat(tidx,:));
        if any(fillidx)
            alongtrack.(extent).lat(tidx,fillidx) = ...
                nanmean(alongtrack.(extent).lat(tidx,:));
            alongtrack.(extent).lon(tidx,fillidx) = ...
                intrplon(groundpath(tidx).lat,groundpath(tidx).lon,...
                alongtrack.(extent).lat(tidx,fillidx));
        end
    end
end

%% Estimate the linear trends using regression for each extent and track
for extidx = 1:nextents
    extent = extents{extidx};
    for tidx = 1:nfiles
        
        mdl = LinearModel.fit(0:10,alongtrack.(extent).lat(tidx,1:11));
        alongtrack.(extent).trend(tidx) = mdl.Coefficients.Estimate(2);
        alongtrack.(extent).trend_pvalue(tidx) = mdl.Coefficients.pValue(2);
        
    end
    
end
%% Plot the results of the regression

pval_limit = 0.05;
colors = {'r','k','b'};
clf
for extidx = [1 3]
    extent = extents{extidx};
    sigidx = alongtrack.(extent).trend_pvalue < pval_limit;
    
    [plotlat plotlon] = meanm(alongtrack.(extent).lat', ...
        alongtrack.(extent).lon');    
    
    plotlon = wrapTo360(plotlon);
    
    subplot(2,1,1); hold on;
    stem(plotlon(sigidx),alongtrack.(extent).trend(sigidx),colors{extidx},'filled')
    xlim([0 360])
    set(gca,'XTick',0:30:360)
    grid on
    box on;    
    ylabel('Meridional Trend (\circS yr^{-1})')
    xlabel('Longitude (\circE)')
    subplot(2,1,2); hold on;
    [sortlon sortidx] = sort(plotlon);
    plot(sortlon,smooth(std(alongtrack.(extent).lat(sortidx,:),0,2),10,'mean') ...
        ,colors{extidx},'LineWidth',2)
    xlim([0 360])
    set(gca,'XTick',0:30:360)
    grid on
    box on;
    xlabel('Longitude (\circE)')
    ylabel('Standard Deviation (\circ)') 
    
end
%% Save results
saveas(gcf,[outpath 'figure_longterm_trends.eps'],'epsc')

