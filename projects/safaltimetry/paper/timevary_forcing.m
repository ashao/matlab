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
    alongtrack.(extent).meanlat = zeros(nfiles,1);
    alongtrack.(extent).meanlon = zeros(nfiles,1);
    
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
        
        [alongtrack.(extent).meanlat(tidx) ...
            alongtrack.(extent).meanlon(tidx)] = meanm( ...
            alongtrack.(extent).lat(tidx,:),alongtrack.(extent).lon(tidx,:));
        alongtrack.(extent).meanlon(tidx) = ...
            wrapTo360(alongtrack.(extent).meanlon(tidx));
    end
end
%% Load files and setup for the analysis
load /home/ashao/uw-apl/projects/saf_altimetry/ecmwf.tauanom.mat
load /ltraid3/ashao/uw-apl/projects/saf_altimetry/ecmwf.sam.mat
load /ltraid3/ashao/uw-apl/projects/saf_altimetry/nino34.annual.mat
year = datevec(sam.annual.time);
year = year(:,1);
yearidx = year >= 1993 & year <= 2012;
samidx = sam.annual.idx(yearidx);

yearidx = nino34.annual.year >= 1993 & nino34.annual.year <= 2012;
ensoidx = nino34.annual.idx(yearidx);


%%
lonranges = [0 60;60 120;120 180;180 240;240 300;300 360];
nranges = length(lonranges);
latrange = [-60 -40];
idxnames = {'sam','enso','wsi','wsc'};
for rangeidx = 1:nranges
    lonrange = lonranges(rangeidx,:);
    wsi = eof_annual(ecmwf.iews_anomaly,ecmwf.time,ecmwf.lat,ecmwf.lon, ...
        latrange,lonrange);    
    wsc = eof_annual(ecmwf.curl_anomaly,ecmwf.time,ecmwf.lat,ecmwf.lon, ...
        latrange,lonrange);
    fprintf('Longitude Range: %f %f\n',min(lonrange),max(lonrange));
    for extidx = 1:nextents
        extent = extents{extidx};
        lonidx = alongtrack.(extent).meanlon >= min(lonrange) & ...
            alongtrack.(extent).meanlon < max(lonrange);
        lattab = alongtrack.(extent).lat(lonidx,:);
        lattab = detrend(lattab);
        [U S V] = svd(lattab);
        pcs = U*S;
        fprintf('Extent: %s\n',extent)
        for nmodes = 1:20
           
            if nmodes == 1
                [C.sam lags.sam] = xcov(pcs(1,:),samidx,2,'coeff');
                [C.enso lags.enso] = xcov(pcs(1,:),ensoidx,2,'coeff');
                [C.wsi lags.wsi] = xcov(pcs(1,:),wsi,2,'coeff');
                [C.wsc lags.wsc] = xcov(pcs(1,:),wsc,2,'coeff');
            else
                [C.sam lags.sam] = xcov(sum(pcs(1:nmodes,:)),samidx,2,'coeff');
                [C.enso lags.enso] = xcov(sum(pcs(1:nmodes,:)),ensoidx,2,'coeff');
                [C.wsi lags.wsi] = xcov(sum(pcs(1:nmodes,:)),wsi,2,'coeff');
                [C.wsc lags.wsc] = xcov(sum(pcs(1:nmodes,:)),wsc,2,'coeff');
            end            
                        
            [null maxidx] = max(abs(C.sam));
            alongtrack.(extent).samcorr(nmodes) = C.sam(maxidx);
            alongtrack.(extent).samlag(nmodes) = lags.sam(maxidx);
            
            [null maxidx] = max(abs(C.enso));
            alongtrack.(extent).ensocorr(nmodes) = C.enso(maxidx);
            alongtrack.(extent).ensolag(nmodes) = lags.enso(maxidx);
            
            [null maxidx] = max(abs(C.wsi));
            alongtrack.(extent).wsicorr(nmodes) = C.wsi(maxidx);
            alongtrack.(extent).wsilag(nmodes) = lags.wsi(maxidx);
            
            [null maxidx] = max(abs(C.wsc));
            alongtrack.(extent).wsccorr(nmodes) = C.wsc(maxidx);
            alongtrack.(extent).wsclag(nmodes) = lags.wsc(maxidx);
            
            
        end
        
        [null maxidx] = max(C.sam);        
        fprintf('SAM Modes: %d Correlation %f Lag %d\n',maxidx, ...
            C.sam(maxidx),lags.sam(maxidx));
        
        [null maxidx] = max(C.enso);        
        fprintf('ENSO Modes: %d Correlation %f Lag %d\n',maxidx, ...
            C.enso(maxidx),lags.enso(maxidx));
        
        [null maxidx] = max(C.wsi);        
        fprintf('WSI Modes: %d Correlation %f Lag %d\n',maxidx, ...
            C.wsi(maxidx),lags.wsi(maxidx));
        
        [null maxidx] = max(C.wsc);        
        fprintf('WSC Modes: %d Correlation %f Lag %d\n',maxidx, ...
            C.wsc(maxidx),lags.wsc(maxidx));
        
    end
    
end

% wsi = eof_annual(ecmwf.iews_anomaly,ecmwf.time,ecmwf.lat,ecmwf.lon,latrange,lonranges(1));



