function [ numempty ] = timevary_find_mean_centers( inpath, outpath )
numempty = 0;
files = dir([inpath 't*.mat']);
nfiles = length(files)

for tidx = 1:nfiles
    clf; hold on;
    fprintf('Track %03d\n',tidx);
    load([inpath filesep files(tidx).name],'track');
    track.tracknum = tidx;   
    opt_track = find_mean_centers( track );
    if isempty(opt_track)
        numempty = numempty+1
    end   
    save([outpath filesep files(tidx).name],'opt_track')
end

end

function  opt_track  = find_mean_centers( intrack )
%
% 
offset_time = datenum(1950,1,1); % Conversion from AVISO to matlab date format
% first_time = offset_time+15608; % Mintime based on first T/P measurement
% first_time = datenum(1993,1,1)-1;
% last_time = offset_time + 23004.96;
% window_length = 365; % Length of window is 1000 days ~~ 100 measurements
% slide_length = 365;
% 
% start_window = first_time;
% end_window = start_window + window_length;
window_counter = 0;
for year = 1993:2012

% window_counter = 0; 
% while end_window <  last_time

    window_counter = window_counter + 1;
    avgtracktime = nanmean(intrack.time,2);    
    winidx = (avgtracktime + offset_time) > datenum(year,1,1)-1 & ...
        (avgtracktime + offset_time) < datenum(year+1,1,1);
    
    [npts nfiles] = size(intrack.sla);
    latrange = intrack.latrange;
    buffer = 5; % How far out to extend the previously selected latrage
    
    minlat = min(latrange)-5;
    maxlat = max(latrange)+5;
    
    trunc.idx = intrack.lat >= minlat & intrack.lat <= maxlat;
    
    intrack.skewness = skewness(intrack.sla(winidx,:),1);
    siglevel = nanmean(abs(intrack.skewness));
    if siglevel > 0.5
        siglevel = 0.5
    end
    
    intrack.std = nanstd(intrack.sla(winidx,:),1);
    intrack.mean = nanmean(intrack.sla(winidx,:),1);
    intrack.smoothskew=smooth(intrack.skewness,5,'mean');
    
    trunc.lat = intrack.lat(trunc.idx);
    trunc.mean = intrack.mean(trunc.idx);
    trunc.skewness = intrack.skewness(trunc.idx);
    trunc.smoothskew =intrack.smoothskew(trunc.idx);
    trunc.std = intrack.std(trunc.idx);
%     plot(trunc.lat,trunc.smoothskew); hold on;
    
    [maxtab mintab] = peakdet(trunc.smoothskew,1.5*siglevel);
    
    if ~isempty(maxtab) & ~isempty(mintab)
        maxtab(:,3) = trunc.lat(maxtab(:,1));
        mintab(:,3) = trunc.lat(mintab(:,1));
        
        
        maxtab = sortrows(maxtab,3);
        mintab = sortrows(mintab,3);
        %
        %     inrange = maxtab(:,3) >= minlat & maxtab(:,3) <= maxlat;
        %     maxtab = maxtab(inrange,:);
        %
        %     inrange = mintab(:,3) >= minlat & mintab(:,3) <= maxlat;
        %     mintab = mintab(inrange,:);
        

        while maxtab(1,3) > mintab(1,3)
            mintab(1,:) = [];
            if isempty(mintab)
                break
            end
        end
        ntrans = min([size(maxtab,1) size(mintab,1)]);
        
        if ntrans > 0
            
            maxtab=maxtab(1:ntrans,:);
            mintab=mintab(1:ntrans,:);
            
            for i = 1:ntrans
                
                transidx = min([maxtab(i,1) mintab(i,1)]):max([maxtab(i,1) mintab(i,1)]);
                opt_track(window_counter).tracknum = intrack.tracknum;
                opt_track(window_counter).skewness{i} = trunc.skewness(transidx);
                opt_track(window_counter).mean{i} = trunc.mean(transidx);
                opt_track(window_counter).std{i} = nanmean(trunc.std(transidx));
                opt_track(window_counter).lat{i} = trunc.lat(transidx);
                opt_track(window_counter).time = datenum(year,1,1) + 365/2;
                
                %             opt_track(counter).lon{i}= trunc.lon(transidx);
%                 plot(opt_track(window_counter).lat{i},opt_track(window_counter).skewness{i},'LineWidth',4)
            end
            
        else
            opt_track = [];
            
        end
        fprintf('Window: %d Transitions %d\n',window_counter,ntrans)
        
    end
    
%     start_window = start_window + slide_length;
%     end_window = start_window + window_length;
end
end

