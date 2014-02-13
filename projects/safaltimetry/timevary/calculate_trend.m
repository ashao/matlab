datapath = 'C:\Users\ashao\Data\2013.16.12\';
trackfiles = dir([datapath '*.mat']);
nwindows = 37;
ntracks = length(trackfiles);



saf.north.trend = zeros(ntracks,1);
saf.center.trend = zeros(ntracks,1);
saf.south.trend = zeros(ntracks,1);
saf.tracknum = zeros(ntracks,1);

fillmat = zeros(nwindows,ntracks);
saf.time = fillmat;
saf.north.lat = fillmat;
saf.north.lon = fillmat;
saf.center.lat = fillmat;
saf.center.lon = fillmat;
saf.south.lat = fillmat;
saf.south.lon = fillmat;



for track = 1:ntracks
   
    load([datapath trackfiles(track).name]);
    fprintf('Track %d ... ',track);
    for window = 1:nwindows
        
        output = cell2mat(opt_track(window).optimal');
        saf.north.lat(window,track) = max(output(:,2));
        saf.south.lat(window,track) = min(output(:,2));
        saf.center.lat(window,track) = mean(output(:,2));                                
        saf.time(window,track) = opt_track(window).time;
    end
    
    coeffs = polyfit(saf.time(:,track),saf.north.lat(:,track),1);
    saf.north.trend(track) = coeffs(1);
    coeffs = polyfit(saf.time(:,track),saf.center.lat(:,track),1);
    saf.center.trend(track) = coeffs(1);
    coeffs = polyfit(saf.time(:,track),saf.south.lat(:,track),1);
    saf.south.trend(track) = coeffs(1);
    saf.tracknum(track) = str2num(trackfiles(track).name(2:4));
    
    
    fprintf('Done!\n');
    
end

% %% Plot track by track
% 
% for track = 1:ntracks
%    
%     clf; hold on
%     plot(saf.north.lat(:,track))
%     plot(saf.center.lat(:,track))
%     plot(saf.south.lat(:,track))
%     
%     
% end