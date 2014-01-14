function [  ] = manual_lag_endpoints(lagpath,outpath,fronts,start_track)

tracks=trackconsec;
% tracks=1:254;
if nargin<4
    start_track=1;
end

if ~exist(outpath,'dir')
    mkdir(outpath)
end

nfronts=length(fronts);



% Plot possible frontal positions
subplot(2,2,[3 4])
%Set mapping projection

% SAF
m_proj('Mercator','lon',[0 360],'lat',[-65 -40]);
m_grid('ytick',[-65:5:-40]);
% Agulhas
% m_proj('Mercator','lon',[0 140],'lat',[-50 -10]);
% m_grid('ytick',[-50:5:-10]);

m_coast('patch',[0 0 0]);
% m_grid;
hold on
trackh=[];

if nargin <5
    centers_init=[];
end
for front=1:nfronts
    m_plot(fronts(front).lon,fronts(front).lat)
end
trackh=[];

for trackn=start_track:length(tracks)
    
    tracknum=tracks(trackn);
    load([lagpath filesep sprintf('t%03d.mat',tracknum)])
    lon_inter=[];
    lat_inter=[];
    notnan=find(~isnan(array.lon));
    for front=1:nfronts
        [templon templat]=intersections(fronts(front).lon,fronts(front).lat,array.lon(notnan),array.lat(notnan));
        lon_inter=[lon_inter ; templon];
        lat_inter=[lat_inter ; templat];
    end
    if isempty(lat_inter)
        lat_inter=-45;
    end
    
    subplot(2,2,[3 4])
    hold on
    title(sprintf('Track %d (%d/%d)',tracknum,trackn,length(tracks)))
    delete(trackh)
    trackh=m_plot(array.lon,array.lat);
    
    %% Positive Indices
    subplot(2,2,1,'replace')
    hold on
    notnan=find(~isnan(array.skewness_pos));
    filtskew=filtfilt(ones(3,1),1,array.skewness_pos(notnan))/3^2;
    plot(array.lat(notnan),filtskew)
    plot(array.lat,array.skewness_pos)
%     lat_inter
    xlim([min(lat_inter)-5 max(lat_inter)+5])
    scatter(lat_inter,zeros(size(lat_inter)),'ko','LineWidth',2)
    
    %% Negative Indices
    subplot(2,2,2,'replace')
    hold on
    notnan=find(~isnan(array.skewness_neg));
    filtskew=filtfilt(ones(3,1),1,array.skewness_neg(notnan))/3^2;
    plot(array.lat(notnan),filtskew)
    plot(array.lat,array.skewness_neg)
    scatter(lat_inter,zeros(size(lat_inter)),'ko','LineWidth',2)
    xlim([min(lat_inter)-5 max(lat_inter)+5])
    flag=true;
    count=0;
    
    while flag
        
        [centerlat centerlon button]=ginput(2);
        if isempty(centerlat) || sum(button==3)>0
            flag=false;
        else
            %% Set Positive Centers
            count=count+1;
            subplot(2,2,1)
            [minlat maxlat]=set_endpoints(array.lat,array.skewness_pos,centerlat(1));
            scatter([minlat maxlat],interp1(array.lat(notnan),array.skewness_pos(notnan),[minlat,maxlat]),'x','LineWidth',2)
            array.optpar_pos(count,:)=[centerlat(1) minlat maxlat];
            %% Set Negative Centers
            subplot(2,2,2)
            [minlat maxlat]=set_endpoints(array.lat,array.skewness_neg,centerlat(2));
            scatter([minlat maxlat],interp1(array.lat(notnan),array.skewness_neg(notnan),[minlat,maxlat]),'x','LineWidth',2)
            array.optpar_neg(count,:)=[centerlat(2) minlat maxlat];
                        
        end
        
    end
   
    if count>0
       save([outpath filesep sprintf('t%03d.mat',tracknum)],'array'); 
    end
    
end
    
    
    function [ tracks ] = trackconsec( )
    %     Sorts ascending/descending tracks in the correct order
    descend=rem(1+(0:126)*76,254);
    ascend=rem(0+(0:126)*76,254);
    ascend(ascend==0)=254;
    ascend=ascend([1:36 37:127]);
    tracks=[descend ascend];
    %     tracks=descend;
    end
    
    
    
end

