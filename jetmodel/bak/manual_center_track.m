function [ centers_init ] = manual_center_track(allarray,outpath,fronts,start_track,centers_init)

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

for trackn=start_track:length(tracks)
    
    [sidx null]=size(centers_init);
    delete(trackh)
    tracknum=tracks(trackn);
    clear array;
    array.lat=allarray.latarray(tracknum,:);
    array.lon=allarray.lonarray(tracknum,:);
    array.skewness=allarray.skewness(tracknum,:);
    arenan=find(isnan(array.skewness));
    array.lat(arenan)=[];
    array.lon(arenan)=[];
    array.skewness(arenan)=[];
    lon_inter=[];
    lat_inter=[];
    subplot(2,2,[3 4])
    trackh=m_plot(array.lon,array.lat);
    
    title(sprintf('Track %d (%d/%d)',tracknum,trackn,length(tracks)))
    
    for front=1:nfronts
        [templon templat]=intersections(fronts(front).lon,fronts(front).lat,array.lon,array.lat);
        lon_inter=[lon_inter ; templon];
        lat_inter=[lat_inter ; templat];
    end
    
    filtskew=filtfilt(ones(2,1),1,array.skewness)/2^2;
    subplot(2,2,1,'replace')
    hold on
    plot(array.lat,filtskew)
    plot(array.lat,array.skewness,'LineWidth',1.5)
    scatter(lat_inter,zeros(size(lat_inter)),'o','filled','LineWidth',2)
    xlim([-50 -20])
    xlim([min(lat_inter)-5 max(lat_inter)+5])
%     ylim([-1.5 1.5])
    grid on
    
    count = 0 ;
    
    flag = true;    
    
    while flag
        clear centerlat
        [centerlat lon button]=ginput(1);
        if isempty(centerlat) || button~=1
            flag = false;
        else            
            count=count+1;
            [minlat maxlat]=set_endpoints(array.lat,array.skewness,centerlat);
            scatter([minlat maxlat] ,interp1(array.lat,array.skewness,[minlat,maxlat]),'x','LineWidth',2)

            subplot(2,2,[3 4])
            hold on;
            front_lat=mean([minlat maxlat]);
            front_lon=interp1(array.lat,array.lon,front_lat);
            m_plot(front_lon,front_lat,'x','LineWidth',2);
            centers_init(sidx+1,:)=[front_lon front_lat];
            
            subplot(2,2,2)
            plot(array.lat,array.skewness)
            grid on
            xlim([minlat maxlat])
            array.optpar( count,: )=[ centerlat minlat maxlat ];
        end
    end
    disp(sprintf('Track %d with %d centers',tracknum,count))
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

