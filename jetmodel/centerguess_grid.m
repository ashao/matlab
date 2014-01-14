function array = centerguess_grid(array,front0)
% Set initial conditions for optcruncher
clf

idx=find(front0.lon>=0 & front0.lon<360);
front0.lon=front0.lon(idx);
front0.lat=front0.lat(idx);

subplot(2,1,1)
plot(front0.lon,front0.lat); hold on;

[nlon nlat]=size(array.skewness);
latweight=0.85;
skewweight=-0.15;
array.initial=zeros(nlon,3)*NaN;

%% Guess center at each longitude
for loni=1:nlon
    
    %% Interpolate mean position of front
    frontlat=interp1(front0.lon,front0.lat,array.lon(loni));
    %% Truncate input skewness to area around guess;
    minlat=frontlat-10;
    maxlat=frontlat+10;
    latidx=find(array.lat>minlat & array.lat<maxlat);
    %     skewgrid=array.skewness(:,latidx);
    latvec=array.lat(latidx);
    
    %% Calculate Peaks
    skewvec=moving_average(array.skewness(loni,latidx),1);
    [maxtab mintab]=peakdet(skewvec,0.1);
    
    if ~isempty(maxtab) & ~isempty(mintab)
        
        peakidx=maxtab(:,1);
        peakskew=maxtab(:,2);
        peaklat=latvec(peakidx);
        
        valleyidx=mintab(:,1);
        valleyskew=mintab(:,2);
        valleylat=latvec(valleyidx);
        
        
        canprop=zeros(length(peakidx),2)*NaN;
        lattab=zeros(length(peakidx),30)*NaN;
        
        %% Calculate property of each candidate peak
        % canprop = [ latdiff skewchange ]
        for i=1:length(peakskew)
            nextvalleyidx=min(find( (valleylat-peaklat(i))>0 ));
            
            if ~isempty(nextvalleyidx)
                skewrange=skewvec(peakidx(i):valleyidx(nextvalleyidx));
                latrange=latvec(peakidx(i):valleyidx(nextvalleyidx));
                [null minskewidx]=min(abs( skewrange ) );
                %             whos
                canprop(i,1)=abs( frontlat-latrange(minskewidx) );
                canprop(i,2)=peakskew(i)-valleyskew(nextvalleyidx);
                
                lattab(i,1:length(latrange))=latrange;
            end
        end
        
        score=latweight*(canprop(:,1)./max(canprop(:,1)))+skewweight*(canprop(:,2)./max(canprop(:,2)));
        [null minidx]=min(score);
        frontguess=nanmean(lattab(minidx,:));
        array.initial(loni,:)=[0.5 frontguess 0.5];
        subplot(2,1,2)
        plot(latvec,skewvec); hold on;
        scatter(frontguess,0,'o')
        hold off;
        subplot(2,1,1)
        scatter(array.lon(loni),array.initial(loni,2));
        drawnow;
    end
    
    
end

end


