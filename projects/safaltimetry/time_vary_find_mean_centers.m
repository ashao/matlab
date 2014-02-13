function  opt_track  = find_mean_centers( track )

[npts nfiles] = size(track.sla);
winlength = 150;
winslide = 25;
% track.dist = [0 ; cumsum(sw_dist(track.lat,track.lon))];
% subplot(2,1,2);
% plot_fronts
% plotm(track.lat,track.lon,'LineWidth',2,'Color','black')
% % m_coast('patch',[0 0 0]);
% % m_grid('ytick',-70:5:-30)
%
% % Choose the range based on all the
% subplot(2,1,1);cla;hold on
% randidx = randi(nfiles,200,1);
% track.skewness = skewness (track.sla(:,randidx),1,2);
% track.smoothskew=smooth(track.skewness,10,'mean')
% plot(track.lat,track.skewness);
% plot(track.lat,track.smoothskew,'LineWidth',2,'Color','black');
%
% xlim([-65 -35])
% ylim([-2 2])
% grid on
% [latrange null]=ginput(2);
latrange = track.latrange;

minlat = min(latrange);
maxlat = max(latrange);

trunc.idx = track.lat >= minlat & track.lat <= maxlat;
counter = 0;

for sidx = 1:winslide:(1060);
    
    eidx = sidx + winlength-1;
    if eidx > 1060
        break
    end
    counter = counter + 1;
    
    track.skewness = skewness(track.sla(:,sidx:eidx),1,2);
    track.std = nanstd(track.sla(:,sidx:eidx),1,2);
    track.smoothskew=smooth(track.skewness,10,'mean');
    
    trunc.lat = track.lat(trunc.idx);
    trunc.skewness = track.skewness(trunc.idx);
    trunc.smoothskew =track.smoothskew(trunc.idx);
    trunc.std = track.std(trunc.idx);
    trunc.dist = track.dist(trunc.idx);
    
    [maxtab mintab] = peakdet(trunc.smoothskew,0.2);
    
    if ~isempty(maxtab) & ~isempty(mintab)
        maxtab(:,3) = trunc.lat(maxtab(:,1));
        mintab(:,3) = trunc.lat(mintab(:,1));
        
        %
        %         scatter(maxtab(:,3),maxtab(:,2),'.');
        %         scatter(mintab(:,3),mintab(:,2),'.');
        
        maxtab = sortrows(maxtab,3);
        mintab = sortrows(mintab,3);
        
        inrange = maxtab(:,3) >= minlat & maxtab(:,3) <= maxlat;
        maxtab = maxtab(inrange,:);
        
        inrange = mintab(:,3) >= minlat & mintab(:,3) <= maxlat;
        mintab = mintab(inrange,:);
        
        
        while maxtab(1,3) > mintab(1,3)
            mintab(1,:) = [];
        end
        
        ntrans = min([size(maxtab,1) size(mintab,1)]);
        if ntrans > 0
            
            maxtab=maxtab(1:ntrans,:);
            mintab=mintab(1:ntrans,:);
            
            opt_track(counter).time = nanmean(nanmean(track.time(:,sidx:eidx)));
            
            for i = 1:ntrans
                
                transidx = maxtab(i,1):mintab(i,1);
                opt_track(counter).tracknum = track.tracknum;
                opt_track(counter).skewness{i} = trunc.skewness(transidx);
                opt_track(counter).std{i} = mean(trunc.std(transidx));
                opt_track(counter).lat{i} = trunc.lat(transidx);
                opt_track(counter).dist{i} = trunc.dist(transidx);
                %         opt_track(counter).lon{i}= trunc.lon(transidx);
                %         plot(opt_track(counter).lat{i},opt_track(counter).skew{i},'LineWidth',4)
            end
        end
    else
        opt_track = [];
        
    end
    
end
%         pause(1)
end

