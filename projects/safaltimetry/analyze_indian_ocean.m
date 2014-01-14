trackpath = 'C:\Users\ashao\data\aviso\';
trackfiles = dir([trackpath '*.mat']);

trackidx = zeros(254,1);
for i = 1:length(trackfiles)
    
    load([trackpath filesep trackfiles(i).name]);
    [x0 y0] = intersections(fronts.saf.lon,fronts.saf.lat,track.lon,track.lat);
    trackidx(i) = any(x0>=80 & x0<=180);
    
end

%%
m_proj('Mercator','lon',[0 360],'lat',[-70 -30]);
mmapcmds = ['m_coast(''patch'',[0 0 0]);' 'm_grid;'];

jet.lon = zeros(length(trackfiles));
jet.lat = zeros(length(trackfiles));
jet.L = zeros(length(trackfiles));
jet.sigma = zeros(length(trackfiles));
counter = 0;
for i=1:length(trackfiles)
    if trackidx(i)
        counter = counter + 1;
        clf
        load([trackpath filesep trackfiles(i).name]);
        [x0 y0] = intersections(fronts.saf.lon,fronts.saf.lat,track.lon,track.lat);
        track.dist = [0 ; cumsum(sw_dist(track.lat,track.lon))];
        subplot(2,1,1); hold on;
        skew = skewness(track.sla,1,2);
        plot(track.lat,skew);
        scatter(y0,interp1(track.lat,skew,y0))
        ylim([-1 1])
        xlim([-60 -40])
        title(sprintf('Track %02d',track.tracknum))
        subplot(2,1,2); hold on;
        m_plot(fronts.saf.lon,fronts.saf.lat,'k')
        m_plot(track.lon,track.lat)
        eval(mmapcmds)
        
        bounds = ginput(2);
        if ~isempty(bounds)
        track.skewness = skewness(track.sla,1,2);
        latidx = track.lat >= min(bounds(:,1)) & track.lat <= max(bounds(:,1));
        datay = track.dist(latidx);
        dataskew = track.skewness(latidx);
        LB = [20 mean(datay)-500 20];
        x0 = [100 mean(datay) 100];
        UB = [250 mean(datay)+500 250];
        
        
        [optimal fval simskew] = jet_model_estimate(datay', ...
            dataskew',LB,UB,x0);
        jet.L(counter) = optimal(1);
        jet.lat(counter) = interp1(track.dist,track.lat,optimal(2));
        jet.lon(counter) = interp1(track.lat,track.lon,jet.lat(counter));
        jet.sigma(counter) = optimal(3);
        end
    end
end