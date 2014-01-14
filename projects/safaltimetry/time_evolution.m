trackpath = '/ltraid4/aviso/alongtrack/tracks/';
load([trackpath 't005.mat']);
%%
plotidx = track.lat < -52.5 & track.lat > -53.5;

startidx = 1;
endidx = 100;
clf; hold on;
counter = 0;
zerocross = [];
timecross = [];
while endidx <= 1060
    counter = counter+1;
    timeidx = startidx:endidx;
    skew = skewness(track.sla(plotidx,timeidx),1,2);
    plot(track.lat(plotidx),skew)
    skew = sort(skew,1,'descend');
    
    zerocross(counter) = interp1(skew,track.lat(plotidx),0,'linear','extrap');
    timecross(counter) = nanmean(nanmean(track.time(plotidx,timeidx)));
    startidx = startidx+20;
    endidx = startidx + 200;   
    pause(1)
end
