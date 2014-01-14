slapath = 'C:\Users\ashao\data\aviso\';
files = dir([slapath '*.mat']);
ntracks = length(files)

for i=[10 12 14 16]
    clf; hold on;
    load([slapath files(i).name]);
    track.mean = nanmean(track.sla,2);
    track.std = nanstd(track.sla,0,2);
    track.skewness = skewness(track.sla,0,2);    
    
    track.skewness(isnan(track.mean))=NaN;
    
    plot(track.lat,track.skewness,'-k','LineWidth',2)
    plot(track.lat,abs(track.mean./track.std),'-b','LineWidth',2)
    xlim([-70 -30])
    ylim([-1 1])
    grid on
    title(files(i).name)
    pause
    
    
   
    
end