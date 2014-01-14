load C:\USers\ashao\data\aviso\t012.mat
track.dist = [0 ; cumsum(sw_dist(track.lat,track.lon))];
subplot(2,1,2);
m_plot(track.lon,track.lat,'LineWidth',2,'Color','black')
plot_fronts
m_coast('patch',[0 0 0]);
m_grid;
for sidx = 1:20:(1060);
    
    eidx = sidx + 100-1;    
    if eidx > 1060
        break
    end
    counter = counter + 1;
    track.skewness = skewness(track.sla(:,sidx:eidx),1,2);
    subplot(2,1,1);cla;hold on
    plot(track.lat,track.skewness);
    plot(track.lat,smooth(track.skewness,10,'mean'),'LineWidth',3,'Color','black');
    xlim([-60 -40])
    ylim([-2 2])
    grid on
    pause
end