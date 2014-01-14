load C:\USers\ashao\data\aviso\t014.mat
track.dist = [0 ; cumsum(sw_dist(track.lat,track.lon))];

idx_mat = reshape(1:1060,106,10);
L = zeros(100,1);
sigma = zeros(100,1);
y0 = zeros(100,1);
time = zeros(100,1);

sidx = 1;
counter = 0;
for sidx = 1:20:(1060);
    
    eidx = sidx + 100-1;    
    if eidx > 1060
        break
    end
    counter = counter +1;
    track.skewness = skewness(track.sla(:,sidx:eidx),1,2);
    plot(track.lat,smooth(track.skewness,6,'mean'));
    xlim([-60 -40])
    ylim([-2 2])
    bounds = ginput(2);
    latidx = track.lat >= min(bounds(:,1)) & track.lat <= max(bounds(:,1));
    datay = track.dist(latidx);
    dataskew = track.skewness(latidx);
    LB = [20 mean(datay)-500 20];
    x0 = [100 mean(datay) 100];
    UB = [250 mean(datay)+500 250];
    
    
    [optimal fval simskew] = jet_model_estimate(datay', ...
        dataskew',LB,UB,x0);    
    L(counter) = optimal(1);
    y0(counter) = interp1(track.dist,track.lat,optimal(2));
    sigma(counter) = optimal(3);
    time(counter) = nanmean(nanmean(track.time(:,sidx:eidx)));
end
%%
delidx = time==0;
L(delidx)=[];
x0(delidx)=[];
sigma(delidx)=[];
time(delidx) = [];

%%
m_proj('Mercator','lon',[80 180],'lat',[-65 -30])
load t012.mat
load t016.mat
clf
subplot(1,2,1)
load C:\USers\ashao\data\aviso\t012.mat
hold on;
m_plot(track.lon,track.lat,'r-x');
load C:\USers\ashao\data\aviso\t016.mat
m_plot(track.lon,track.lat,'b-x');
m_coast('patch',[ 0 0 0])
m_grid;
legend({'t012','t016'})

subplot(1,2,2); hold on;
time = t012.time+datenum(1950,1,1);
plot(time,t012.y0,'r-x')
time = t016.time + datenum(1950,1,1);
plot(time,t016.y0)
legend('t012','t016')
datetick
title('Jet center over time')


