%% Test on synthetic data
[datay datah dataskew argument] = simulate_SSH([100 0 100]);
idx = datay < 100 & datay > -100;
LB = [0.01 -100  0.01];
UB = [1000 100 1000];
x0 = [50 50 50];
[optimal fval simskew] = jet_model_estimate(datay(idx),dataskew(idx),LB,UB,x0);
disp(optimal)

clf
subplot(1,2,1)
plot(datay,datah);
subplot(1,2,2)
hold on;
plot(datay,dataskew,'k-x');
plot(datay(idx),simskew,'b-x');
xlim([-500 500])

%% Test with real data from t001
load t001.mat
track.dist = [0 ; cumsum(sw_dist(track.lat,track.lon))];
track.skewness = skewness(track.sla,1,2);
latidx = track.lat >=-46.87 & track.lat <= -44.96;
% datay = linspace(min(track.dist(latidx)), max(track.dist(latidx)));
datay = track.dist(latidx);
% dataskew = interp1(track.dist,track.skewness,datay);
dataskew = track.skewness(latidx);
[optimal fval simskew] = jet_model_estimate(datay', ...
    dataskew',[50 1800 0.001],[200 2400 500],[75 2000 75]);

clf
hold on;
plot(track.dist(latidx),track.skewness(latidx)+1,'-x')
plot(datay,simskew+1,'k-x')
disp(optimal)