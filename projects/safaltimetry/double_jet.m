jet1.param = [1 -0.5 1];
jet2.param = [1 0.5 1];
simy = linspace(-5,5,2000);

jet1.trackheight = simulate_SSH(jet1.param,simy);
jet2.trackheight = simulate_SSH(jet2.param,simy);

simheight = jet1.trackheight+jet2.trackheight;
simskew = skewness(simheight);
clf; hold on;
plot(simy,simskew)
plot(simy,skewness(jet1.trackheight),'g')
plot(simy,skewness(jet2.trackheight),'r')
ylim([-1.5 1.5])
