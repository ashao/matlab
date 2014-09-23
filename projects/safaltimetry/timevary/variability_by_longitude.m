load ~/uw-apl/projects/saf_altimetry/timevary.mat
%%
subplot(2,1,1)
varidx = std(timevary.south.lat');
[null sortidx] = sort(timevary.south.lon(:,1));
plot(timevary.south.lon(sortidx,1),smooth(varidx(sortidx),10,'mean'))
subplot(2,1,2)
worldmap([-60 60],[0 360])

geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])
%%
