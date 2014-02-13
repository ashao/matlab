
parameters = [1 -60 1];
[simy track_heights simskew argument] = simulate_SSH(parameters);
outpath = 'C:\Users\ashao\Documents\GitHub\matlab\projects\safaltimetry\poster\';

%%
clf
subplot(3,3,[1 2 3])
plot(simy,track_heights(1:10,:),'-kx','LineWidth',2)
xlim([-66 -54])
grid on;
xlabel('Latitude')
ylabel('SSH')
title('Measured SSH')


sidx = 43;
nidx = 57;
midx = 50;

subplot(3,3,4)
hist(track_heights(:,sidx))
xlabel('SSH')
ylabel('Count')
title('SSH Distribution at -61')
xlim([-2 2])
grid on;


subplot(3,3,5)
hist(track_heights(:,midx))
xlabel('SSH')
ylabel('Count')
xlim([-2 2])
grid on
title('SSH Distribution at -60')


subplot(3,3,6)
hist(track_heights(:,nidx))
xlabel('SSH')
ylabel('Count')
xlim([-2 2])
grid on
title('SSH Distribution at -59')

subplot(3,3,7:9)
plot(simy,skewness(track_heights),'LineWidth',2)
xlabel('Latitude')
ylabel('Skewness')
ylim([-3 3]);xlim([-66 -54]); grid on;

export_fig(gcf,[outpath 'schematic.eps'],'-eps','-transparent','-cmyk')
