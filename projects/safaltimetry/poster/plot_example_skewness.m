inpath = 'C:\Users\ashao\Data\vxxc_matlab\';
files = dir([inpath filesep 't*.mat']);
load([inpath files(1).name])

%%
clf; hold on;
box on
h=patch([-58.5 -57.5 -57.5 -58.5],[-1 -1 1 1],[0 0 0]+0.5);
alpha(h,0.2)
h=patch([-47.5 -46 -46 -47.5],[-1 -1 1 1],[0 0 0]+0.5);
alpha(h,0.2)
plot(track.lat,skewness(track.sla),'k','LineWidth',2)
xlim([-60 -40])
grid on;
ylim([-1 1])
% line([-48 -48],[-1 1],'Color','k','LineWidth',2)
title('Altimeter Track 1')
xlabel('Latitude')
ylabel('Skewness')


pbaspect([3 1 1])
%%
outpath = 'C:\Users\ashao\Documents\GitHub\matlab\projects\safaltimetry\poster\';
export_fig([outpath 't001_skewness.eps'],'-eps','-rgb')