
parameters = [1 -60 1];
[simy track_heights simskew argument] = simulate_SSH(parameters);
outpath = 'C:\Users\ashao\Documents\GitHub\matlab\projects\safaltimetry\paper\';

%%
clf
subplot(3,1,1)
plot(simy,track_heights(1:10,:),'-kx','LineWidth',2)
xlim([-66 -54])
grid on;
xlabel('Latitude')
ylabel('SSH')
title('(a)')
% title('Synthetic SSH from Gaussian meandering jet model')
box on;
sidx = 43;
nidx = 57;
midx = 50;

[Ns Xs] = hist(track_heights(:,sidx),20);
[Nm Xm] = hist(track_heights(:,midx),20);
[Nn Xn] = hist(track_heights(:,nidx),20);
subplot(3,1,2);hold on;
plot(Xs,Ns./sum(Ns),'r','LineWidth',2)
plot(Xm,Nm./sum(Nm),'k','LineWidth',2)
plot(Xn,Nn./sum(Nn),'b','LineWidth',2)
xlim([-1.5 1.5])
ylim([0 0.25])
xlabel('SSH')
ylabel('PDF')
grid on
legend('61\circ S','60\circ S','59\circ S')
title('b')
box on;

subplot(3,1,3)
plot(simy,skewness(track_heights),'k','LineWidth',2)
xlabel('Latitude')
ylabel('Skewness')
ylim([-2 2]);xlim([-66 -54]); grid on;
title('c')
box on;
export_fig(gcf,[outpath 'schematic.eps'],'-eps','-transparent','-cmyk')
