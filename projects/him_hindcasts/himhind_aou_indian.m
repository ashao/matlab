load ~/uw-apl/models/HIM/hindcast/hindcast.PV.mat
hindcast.temp = load('~/uw-apl/models/HIM/hindcast/hindcast.temp.mat');
hindcast.salt = load('~/uw-apl/models/HIM/hindcast/hindcast.salt.mat');
load metrics
%% 32S properties through core of mode water
latidx = abs(metrics.lath.data-(-32.5))<0.001;
layidx = 20; % Core of mode water
years = 1948:2007;

subplot(2,2,1)
data= double(squeeze(hindcast.PV(:,layidx,latidx,:)));
data = [data data(:,1:180)];
lonplot = metrics.geolon.data(latidx,:);
lonplot = [lonplot wrapTo360(lonplot(1:180))];

contourf(lonplot,years,log10(-data),-10.5:.1:-9.5,'LineColor','none')
caxis([-10.5 -9.5])
xlim([20 120])
cax=colorbar;
ylabel(cax,'log10(PV)')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
% title('SAMW 32S')

subplot(2,2,2)
data= double(squeeze(hindcast.depth(:,layidx,latidx,:)));
data = [data data(:,1:180)];
contourf(lonplot,years,data,'LineColor','none')
% caxis([11 13])
xlim([20 120])
cax=colorbar;
ylabel(cax,'depth')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
% title('SAMW 32S')

subplot(2,2,3)
data= double(squeeze(hindcast.temp.temp(:,layidx,latidx,:)));
data = [data data(:,1:180)];
contourf(lonplot,years,data,11:.1:13,'LineColor','none')
caxis([11 13])
xlim([20 120])
cax=colorbar;
ylabel(cax,'Temperature')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
% title('SAMW 32S')

subplot(2,2,4)
data= double(squeeze(hindcast.salt.salt(:,layidx,latidx,:)));
data = [data data(:,1:180)];
contourf(lonplot,years,data,34.5:.05:35.5,'LineColor','none')
% caxis([11 13])
xlim([20 120])
cax=colorbar;
ylabel(cax,'Salinity')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
% title('SAMW 32S')
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/him_hindcasts/SAMW.32S.hovmoeller.eps','epsc')
%%
%% 32S properties through core of intermediate water
latidx = abs(metrics.lath.data-(-32.5))<0.001;
layidx = 27; % Core of mode water
years = 1948:2007;

subplot(2,2,1)
data= double(squeeze(hindcast.PV(:,layidx,latidx,:)));
data = [data data(:,1:180)];
contourf(lonplot,years,log10(-data),-10:.025:-9,'LineColor','none')
caxis([-9.8 -9.5])
xlim([20 120])
cax=colorbar;
ylabel(cax,'log10(PV)')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
title('SAMW 32S')

subplot(2,2,2)
data= double(squeeze(hindcast.depth(:,layidx,latidx,:)));
data = [data data(:,1:180)];
contourf(lonplot,years,data,500:50:1000,'LineColor','none')
% caxis([11 13])
xlim([20 120])
cax=colorbar;
ylabel(cax,'Depth')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
% title('SAMW 32S')

subplot(2,2,3)
data= double(squeeze(hindcast.temp.temp(:,layidx,latidx,:)));
data = [data data(:,1:180)];
contourf(lonplot,years,data,3:0.25:8,'LineColor','none')
% caxis([11])
xlim([20 120])
cax=colorbar;
ylabel(cax,'Temperature')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
% title('SAMW 32S')

subplot(2,2,4)
data= double(squeeze(hindcast.salt.salt(:,layidx,latidx,:)));
data = [data data(:,1:180)];
contourf(lonplot,years,data,34:.01:34.5,'LineColor','none')
% caxis([11 13])
xlim([20 120])
cax=colorbar;
ylabel(cax,'Salinity')
colormap(flipud(othercolor('Blues9')))
% shading flat
xlabel('Longitude')
% ylabel('Time')
% title('SAMW 32S')
saveas(gcf,'/ltraid3/ashao/uw-apl/figs/him_hindcasts/AAIW.32S.hovmoeller.eps','epsc')