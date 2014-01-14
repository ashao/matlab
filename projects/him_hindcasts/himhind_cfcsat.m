inpath = '/ltraid3/ashao/uw-apl/models/offtrac/active_runs/hindcast/';
syear = 1936;
eyear = 1980;
years = syear:eyear;
nyears = length(years);


%% Extract data
hindcast.cfc11sat = zeros(nyears,210,360);

count4d = [12 1 inf inf];
count3d = [12 inf inf];
infile = [inpath 'cfcs.hindcast.geonly.0899.nc'];


for i=nyears:nyears
    start4d = [(i-1)*12 0 0 0];
    start3d = [(i-1)*12 0 0];
    fprintf('Year %d\n',years(i))    
    hindcast.cfc11sat(i,:,:)=squeeze(min( ...
        nc_varget(infile,'mn_cfc11',start4d,count4d)./...
        nc_varget(infile,'mn_cfc11_sat',start3d,count3d)));    
    
end

%%
outpath = '/ltraid3/ashao/uw-apl/figs/him_hindcasts/cfcs/';
% mkdir(outpath)
load metrics
m_proj('Lambert Conformal Conic','lon',[-240 -120],'lat',[20 60]);
m_proj('Mercator','lon',[-280 80],'lat',[-65 65]);
colormap(othercolor('BuDRd_12'))
for i = nyears
    clf; hold on;
%     m_pcolor(metrics.geolon.data,metrics.geolat.data,(1-hindcast.cfc12sat(i,:,:))*100);
%     shading flat
    title(sprintf('Year %d',years(i)));
    h = colorbar('SouthOutside');
    m_grid;
    m_coast('patch',[0.5 0.5 0.5]);
    xlabel(h,'Minimum CFC-12 Saturation')
%     pause
    print(gcf,'-dpng','-r200',[outpath sprintf('cfc11sat.%d.png',years(i))])    
end

    m_contourf(metrics.geolon.data,metrics.geolat.data, ...
        (1-hindcast.cfc12sat(i,:,:))*100,[-10 0:0.5:15],'LineColor','None')
    caxis([0 15]);
    title(sprintf('Year %d',years(i)));
    h = colorbar('SouthOutside');
    m_grid;
    m_coast('patch',[0.5 0.5 0.5]);
    xlabel(h,'Minimum CFC-12 Saturation')
%     pause
    print(gcf,'-dpng','-r200',[outpath sprintf('cfc11sat.%d.png',years(i))])    
end

%% Volume of mixed layer depth (NOT DONE CAREFULLY)
hindcast.cfc11sat_vol = zeros(nyears,1);
maxlon = -120;
minlon = -240;
maxlat = 60;
minlat = 20;

for i=1:nyears
    
    idx = squeeze(hindcast.cfc11sat(i,:,:))>200 & ...
        metrics.geolon.data < maxlon & metrics.geolon.data > minlon &...
        metrics.geolat.data < maxlat & metrics.geolat.data > minlat;
    hindcast.cfc11sat_vol(i) = sum(metrics.Ah.data(idx)'.* ...
        hindcast.cfc11sat(i,idx));
    
end

subplot(2,1,1)
plot(years,hindcast.cfc11sat_vol,'-x','LineWidth',2)
xlabel('Year')
ylabel('Volume of N. Pac. > 200m');
xlim([min(years) max(years)]);
grid on;

subplot(2,1,2)
plot(years,...
    (hindcast.cfc11sat_vol./mean(hindcast.cfc11sat_vol)-1)*100,'-x','LineWidth',2)
xlabel('Year')
ylabel('% differenc from average');
xlim([min(years) max(years)]);
grid on;

saveas(gcf,[outpath 'mldepth_vol.1948.2007.eps'],'epsc')