inpath = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/';
syear = 1948;
eyear = 2007;
years = syear:eyear;
nyears = length(years);


%% Extract data
hindcast.mldepth = zeros(nyears,210,360);
start4d = [0 0 0 0];
count4d = [inf 2 inf inf];

for i=1:nyears
   
    fprintf('Year %d\n',years(i))
    infile = [inpath sprintf('ocean_month.%d.nc',years(i))];
    hindcast.mldepth(i,:,:)=squeeze(max(sum( ...
        nc_varget(infile,'h',start4d,count4d),2)));    
    
end

%%
outpath = '/ltraid3/ashao/uw-apl/figs/him_hindcasts/mldepth/';
load metrics
m_proj('Mercator','lon',[-240 -120],'lat',[20 65]);
colormap(othercolor('BuDRd_12'))
for i = 1:nyears
    clf; hold on;
    m_pcolor(metrics.geolon.data,metrics.geolat.data,hindcast.mldepth(i,:,:));
    shading flat
    m_contour(metrics.geolon.data,metrics.geolat.data, ...
        hindcast.mldepth(i,:,:),[150 200],'LineColor','Black','LineWidth',2)
    caxis([0 300]);
    title(sprintf('Year %d',years(i)));
    h = colorbar('SouthOutside');
    m_grid;
    m_coast('patch',[0.5 0.5 0.5]);
    xlabel(h,'Maximum mixed layer depth')
    print(gcf,'-dpng','-r200',[outpath sprintf('mldepth.%d.png',years(i))])    
end

%% Volume of mixed layer depth (NOT DONE CAREFULLY)
hindcast.mldepth_vol = zeros(nyears,1);
maxlon = -120;
minlon = -240;
maxlat = 60;
minlat = 20;

for i=1:nyears
    
    idx = squeeze(hindcast.mldepth(i,:,:))>200 & ...
        metrics.geolon.data < maxlon & metrics.geolon.data > minlon &...
        metrics.geolat.data < maxlat & metrics.geolat.data > minlat;
    hindcast.mldepth_vol(i) = sum(metrics.Ah.data(idx)'.* ...
        hindcast.mldepth(i,idx));
    
end

subplot(2,1,1)
plot(years,hindcast.mldepth_vol,'-x','LineWidth',2)
xlabel('Year')
ylabel('Volume of N. Pac. > 200m');
xlim([min(years) max(years)]);
grid on;

subplot(2,1,2)
plot(years,...
    (hindcast.mldepth_vol./mean(hindcast.mldepth_vol)-1)*100,'-x','LineWidth',2)
xlabel('Year')
ylabel('% differenc from average');
xlim([min(years) max(years)]);
grid on;

saveas(gcf,[outpath 'mldepth_vol.1948.2007.eps'],'epsc')