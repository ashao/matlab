%% Plot annually averaged temperature from year 456 of HIM spinup along 30.5N

load metrics.mat
inpath.ocean_month = '/ltraid4/ashao/HIM/himw/him_sis/ocean_month.nc';
himclim.temp = squeeze(mean( nc_varget(inpath.ocean_month,'temp')));
himclim.salt = squeeze(mean( nc_varget(inpath.ocean_month,'salt')));
himclim.depth = cumsum(squeeze(mean( nc_varget(inpath.ocean_month,'h'))));


%% 
latidx = find(metrics.lath.data==30.5);
longrid = repmat(metrics.geolon.data(latidx,:),[49 1]);
colormap(flipud(othercolor('RdYlBu9')))

% Set some plotting commands common to both subplots
toplim = 1000;
labels.x = 'Longitude';
labels.y = 'Depth (m)';
labels.cax = '\theta_{2000} (\circ C)';
plotcmds = ['caxis([0 25]);' ...
    'ylabel(labels.y);xlabel(labels.x);' ...
    'set(gca,''ydir'',''reverse'');' ...
    'xlim([-240 -110]);'];

% Actually make the plots
% NOTE: pcolor shows the actual grid point by grid point temperature
subplot(2,1,1)
pcolor(longrid, squeeze(himclim.depth(:,latidx,:)), ...
    squeeze(himclim.temp(:,latidx,:)))
ylim([0 toplim])
eval(plotcmds)

subplot(2,1,2)
pcolor(longrid, squeeze(himclim.depth(:,latidx,:)), ...
    squeeze(himclim.temp(:,latidx,:)))
ylim([toplim 6000])
eval(plotcmds)

cax = colorbar('SouthOutside');
xlabel(cax,labels.cax)