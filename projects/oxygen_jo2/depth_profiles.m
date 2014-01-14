%% Create structure with paths to the various model runs
infile.cfc = '/home/ashao/uw-apl/models/offtrac/run_archive/cfcsf6.normalyear/cfcs.sf6.control.0899.nc';
infile.age = '/home/ashao/uw-apl/models/offtrac/run_archive/age.23999.nc';
infile.oxygen = '/home/ashao/uw-apl/models/offtrac/active_runs/oxygen/oxygen.saturation.1199.nc';

%% Files dealing with things like temp, salinity, etc.
infile.ts = '/home/ashao/uw-apl/data/offtrac/input/normalyear/ts.nc';
infile.h = '/home/ashao/uw-apl/data/offtrac/input/normalyear/H-clim.nc';

% Calculate annual averages
offtrac.temp = mean(nc_varget(infile.ts,'temp'));
offtrac.salt = mean(nc_varget(infile.ts,'salt'));
offtrac.depth = cumsum(squeeze(mean(nc_varget(infile.h,'h'))));
offtrac.depth(1,:,:) = 0; % Set the first depth to 0

% Calcuate density
offtrac.pden2000 = sw_pden(offtrac.salt,offtrac.temp,2000,2000);

%% Get the last year of oxygen data

% Note netcdf uses C format for indexing e.g. the first element of an 
% array has an index of zero
start4d = [1200-12-1 0 0 0]; 
count4d = [12 inf inf inf]; % inf means read all along that dimension

offtrac.o2 = squeeze(mean( ... 
    nc_varget(infile.oxygen,'mn_oxygen',start4d,count4d)));
offtrac.o2_sat = squeeze(mean( ... 
    nc_varget(infile.oxygen,'mn_o2sat',start4d,count4d)));
offtrac.jo2 = squeeze(mean( ... 
    nc_varget(infile.oxygen,'mn_jo2',start4d,count4d)));
% Calculate AOU
offtrac.aou = offtrac.o2_sat-offtrac.o2;

%% Last year of CFC data
start4d = [888-12-1 0 0 0];
count4d = [12 inf inf inf];
offtrac.cfc11 = squeeze(mean(nc_varget(infile.cfc,'mn_cfc11',start4d,count4d)));
offtrac.cfc12 = squeeze(mean(nc_varget(infile.cfc,'mn_cfc12',start4d,count4d)));

%% Lastly load ideal age tracer (Already annually averaged)
start4d = [1200-1 0 0 0];
count4d = [1 inf inf inf];
offtrac.age = nc_varget(infile.age,'mn_age',start4d,count4d);

%% Plot for N.Pacific 152.5W and 30.5N
load /home/ashao/matlab/metrics.mat % HIM grid information
plotlat = 30.5;
plotlon = -152.5;
plotidx = floor(metrics.geolat.data*10)/10==30.5 & ...
    floor(metrics.geolon.data*10)/10==-152.5;
voltogravconv = 1000./offtrac.pden2000; % Convert mol/L to mol/kg
figure(1) % Oxygen related plots

subplot(1,3,1)
moltomicromol = 1e3; % NEED TO CHECK CONVERSION
plot(offtrac.o2(:,plotidx).*voltogravconv(:,plotidx)*moltomicromol,offtrac.depth(:,plotidx))
xlabel('Oxygen (10^{-6} mol kg^{-1})')
ylabel('Depth')
set(gca,'ydir','reverse')

subplot(1,3,2)
moltomicromol = 1e3; % NEED TO CHECK CONVERSION
plot(offtrac.jo2(:,plotidx).*voltogravconv(:,plotidx)*moltomicromol,offtrac.depth(:,plotidx))
xlabel('Oxygen (10^{-6} mol kg^{-1} yr^{-1})')
ylabel('Depth')
set(gca,'ydir','reverse')

subplot(1,3,3)
moltomicromol = 1e3; % NEED TO CHECK CONVERSION
plot(offtrac.aou(:,plotidx).*voltogravconv(:,plotidx)*moltomicromol,offtrac.depth(:,plotidx))
xlabel('AOU (10^{-6} mol kg^{-1})')
ylabel('Depth')
set(gca,'ydir','reverse')

figure(2) % CFC/Age plots
subplot(1,3,1)
plot(offtrac.cfc11(:,plotidx).*voltogravconv(:,plotidx),offtrac.depth(:,plotidx))
xlabel('CFC-11 (10^{-6} mol kg^{-1})')
ylabel('Depth')
ylim([0 3000])
set(gca,'ydir','reverse')

subplot(1,3,2)
plot(offtrac.cfc12(:,plotidx).*voltogravconv(:,plotidx),offtrac.depth(:,plotidx))
xlabel('CFC-12 (10^{-6} mol kg^{-1}')
ylabel('Depth')
set(gca,'ydir','reverse')
ylim([0 3000])

subplot(1,3,3)
plot(offtrac.age(:,plotidx),offtrac.depth(:,plotidx))
xlabel('Ideal Age (yr)')
ylabel('Depth')
set(gca,'ydir','reverse')
ylim([0 3000])