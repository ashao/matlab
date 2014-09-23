hinddir = '/ltraid4/ashao/HIM/hyak_store/HINDCAST/';
flux_files = dir([hinddir 'flux_month*.nc']);
nyears = length(flux_files);
nmonths = nyears*12;
flux.lat = nc_varget([hinddir flux_files(1).name],'lat');
flux.lon = nc_varget([hinddir flux_files(1).name],'lon');
nlat = length(flux.lat);
nlon = length(flux.lon);
flux.landmask = ~(logical(nc_varget([hinddir flux_files(1).name],'land_mask')));
flux.taux = zeros(nyears*12,nlat,nlon);
flux.tauy = zeros(nyears*12,nlat,nlon);
%%
for year = 1:nyears
   
    fprintf('Year %d\n',year)
    infile = [hinddir flux_files(year).name];
    sidx = (year-1)*12+1;
    eidx = sidx+11;
    flux.taux(sidx:eidx,:,:) = nc_varget(infile,'tau_x');
    flux.tauy(sidx:eidx,:,:) = nc_varget(infile,'tau_y');
    
end
%%
flux.curltau = zeros(size(flux.taux));
flux = structfun(@double,flux,'UniformOutput',false);
[latgrid longrid] = meshgrat(flux.lat,flux.lon);
for mon = 1:nmonths

    fprintf('Month %d/%d\n',mon,nmonths)

    [null null grady null] = gradientm(latgrid,longrid,squeeze(flux.taux(mon,:,:))./1025./sw_f(latgrid).*flux.landmask);
    [null null null gradx] = gradientm(latgrid,longrid,squeeze(flux.tauy(mon,:,:))./1025./sw_f(latgrid).*flux.landmask);
    flux.curltau(mon,:,:) = grady-gradx;
            
end

%% Truncate for just the Indian ocean basin
indian.latrange = [-60 -20];
indian.lonrange = [20 160];
indian.mask = latgrid > min(indian.latrange) & latgrid < max(indian.latrange) & ...
    longrid < max(indian.lonrange) & longrid > min(indian.lonrange);
indian.landmask = logical(flux.landmask(indian.mask));
lat = unique(latgrid(indian.mask));
lon = unique(longrid(indian.mask));
[indian.latgrid indian.longrid] = meshgrat(lat,lon);

indian.curltau=zeros([720 size(indian.latgrid)]);
indian.taux=zeros([720 size(indian.latgrid)]);
indian.curltau(:) = flux.curltau(:,indian.mask);
indian.taux(:) = flux.taux(:,indian.mask);
load('~/uw-apl/data/ncep/sam_from_ncep.mat');
%% Calculate EOF for Ekman pumping

indian.annual.curltau = nan([nyears size(indian.latgrid)]);
for year = 1:nyears
    sidx = (year-1)*12+1;
    eidx = sidx+11;
    indian.annual.curltau(year,:,:) = mean(indian.curltau(sidx:eidx,:,:));
    
end
indian.annual.curltau = indian.annual.curltau - repmat(mean(indian.annual.curltau),[nyears 1 1]);

data = indian.annual.curltau(:,:);
data = data(:,logical(indian.landmask));
data = detrend(data')';
[U S V] = svd(data);


mode = 1;
clf
eof = nan(size(indian.latgrid));
eof(indian.landmask) = S(mode,mode)*V(:,mode);
subplot(3,2,[1 2 3 4])
worldmap(indian.latrange,indian.lonrange)
contourfm(indian.latgrid,indian.longrid,eof*1e6,-2:0.1:2,'LineStyle','none')
caxis([-2 2])
cax = colorbar;
ylabel(cax,'Ekman pumping [10^-6 1/(m s)]')
geoshow('landareas.shp', 'FaceColor', [0.5 0.5 0.5])
pcs = U*S;

subplot(3,2,[5 6]); hold on

area(datenum(1948:2007,6,15),pcs(:,mode)./std(pcs(:,mode)))
plot(sam_ncep.time,smooth(sam_ncep.pc./std(sam_ncep.pc),12,'mean'),'k-','LineWidth',2)

ylim([-2 2]);
xlim(datenum([1948 2008],1,1))
datetick('KeepTicks','KeepLimits');
legend('PC','SAM')
