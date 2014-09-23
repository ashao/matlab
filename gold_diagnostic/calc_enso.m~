load metrics

%% Nino 3.4 region
ninomask = metrics.geolon.data <= -120 & metrics.geolon.data >= -170 & ...
    abs(metrics.geolat.data)<=5;
areawt = metrics.Ah.data;
areawt(~ninomask)=0;
areawt = areawt/sum(makevec(areawt));

%% SST from ice_month
icefile = '/ltraid1/ashao/gold/HINDCAST/ice_month.hindcast.nc';
sst = nc_varget(icefile,'SST');
%%
areawttemp(1,:,:)=areawt;
areawtgrid=repmat(areawttemp,[length(sst) 1 1]);
nino34 = sum(squeeze(sst(:,ninomask)).*areawtgrid(:,ninomask),2);
clf
plot((1:length(sst))/12+1948,(nino34-mean(nino34))/std(nino34))
ylim([-3 3]);
grid on