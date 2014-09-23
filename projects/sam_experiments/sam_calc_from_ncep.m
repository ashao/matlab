infiles.z700mb = '/ltraid4/ncep/NCEPV1/pres/hgt.mon.mean.nc';
ncep.z700mb = nc_varget(infiles.z700mb,'hgt',[0 3 0 0],[-1 1 -1 -1]);
ncep.time = nc_varget(infiles.z700mb,'time')/24+datenum(1,1,1);
ncep.lon = nc_varget(infiles.z700mb,'lon');
ncep.lat = nc_varget(infiles.z700mb,'lat');
%%
% Calculate the Southern Annular Mode from 700mb surface as the first EOF
% between 1979-2000

[ncep.longrid ncep.latgrid] = meshgrid(ncep.lon,ncep.lat);
sidx = ncep.latgrid<-20 & ncep.latgrid > -90;
samtidx = ncep.time > datenum(1979,1,1) & ncep.time < datenum(2001,1,1);

sam.zhght = ncep.z700mb(samtidx,sidx);
[ntime npts] = size(sam.zhght);

%% Filter out the seasonal cycle
sam.filt.zhght = detrend(sam.zhght')';

fprintf('Filtering out seasonal cycle\n')
for i=1:npts
    for mon=1:12
       tidx = mon:12:ntime;
       sam.filt.zhght(tidx,i)=sam.filt.zhght(tidx,i)- ...
           mean(sam.filt.zhght(tidx,i));       
    end    
end

wts = repmat(ncep.latgrid(sidx)',[ntime 1]);
sam.filt.zhght=sam.filt.zhght.*wts;
%%
[U S V] = svds(double(sam.filt.zhght),1);
sam.pc = U*S;;
sam.eof = nan(size(ncep.longrid));
sam.eof(sidx) = V;

%% Filter all the data and project onto the mode
sozhght=ncep.z700mb(:,sidx);
[ntime npts] = size(sozhght);
for i=1:npts
    for mon=1:12
        tidx = mon:12:ntime;
        sozhght(tidx,i)=sozhght(tidx,i)-mean(sozhght(tidx,i));
    end
end
wts = repmat(ncep.latgrid(sidx)',[ntime 1]);
sozhght=sozhght.*wts;
%%
sam.pc = -sozhght*V;
colormap(othercolor('BuDRd_12'))
m_proj('stereographic','lon',0,'lat',-90,'radius',70)
subplot(3,1,[1 2])
m_contourf(ncep.lon,ncep.lat,sam.eof)
m_coast('line','Color','black','linewidth',2);
m_grid;
caxis([-1 1]*2.5e-2)
colorbar

subplot(3,1,3)
area(ncep.time,smooth(sam.pc./std(sam.pc(samtidx)),3,'mean'))
xlim([datenum(1948,1,1),datenum(2008,1,1)])
datetick('KeepLimits')
ylim([-2 2])
grid on
