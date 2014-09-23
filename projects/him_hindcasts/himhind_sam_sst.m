infile.slp = '/ltraid4/ncep/hgt.mon.mean.nc';
infile.sst = '/ltraid3/ashao/uw-apl/models/HIM/hindcast/hindcast.mldepth.sst.mat';

years = 1979:2007;
nyears=length(years);
ncep.slp = nc_varget(infile.slp,'hgt',[0 3 0 0],[nyears*12 1 inf inf]);
ncep.lat = nc_varget(infile.slp,'lat');
ncep.lon = nc_varget(infile.slp,'lon');
[longrid latgrid] = meshgrid(ncep.lon,ncep.lat);
ncep.so.idx = latgrid <= -20 & latgrid>-90;
ncep.so.slp = ncep.slp(:,ncep.so.idx);
% ncep.so.slp = ncep.so.slp(1:348,:);
ncep.time = nc_varget(infile.slp,'time',0,nyears*12)/24;

%%
load(infile.sst)
%%
load metrics
him.so.sidx = metrics.geolat.data < -20 & logical(metrics.wet.data);
him.so.tidx = 373:720;
%%
him.so.sst = him.global.sst(him.so.tidx,him.so.sidx);
him.so.wts = metrics.Ah.data(him.so.sidx)';
him.so.wts = repmat(him.so.wts,[348 1])./sum(him.so.wts(:));

%%
[ntime nhimpts] = size(him.so.sst);

% him.so.sst_filt=detrend(him.so.sst_filt')';
for i=1:nhimpts
    if mod(i,500)==0
    fprintf('HIM %d/%d\n',i,nhimpts);
    end
%     him.so.sst_filt(:,i) = annual_harms(him.so.sst(:,i),him.time(him.so.tidx),6,0,365.25);
    
    for t=1:12
            tidx = t:12:ntime;
            him.so.sst_filt(tidx,i) = him.so.sst_filt(tidx,i)-mean(him.so.sst_filt(tidx,i));
    end
end
him.so.sst_filt = him.so.sst_filt.*him.so.wts;

%%
[ntime nnceppts] = size(ncep.so.slp);
% ncep.so.slp_filt=detrend(ncep.so.slp')';
for i=1:nnceppts   
    for t=1:12        
        tidx = t:12:ntime;
        ncep.so.slp_filt(tidx,i) = ncep.so.slp_filt(tidx,i)-mean(ncep.so.slp_filt(tidx,i));
    end
end
ncep.so.slp_filt=ncep.so.slp_filt./repmat(sqrt(cosd(latgrid(ncep.so.idx)))',[348 1]);
% ncep.so.slp_filt = detrend(ncep.so.slp_filt')';


%%
R_xy = double(ncep.so.slp_filt'*him.so.sst_filt./(ntime-1));

[U S V] = svds(R_xy);

% save CCPA_so.mat U S V so -v7.3

%%

% testeof = V*S;
colormap(othercolor('BuDRd_12'))
mode = 1;

subplot(3,4,[1 2 5 6])
himfillmat = nan(210,360);
himfillmat(him.so.sidx)=V(:,mode);

m_proj('Stereographic','lat',-90,'lon',0,'radius',70)
levels = -0.1:0.0025:0.1;
m_contourf(metrics.geolon.data,metrics.geolat.data,himfillmat,levels,'LineColor','None')
% shading flat
m_coast('line','LineWidth',2,'Color','Black');
m_grid;
caxis([-1 1]*2e-2)
cbarf(himfillmat,-2e-2:0.005:2e-2);
title('SST Coupled Mode 1')

subplot(3,4,[3 4 7 8])
ncepfillmat = nan(length(ncep.lat),length(ncep.lon));
ncepfillmat(ncep.so.idx)=U(:,mode);

m_contourf(ncep.lon,ncep.lat,ncepfillmat,levels,'LineColor','None')
shading flat
m_coast('line','LineWidth',2,'Color','Black');
caxis([-2 2]*1e-2)
m_grid;
cbarf(ncepfillmat,-0.02:0.005:0.02);
title('700mb Z-Height Coupled Mode 1')

SL = ncep.so.slp_filt*U;
SR = him.so.sst_filt*V;

subplot(3,4,[9 10])
pcR = SR(:,mode);
pcR = pcR./std(pcR(:,mode));
plot(ncep.time + datenum(1800,1,1),pcR,'LineWidth',1.5,'Color','Black')
datetick('KeepTicks')
ylim([-3 3])
grid on

subplot(3,4,[11 12])
pcL = smooth(SL(:,mode),6,'mean');
pcL = pcL./std(pcL);
plot(ncep.time + datenum(1800,1,1),pcL,'LineWidth',1.5,'Color','Black')
ylim([-3 3])
datetick('KeepTicks')
grid on

% cbarf(

%% Compare SST to actual SAM
samurl= 'http://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/aao/monthly.aao.index.b79.current.ascii';
urlwrite(samurl,'temp.txt');
[year mon index] = csvimport('temp.txt','delimiter',' ','column',[1 2 3],'noheader',true);

%%
subplot(2,1,1)
pcR = SR(:,mode);
pcR = pcR./std(pcR(:,mode));
plot(ncep.time + datenum(1800,1,1),pcR,'LineWidth',1.5,'Color','Black')
datetick('KeepTicks')
ylim([-3 3])
grid on

subplot(2,1,2)
pcL = smooth(SL(:,mode),8,'mean');
% pcL = SL(:,mode);
pcL = pcL./std(pcL);
plot(ncep.time + datenum(1800,1,1),pcL,'LineWidth',1.5,'Color','Black')
ylim([-3 3])
datetick('KeepTicks')
grid on

max(xcorr(pcL,pcR,'biased'))