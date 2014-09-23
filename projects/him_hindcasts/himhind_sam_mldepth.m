infile.slp = '/ltraid4/ncep/NCEPV1/pres/hgt.mon.mean.nc';
infile.mldepth = '/ltraid3/ashao/uw-apl/models/HIM/hindcast/hindcast.mldepth.sst.mat';

%%
ncep.lat = nc_varget(infile.slp,'lat');
latidx = find(ncep.lat<-20 & ncep.lat>-90);
start4d = [0 3 min(latidx)-1 0];
count4d = [-1 1 length(latidx) -1];
ncep.so.slp = nc_varget(infile.slp,'hgt',start4d,count4d);
ncep.so.lat = nc_varget(infile.slp,'lat',min(latidx-1),length(latidx));
ncep.so.lon = nc_varget(infile.slp,'lon');
%%
[longrid latgrid] = meshgrid(ncep.so.lon,ncep.so.lat);
ncep.time = nc_varget(infile.slp,'time')/24+datenum(1,1,1);
% Only get the times betwen 1979 and 2000 for the EOF anlaysis
ncep.sam.tidx = ncep.time > datenum(1979,1,1) & ncep.time < datenum(2001,1,1); 

%%
load(infile.mldepth)
%%
load metrics
him.so.sidx = metrics.geolat.data < -20 & logical(metrics.wet.data);
him.sam.tidx = ((1979-1948)*12+1):((2001-1948)*12)
%%
him.so.mldepth = him.global.h(:,him.so.sidx);
him.so.wts = metrics.Ah.data(him.so.sidx)';
him.so.wts = repmat(him.so.wts,[length(him.sam.tidx) 1])./sum(him.so.wts(:));

%%
[ntime nhimpts] = size(him.so.mldepth(him.sam.tidx,:));
% him.so.mldepth_filt = him.so.mldepth;

him.so.mldepth_filt=detrend(him.so.mldepth(him.sam.tidx,:)')';
for i=1:nhimpts
    if mod(i,500)==0
    fprintf('HIM %d/%d\n',i,nhimpts);
    end
%     him.so.sst_filt(:,i) = annual_harms(him.so.sst(:,i),him.time(him.so.tidx),6,0,365.25);
    
    for t=1:12
            tidx = t:12:ntime;
            him.so.mldepth_filt(tidx,i) = him.so.mldepth_filt(tidx,i)-mean(him.so.mldepth_filt(tidx,i));
    end
end
him.so.mldepth_filt = him.so.mldepth_filt.*him.so.wts;

%%

% ncep.so.slp_filt = ncep.so.slp;
ncep.so.slp_filt=detrend(ncep.so.slp(ncep.sam.tidx,:)')';
[ntime nnceppts] = size(ncep.so.slp_filt);
for i=1:nnceppts   
    for t=1:12        
        tidx = t:12:ntime;
        ncep.so.slp_filt(tidx,i) = ncep.so.slp_filt(tidx,i)-mean(ncep.so.slp_filt(tidx,i));
    end
end

wts = repmat(sqrt(cosd(latgrid(:)))',[sum(ncep.sam.tidx) 1]);
wts = wts./sum(wts(1,:));
ncep.so.slp_filt=ncep.so.slp_filt.*wts;
% ncep.so.slp_filt = detrend(ncep.so.slp_filt')';


%%
R_xy = double(ncep.so.slp_filt'*him.so.mldepth_filt./(ntime-1));

[U S V] = svds(R_xy);

% save CCPA_so.mat U S V so -v7.3

%%

% testeof = V*S;
colormap(othercolor('BuDRd_12'))
mode = 1;
clf
subplot(3,4,[1 2 5 6])
himfillmat = nan(210,360);
himfillmat(him.so.sidx)=V(:,mode);

m_proj('Stereographic','lat',-90,'lon',0,'radius',70)
levels = -0.1:0.0025:0.1;
m_contourf(metrics.geolon.data,metrics.geolat.data,himfillmat,levels,'LineColor','None')
% shading flat
m_coast('line','LineWidth',2,'Color','Black');
m_grid;
caxis([-1 1]*2.5e-2)
% cbarf(himfillmat,-5e-2:0.005:5e-2);
colorbar
title('Mixed Layer Depth Mode 1')

subplot(3,4,[3 4 7 8])
% ncepfillmat = nan(size(ncep.so.slp));
ncepfillmat=reshape(U(:,mode),size(latgrid));

m_contourf(ncep.so.lon,ncep.so.lat,ncepfillmat,levels,'LineColor','None')
shading flat
m_coast('line','LineWidth',2,'Color','Black');
caxis([-1 1]*2.5e-2)
m_grid;
% cbarf(ncepfillmat,-0.05:0.005:0.05);
colorbar;
title('700mb Z-Height Coupled Mode 1')


temparray = ncep.so.slp(1:end,:);
% temparray = detrend(temparray')';
[ntime npts ] =size(temparray);
for i= 1:npts
    for t=1:12        
        tidx = t:12:ntime;
        temparray(tidx,i) = temparray(tidx,i)-mean(temparray(tidx,i));
    end
end
SL = temparray*U;

temparray = him.so.mldepth(1:end,:);
% temparray = detrend(temparray')';
[ntime npts ] =size(temparray);
for i= 1:npts
    for t=1:12        
        tidx = t:12:ntime;
        temparray(tidx,i) = temparray(tidx,i)-mean(temparray(tidx,i));
    end
end
SR = temparray*V;
subplot(3,4,[9 10])
pcR = smooth(SR(:,mode),6,'mean');
pcR_norm = pcR./std(pcR(:,mode)); hold on;
plot(him.time/365+1948,pcR_norm,'LineWidth',1.5,'Color','Black')
plot(him.time/365+1948,polyval(polyfit(him.time/365+1948,pcR_norm,1),him.time/365+1948))
% datetick('KeepTicks')
ylim([-3 3])
xlim([1940 2020])
grid on

subplot(3,4,[11 12])
pcL = smooth(SL(:,mode),6,'mean');
pcL_norm = pcL./std(pcL);
hold on
plot(ncep.time ,pcL_norm,'LineWidth',1.5,'Color','Black')
plot(ncep.time,polyval(polyfit(ncep.time,pcL_norm,1),ncep.time))
ylim([-3 3])
xlim([datenum(1940,1,1) datenum(2020,1,1)])
datetick('KeepLimits')
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