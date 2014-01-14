function [ mean_age ] = off_calc_ttd_moments( ncfile, hfile, layer )
%% function [ mean_age ttd_width ] = calc_ttd_moments( ncfile, hfile )
% Calculate the first two moments of the TTD
% Extrapolate the tail by fitting an inverse gaussian
%


minh=1e-3;
load metrics.mat
time=single(nc_varget(ncfile,'Time'));
ntime=length(time);
lath=nc_varget(ncfile,'lath');
nlat=length(lath);
lonh=nc_varget(ncfile,'lonh');
nlon=length(lonh);
mean_age=zeros(nlat,nlon,'single');
ttd_width=zeros(nlat,nlon,'single');
warning off all
tailidx=1;
nyear=floor(ntime/12);
years=1:nyear;
start_4d=[1 layer 1 1]-1;
count_4d=[Inf 1 Inf Inf];

fprintf('Data loading...')
ttd=single(nc_varget(ncfile,'TTD',start_4d,count_4d));
% size(ttd)
fprintf('Loaded\n')
h_layer=squeeze(min(single(nc_varget(hfile,'h',start_4d,count_4d))));
timegrid(:,1)=time;
timegrid=repmat(timegrid,[1 nlon]);

for lat=1:nlat
   
    subttd=squeeze(ttd(:,lat,:));
%     whos
    normttd=subttd./repmat(trapz(time,subttd),[ntime 1]);
    mean_age(lat,:)=trapz(time,timegrid.*normttd);       
    
end

% h_layer=squeeze(min(single(ncread(hfile,'h',start_4d,count_4d))));
h_layer(h_layer<=1.0e-9)=0;
thin=~logical(h_layer);
mean_age(~metrics.wet.data)=NaN;
mean_age(thin)=NaN;

% depth=cumsum(nc_varget(hfile,'h'),2);