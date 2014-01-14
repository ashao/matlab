goldpath = '/ltraid3/ashao/gold_aabw_diagnose/albedo/mid-albedo/';
% goldpath = '/ltraid2/darr/GOLDruns_hyak/200yr-nz63/';
goldfile.albedo.z = [goldpath 'ocean_month_z.nc'];
goldfile.albedo.isopycnal = [goldpath 'ocean_month.nc'];
goldfile.albedo.ice= [goldpath 'ice_month.nc'];
goldfile.albedo.static = [goldpath 'ocean_geometry.nc'];
goldfile.albedo.vertical = [goldpath 'Vertical_coordinate.nc'];

gold.albedo.time = nc_varget(goldfile.albedo.isopycnal,'time');
gold.albedo.layer = nc_varget(goldfile.albedo.vertical,'Layer');
gold.albedo.geolon = nc_varget(goldfile.albedo.static,'geolon');
gold.albedo.geolat = nc_varget(goldfile.albedo.static,'geolat');
gold.albdeo.wet = logical(nc_varget(goldfile.albedo.static,'wet'));
%%
% 
% gold.temp = mean(nc_varget(goldfile.isopycnal,'temp',[120 0 0 0],[12 -1 -1 -1]));
% gold.salt = mean(nc_varget(goldfile.isopycnal,'salt',[120 0 0 0],[12 -1 -1 -1]));

%  gold.temp = nc_varget(goldfile.isopycnal,'temp');
%  gold.salt = nc_varget(goldfile.isopycnal,'salt');
%%
% % gold.depth = mean(nc_varget(goldfile.isopycnal,'h',[120 0 0 0],[12 -1 -1 -1]));
% gold.depth = mean(nc_varget(goldfile.isopycnal,'h');
% gold.depth = cumsum(gold.depth,2);
% %
%%
load metrics
dim.albedo.ntime = length(gold.albedo.time);
dim.albedo.nlayer =length(gold.albedo.layer);
[dim.nlon, dim.nlat]=size(gold.albedo.geolon);
outpath= [goldpath 'figs'];
mkdir(outpath)

%% Adjusted Albedo

for yidx = 1:(dim.albedo.ntime-1)
    
    gold.albedo.temp = nc_varget(goldfile.albedo.isopycnal,'temp',[yidx 0 0 0],[1 -1 -1 -1]);
    gold.albedo.salt = nc_varget(goldfile.albedo.isopycnal,'salt',[yidx 0 0 0],[1 -1 -1 -1]);
    gold.albedo.depth = cumsum(nc_varget(goldfile.albedo.isopycnal,'h',...
        [yidx 0 0 0],[1 -1 -1 -1]));    
    
    lonidx = find(metrics.lonh.data== -240.5);
    latarray = repmat(gold.albedo.geolat(:,lonidx)',[63 1]);
    subplot(1,2,1)
    colormap(othercolor('BuDRd_12'))
    contourf(squeeze(latarray),squeeze(gold.albedo.depth(:,:,lonidx)),...
        double(squeeze(gold.albedo.temp(:,:,lonidx))),-2:0.2:8)
    set(gca,'ydir','reverse')
    xlim([-80 -30])
%     ylim([0 2000])
    caxis([-2 8])
    colorbar
%     shading flat
    axis square
    title(sprintf('Temperature Model Year %02d',yidx))
    
    
    subplot(1,2,2)
    contourf(squeeze(latarray),squeeze(gold.albedo.depth(:,:,lonidx)),...
        double(squeeze(gold.albedo.salt(:,:,lonidx))),33:0.05:34.8)
    set(gca,'ydir','reverse')
    colormap(othercolor('BuDRd_12'))
    xlim([-75 -30])
%     ylim([0 2000])
    caxis([34.1 34.8])
%     shading flat
colorbar
    axis square
    title(sprintf('Salinity Model Year %02d',yidx))
    drawnow;
    filename=sprintf('s3.salt.temp.y%02d.png',yidx);
    print(gcf,'-dpng',[outpath filename])
    
end
%% Control
goldpath = '/ltraid2/darr/GOLDruns_hyak/200yr-nz63/';
goldfile.control.z = [goldpath 'ocean_month_z.nc'];
goldfile.control.isopycnal = [goldpath 'ocean_month.nc'];
goldfile.control.ice= [goldpath 'ice_month.nc'];
goldfile.control.static = [goldpath 'ocean_geometry.nc'];
goldfile.control.vertical = [goldpath 'Vertical_coordinate.nc'];

gold.control.time = nc_varget(goldfile.control.isopycnal,'time');
gold.control.layer = nc_varget(goldfile.control.vertical,'Layer');
gold.control.geolon = nc_varget(goldfile.control.static,'geolon');
gold.control.geolat = nc_varget(goldfile.control.static,'geolat');
gold.control.wet = logical(nc_varget(goldfile.control.static,'wet'));
load metrics
dim.control.ntime = length(gold.control.time);
dim.control.nlayer =length(gold.control.layer);
[dim.nlon, dim.nlat]=size(gold.control.geolon);

%%
outpath = '/ltraid2/ashao/uw-apl/figs/gold_aabw_diagnose/';
plotcmds = [ 'xlim([-80 -30]);caxis([-2 8]);ylim([0 5250]);colorbar;axis(''square''); '];
clevels = -1.8:.5:7.7;
for yidx = 0:(dim.albedo.ntime)
%     yidx = 0;
    gold.control.temp = squeeze(mean(nc_varget(...
        goldfile.control.isopycnal,'temp',[yidx*12 0 0 0],[12 -1 -1 -1])));
    gold.control.depth = squeeze(mean(nc_varget(...
        goldfile.control.isopycnal,'h',...
        [yidx*12 0 0 0],[12 -1 -1 -1])));    
    gold.control.depth = cumsum(gold.control.depth);
    
    gold.albedo.temp = squeeze(nc_varget(...
        goldfile.albedo.isopycnal,'temp',[yidx 0 0 0],[1 -1 -1 -1]));
    gold.albedo.depth = squeeze(nc_varget(...
        goldfile.albedo.isopycnal,'h',...
        [yidx 0 0 0],[1 -1 -1 -1]));    
    gold.albedo.depth = cumsum(gold.albedo.depth);
    
    
    lonidx = find(metrics.lonh.data== -221.5);
    latarray = repmat(gold.control.geolat(:,lonidx)',[63 1]);
    
    
    subplot(1,2,1)
%     colormap(othercolor('BuDRd_12'))
    colormap('jet')
    [cs h]=contourf(squeeze(latarray),squeeze(gold.control.depth(:,:,lonidx)),...
        double(squeeze(gold.control.temp(:,:,lonidx))),clevels);
    clabel(cs,h);
    set(gca,'ydir','reverse')
    eval(plotcmds);
    title(sprintf('\theta Control Model Year %02d',yidx))
    
    
    subplot(1,2,2)
    colormap('jet')
%         colormap(othercolor('BuDRd_12'))
    [cs h]=contourf(squeeze(latarray),squeeze(gold.albedo.depth(:,:,lonidx)),...
        double(squeeze(gold.albedo.temp(:,:,lonidx))),clevels);
    clabel(cs,h);
    set(gca,'ydir','reverse')
    eval(plotcmds);
%     axis square
    title(sprintf('\theta Albedo Model Year %02d',yidx))    
    
    filename = [outpath sprintf('s03.ptemp.compare.y%02d.png',yidx)];
    print(gcf,'-dpng',filename);
%     pause
end
