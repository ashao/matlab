infiles.control = '/ltraid2/darr/GOLDruns_hyak/50yr-nz63/ocean_month.nc';
% infiles.saltflux50 = '/ltraid3/ashao/GOLD/gold_aabw_diagnose/restoring/temp.s1.nc';
% infiles.saltflux100 = '/ltraid3/ashao/GOLD/gold_aabw_diagnose/restoring/temp.s0.5.nc';
% infiles.saltbuoyflux50 = '/ltraid3/ashao/GOLD/gold_aabw_diagnose/restoring/temp.newgold.nc';
infiles.newgold = '/ltraid3/ashao/GOLD/gold_aabw_diagnose/newgold/temp.newgold.00.40.nc';
infiles.himalbedo = '/ltraid3/ashao/GOLD/gold_aabw_diagnose/newgold/albedo.himvals.00.40.nc';
infiles.ic = '/ltraid3/ashao/GOLD/gold_aabw_diagnose/newgold/GOLD_IC.nc'

% infiles.ic = '/home/ashao/uw-apl/data/gold/nz63/GOLD_IC.2010.11.15.nc';
%% Calculate global average
nyears =40;
globavgtmp.control =zeros(nyears,1);
globavgtmp.saltflux100=zeros(nyears,1);
globavgtmp.saltflux50 =zeros(nyears,1);
globavgtmp.saltbuoyflux50 =zeros(nyears,1);

load metrics
Ah = metrics.Ah.data.*metrics.wet.data;
Ahgrid = [];
Ahgrid(1,:,:)=Ah;
Ahgrid=repmat(Ahgrid,[63,1,1]);

for t=1:nyears
    fprintf('\n Year %d: ',t)
    fprintf('Control...')
    % average control run
    temp = nc_varget(infiles.control,'temp',[(t-1)*12 0 0 0],[12 -1 -1 -1]);
    h = nc_varget(infiles.control,'h',[(t-1)*12 0 0 0],[12 -1 -1 -1]);
    
    temp = squeeze(mean(temp));
    h = squeeze(mean(h));
    wts = Ahgrid.*h;
    wts = wts./sum(makevec(wts));
    globavgtmp.control(t) = sum(makevec(wts.*temp));
    
    
    fprintf('Salt Flux 100...')
    % Restoring constant runs
    temp = nc_varget(infiles.saltflux100,'temp',[t-1 0 0 0],[1 -1 -1 -1]);
    h = nc_varget(infiles.saltflux100,'h',[t-1 0 0 0],[1 -1 -1 -1]);
    wts = Ahgrid.*h;
    wts = wts./sum(makevec(wts));
    globavgtmp.saltflux100(t) = nansum(makevec(wts.*temp));
    
    fprintf('Salt Flux 50...')
    temp = nc_varget(infiles.saltflux50,'temp',[(t-1) 0 0 0],[1 -1 -1 -1]);
    h = nc_varget(infiles.saltflux50,'h',[(t-1) 0 0 0],[1 -1 -1 -1]);
    wts = Ahgrid.*h;
    wts = wts./sum(makevec(wts));
    globavgtmp.saltflux50(t) = nansum(makevec(wts.*temp));
    
    fprintf('Salt Buoy Flux 50...')
    temp = nc_varget(infiles.newgold,'temp',[(t-1) 0 0 0],[1 -1 -1 -1]);
    h = nc_varget(infiles.newgold,'h',[(t-1) 0 0 0],[1 -1 -1 -1]);
    wts = Ahgrid.*h;
    wts = wts./sum(makevec(wts));
    globavgtmp.newgold(t) = nansum(makevec(wts.*temp));
end

save restoringexperiments.mat globavgtmp
%% Extract 1 year of data

year = 10;
% temp.control = nc_varget(infiles.control,'temp',[(year-1)*12 0 0 0],[12 -1 -1 -1]);
% h.control = nc_varget(infiles.control,'h',[(year-1)*12 0 0 0],[12 -1 -1 -1]);
% temp.control = squeeze(mean(temp.control));
% h.control = squeeze(mean(h.control));

% Restoring constant runs
% temp.saltflux100 = nc_varget(infiles.saltflux100,'temp',[year-1 0 0 0],[1 -1 -1 -1]);
% h.saltflux100 = nc_varget(infiles.saltflux100,'h',[year-1 0 0 0],[1 -1 -1 -1]);

% temp.saltflux50 = nc_varget(infiles.saltflux50,'temp',[year-1 0 0 0],[1 -1 -1 -1]);
% h.saltflux50 = nc_varget(infiles.saltflux50,'h',[year-1 0 0 0],[1 -1 -1 -1]);

% temp.saltbuoyflux50 = nc_varget(infiles.saltbuoyflux50,'temp',[year-1 0 0 0],[1 -1 -1 -1]);
% h.saltbuoyflux50 = nc_varget(infiles.saltbuoyflux50,'h',[year-1 0 0 0],[1 -1 -1 -1]);
% temp.newgold= nc_varget(infiles.newgold,'temp',[year-1 0 0 0],[1 -1 -1 -1]);
% h.newgold = nc_varget(infiles.newgold,'h',[year-1 0 0 0],[1 -1 -1 -1]);
% 
% temp.himalbedo= nc_varget(infiles.himalbedo,'temp',[year-1 0 0 0],[1 -1 -1 -1]);
% h.himalbedo = nc_varget(infiles.himalbedo,'h',[year-1 0 0 0],[1 -1 -1 -1]);
temp.ic = nc_varget(infiles.ic,'Temp');
h.ic =nc_varget(infiles.ic,'h');
salt.ic  =nc_varget(infiles.ic,'Salt');
depth.ic = cumsum(h.ic);
%%
ic.temp = nc_varget(infiles.control,'temp');
ic.h =nc_varget(infiles.control,'h');
ic.depth = cumsum(ic.h);

%%
fields = fieldnames(temp);
for i=1:length(fields)
    depth.(fields{i})=cumsum(h.(fields{i}));
    
end
%% Plot Southern Ocean lines
titles = { 'GOLD (2011) Fluxconst=0.5',...
    'GOLD (2013) Fluxconst =1',...
    'HIM Albedo Values','Initial Condition'};
wocelabels = {'S1 68W','A23 30W','S2 0E','I6 30E','I8 270W','I9 245W', ...
    'S3 215W','P11 205W','P14 190W','P15 170W','P16 150W','P17 135W', ...
    'P18 105W','P19 90W'};
wocenames = {'S1','A23','S2','I6','I8','I9','S3','P11','P14','P15', ...
    'P16','P17','P18','P19'};
wocelons = [-68.5 -30.5 0.5 30.5 -270.5 -245.5 -215.5 -205.5 -190.5 ...
    -170.5 -150.5 -135.5 -105.5 -90.5];

%%
load metrics
for lineno = 1:length(wocelons);
    lineno
    s3idx = find(abs(metrics.lonh.data- wocelons(lineno))<0.2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    colormap(othercolor('BuDRd_12'))
    counter = 0;
    for i=1:length(fields)
        counter = counter+1;
        subplot(2,2,counter)
        plotdepth = depth.(fields{i});
        plotdepth = squeeze(plotdepth(:,:,s3idx));
        plottemp = temp.(fields{i});
        plottemp= double(squeeze(plottemp(:,:,s3idx)));
        contourf(latgrid,plotdepth,plottemp,-2:.4:10);
        %     shading flat;
        xlim([-75. -45])
        set(gca,'ydir','reverse')
        cax=colorbar;
        ylabel(cax,'Temp')
        caxis([-2 2])
        title(titles{i})
        xlabel(wocelabels{lineno})
    end
    
    saveas(gcf,['/ltraid3/ashao/GOLD/gold_aabw_diagnose/figs/' wocenames{lineno} '.eps'],'epsc');
end
% contourf(
%%
year = 40;
infiles.default='/ltraid3/ashao/GOLD/default/default.newgold.00.40.nc';
temp.default= nc_varget(infiles.default,'temp',[year-1 0 0 0],[1 -1 -1 -1]);
salt.default= nc_varget(infiles.default,'salt',[year-1 0 0 0],[1 -1 -1 -1]);
h.default = nc_varget(infiles.default,'h',[year-1 0 0 0],[1 -1 -1 -1]);
depth.default=cumsum(h.default);
%%
templotcmds = ['plotdepth = squeeze(plotdepth(:,:,s3idx));' ...
    'plottemp= double(squeeze(plottemp(:,:,s3idx)));' ...    
    'contourf(latgrid,plotdepth,plottemp,-2:.4:10);' ...        
    'xlim([-75. -45]);ylim([0 5000]);' ...
    'set(gca,''ydir'',''reverse'');' ...
    'cax=colorbar;' ... 
    'ylabel(cax,''Temp'');' ...
    'caxis([-2 2]);' ...    
    'xlabel(wocelabels{lineno});'];

saltplotcmds = ['plotdepth = squeeze(plotdepth(:,:,s3idx));' ...
    'plotsalt= double(squeeze(plotsalt(:,:,s3idx)));' ...    
    'contourf(latgrid,plotdepth,plotsalt,33.5:0.05:35);' ...        
    'xlim([-75. -45]);ylim([0 5000]);' ...
    'set(gca,''ydir'',''reverse'');' ...
    'cax=colorbar;' ... 
    'ylabel(cax,''Salt'');' ...
    'caxis([34 35]);' ...    
    'xlabel(wocelabels{lineno});'];
 
% densplotcmds= ['plotdepth = squeeze(plotdepth(:,:,s3idx));' ...
%     'plotdens= double(squeeze(plotsalt(:,:,s3idx)));' ...    
%     'contourf(latgrid,plotdepth,plotsalt,34:0.02:35);' ...    
%     'xlim([-75. -45])' ...
%     'set(gca,''ydir'',''reverse'')' ...
%     'cax=colorbar;' ... 
%     'ylabel(cax,''Density'')' ...
%     'caxis([1015:0.5:1030])' ...    
%     'title(wocelabels{lineno})'];
%% Initial Condition and GOLD run Compariosn
load metrics
for lineno = 2:length(wocelons);

    figure('Visible','Off')
    s3idx = find(abs(metrics.lonh.data- wocelons(lineno))<0.2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    colormap(othercolor('BuDRd_12'))
    counter = 0;

    % Temperature
    subplot(2,2,1)
    plotdepth = depth.ic;
    plottemp = temp.ic;
    eval(templotcmds);
        xlabel('Initial Conditions')
    
    subplot(2,2,3)
    plotdepth = depth.default;
    plottemp = temp.default;
    eval(templotcmds);
        xlabel('Default GOLD (40yr)')
    
    % Salt
    subplot(2,2,2)
    plotdepth = depth.ic;
    plotsalt = salt.ic;
    eval(saltplotcmds);
    xlabel('Initial Conditions')
    
    subplot(2,2,4)
    plotdepth = depth.default;
    plotsalt= salt.default;
    eval(saltplotcmds);
    xlabel('Default GOLD (40yr)')
    
    
%     saveas(gcf,['/ltraid3/ashao/GOLD/gold_aabw_diagnose/figs/default/' wocenames{lineno} '.eps'],'epsc');
    
end

%% ESM2G and GOLD comparison
infiles.esm2g.salt= '/ltraid4/cmip5/noaa/gfdl/so_Omon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc';
infiles.esm2g.temp = '/ltraid4/cmip5/noaa/gfdl/thetao_Omon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc';
temp.esm2g = squeeze(mean(nc_varget(infiles.esm2g.temp,'thetao')));
salt.esm2g = squeeze(mean(nc_varget(infiles.esm2g.salt,'so')));
depth.esm2g= nc_varget(infiles.esm2g.salt,'lev');
depth.esm2g=repmat(depth.esm2g,[1 210 360]);
%%
load metrics
for lineno = 1:length(wocelons);

    figure
    s3idx = find(abs(metrics.lonh.data- wocelons(lineno))<0.2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    colormap(flipud(othercolor('PuOr10')))
    counter = 0;

    % Temperature
    subplot(3,2,1)
    latgrid = repmat(metrics.lath.data',[50 1]);
    plotdepth = depth.esm2g;
    plottemp = temp.esm2g-273.15;
    eval(templotcmds);
        title('ESM2G')
    
    subplot(3,2,3)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = depth.default;
    plottemp = temp.default;
    eval(templotcmds);
        title('Default GOLD (40yr)')
    
    % Salt
    subplot(3,2,2)
    latgrid = repmat(metrics.lath.data',[50 1]);
    plotdepth = depth.esm2g;
    plotsalt = salt.esm2g;
    eval(saltplotcmds);
    title('ESM2G')
    
    subplot(3,2,4)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = depth.default;
    plotsalt= salt.default;
    eval(saltplotcmds);
    title('Default GOLD (40yr)')
    
        subplot(3,2,5)
    plotdepth = depth.ic;
    plottemp = temp.ic;
    eval(templotcmds);
        title('Initial Conditions')
        
                subplot(3,2,6)
    plotdepth = depth.ic;
    plotsalt = salt.ic;
    eval(saltplotcmds);
        title('Initial Conditions')
    
    
    saveas(gcf,...
        sprintf('/ltraid3/ashao/GOLD/gold_aabw_diagnose/figs/esm2g_gold/%02d_%s.eps',lineno,wocenames{lineno}),'epsc');
    
end
%% Kd plot

kd = nc_varget('/ltraid3/ashao/temp.salt.Kd.00.05.nc','Kd',[4 0 0 0],[1 -1 -1 -1]);
depth.kd = cumsum(nc_varget('/ltraid3/ashao/temp.salt.Kd.00.05.nc','h',[4 0 0 0],[1 -1 -1 -1]));
temp.kd = nc_varget('/ltraid3/ashao/temp.salt.Kd.00.05.nc','temp',[4 0 0 0],[1 -1 -1 -1]);

for lineno = 1:length(wocelons);
    figure
    subplot(2,1,1);
    s3idx = find(abs(metrics.lonh.data- wocelons(lineno))<0.2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    colormap(othercolor('BuDRd_12'))
    contourf(latgrid,depth.kd(:,:,s3idx),log10(kd(:,:,s3idx)),-5:0.25:-1)
    xlim([-75 -45]);ylim([0 1000])
    caxis([-5 -2])
    cax=colorbar;
    ylabel(cax,'log10(K_d)')
    set(gca,'ydir','reverse')    
    xlabel(wocelabels{lineno})
    
    
    subplot(2,1,2);
    s3idx = find(abs(metrics.lonh.data- wocelons(lineno))<0.2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    colormap(othercolor('BuDRd_12'))
    contourf(latgrid,depth.kd(:,:,s3idx),temp.kd(:,:,s3idx),-2:0.25:2)
    xlim([-75 -45]);ylim([0 1000]);
    cax=colorbar;
    ylabel(cax,'Temp')
    set(gca,'ydir','reverse')    
    xlabel(wocelabels{lineno})
    
    saveas(gcf,...
        sprintf('/ltraid3/ashao/GOLD/gold_aabw_diagnose/figs/%02d_kd.temp.%s.eps',lineno,wocenames{lineno}),'epsc');
end