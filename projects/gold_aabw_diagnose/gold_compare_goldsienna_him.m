infiles.gold= '/ltraid4/ashao/gold/HINDCAST-600yr/ocean_month.1948.nc';
infiles.ic = '/ltraid4/ashao/save/GOLD/gold_aabw_diagnose/newgold/GOLD_IC.nc';
infiles.esm2g.temp = '/ltraid4/cmip5/noaa/gfdl_esm2g_r1i1p1/incoming/thetao_Omon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc';
infiles.esm2g.salt = '/ltraid4/cmip5/noaa/gfdl_esm2g_r1i1p1/incoming/so_Omon_GFDL-ESM2G_historical_r1i1p1_200101-200512.nc';
infiles.him = '/ltraid4/ashao/HIM/hyak_store/NORMALYEAR/month/ocean_month.nc'

wocelabels = {'S1 68W','A23 30W','S2 0E','I6 30E','I8 270W','I9 245W', ...
    'S3 215W','P11 205W','P14 190W','P15 170W','P16 150W','P17 135W', ...
    'P18 105W','P19 90W'};
wocenames = {'S1','A23','S2','I6','I8','I9','S3','P11','P14','P15', ...
    'P16','P17','P18','P19'};
wocelons = [-68.5 -30.5 0.5 30.5 -270.5 -245.5 -215.5 -205.5 -190.5 ...
    -170.5 -150.5 -135.5 -105.5 -90.5];

% Load model output

%% GOLD (Default)
gold.temp = squeeze(mean(nc_varget(infiles.gold,'temp')));
gold.salt = squeeze(mean(nc_varget(infiles.gold,'salt')));
gold.depth = cumsum(squeeze(mean( ...
    nc_varget(infiles.gold,'h'))));

%% Initial Conditions
ic.temp = nc_varget(infiles.ic,'Temp');
ic.salt = nc_varget(infiles.ic,'Salt');
ic.depth = cumsum( nc_varget(infiles.ic,'h'));
%% ESM2G

esm2g.temp = squeeze(mean(nc_varget(infiles.esm2g.temp,'thetao')));
esm2g.salt = squeeze(mean(nc_varget(infiles.esm2g.salt,'so')));
esm2g.depth = nc_varget(infiles.esm2g.temp,'lev' );
esm2g.depth=repmat(esm2g.depth,[1 210 360]);

%% HIM
start4d = [0 0 0 0];
count4d = [60 inf inf inf];
him.temp = squeeze(mean(nc_varget(infiles.him,'temp',start4d,count4d)));
him.salt = squeeze(mean(nc_varget(infiles.him,'salt',start4d,count4d)));
him.depth = cumsum(squeeze(mean( ...
    nc_varget(infiles.him,'h',start4d,count4d))));

%% Set plot commands



templotcmds = ['plotdepth = squeeze(plotdepth(:,:,s3idx));' ...
    'plottemp= double(squeeze(plottemp(:,:,s3idx)));' ...
    'contourf(latgrid,plotdepth,plottemp,-2:.4:10);' ...
    'xlim([-75. -45]);ylim([0 5000]);' ...
    'set(gca,''ydir'',''reverse'');' ...
    'cax=colorbar;' ...
    'ylabel(cax,''Temp'');' ...
    'caxis([-2 2]);' ...    
    'colormap(flipud(othercolor(''PuOr10'')))'];

saltplotcmds = ['plotdepth = squeeze(plotdepth(:,:,s3idx));' ...
    'plotsalt= double(squeeze(plotsalt(:,:,s3idx)));' ...
    'contourf(latgrid,plotdepth,plotsalt,33.5:0.05:35);' ...
    'xlim([-75. -45]);ylim([0 5000]);' ...
    'set(gca,''ydir'',''reverse'');' ...
    'cax=colorbar;' ...
    'ylabel(cax,''Salt'');' ...
    'caxis([34 35]);' ...
    'xlabel(wocelabels{lineno});'];
%%

load metrics


for lineno = 1:length(wocelons);
    
    figure
    s3idx = find(abs(metrics.lonh.data- wocelons(lineno))<0.2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    
    counter = 0;
    
    % Temperature
    subplot(4,2,1)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = ic.depth;
    plottemp = ic.temp;
    eval(templotcmds);
    title('Initial Conditions')
    
    subplot(4,2,7)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = gold.depth;
    plottemp = gold.temp;
    eval(templotcmds);
    title('Default GOLD (600yr)')
    
    subplot(4,2,5)
    latgrid = repmat(metrics.lath.data',[49 1]);
    plotdepth = him.depth;
    plottemp = him.temp;
    eval(templotcmds);
    title('HIM (500yr)')
    
    subplot(4,2,3)
    latgrid = repmat(metrics.lath.data',[50 1]);
    plotdepth = esm2g.depth;
    plottemp = esm2g.temp-273.15;
    eval(templotcmds);
    title('ESM2G')
    xlabel(wocelabels{lineno});
    
    % Salt
    subplot(4,2,2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = ic.depth;
    plotsalt = ic.salt;
    eval(saltplotcmds);
    title('Initial Conditions')
    
    subplot(4,2,8)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = gold.depth;
    plotsalt= gold.salt;
    eval(saltplotcmds);
    title('Default GOLD (600 yr)')   
    
    subplot(4,2,6)
    latgrid = repmat(metrics.lath.data',[49 1]);
    plotdepth = him.depth;
    plotsalt = him.salt;
    eval(saltplotcmds);
    title('HIM (500 yr)')
    
    subplot(4,2,4)
    latgrid = repmat(metrics.lath.data',[50 1]);
    plotdepth = esm2g.depth;
    plotsalt = esm2g.salt;
    eval(saltplotcmds);
    title('ESM2G')
    xlabel(wocelabels{lineno});
    
%     pause
    saveas(gcf,...
        sprintf('/ltraid3/ashao/uw-apl/figs/gold_aabw_diagnose/esm2g_gold_him/%02d_%s.esm2g.gold.him.eps',lineno,wocenames{lineno}),'epsc');
    
end