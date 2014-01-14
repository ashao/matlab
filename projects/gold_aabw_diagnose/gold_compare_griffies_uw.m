infiles.default= '/ltraid4/ashao/gold/default/default.newgold.00.40.nc';
infiles.ic = '/ltraid4/ashao/save/GOLD/gold_aabw_diagnose/newgold/GOLD_IC.nc';
infiles.griffies = '/ltraid4/ashao/from_GFDL/GOLD_for_Andrew/ocean_annual.1988-2007.ann.nc';

wocelabels = {'S1 68W','A23 30W','S2 0E','I6 30E','I8 270W','I9 245W', ...
    'S3 215W','P11 205W','P14 190W','P15 170W','P16 150W','P17 135W', ...
    'P18 105W','P19 90W'};
wocenames = {'S1','A23','S2','I6','I8','I9','S3','P11','P14','P15', ...
    'P16','P17','P18','P19'};
wocelons = [-68.5 -30.5 0.5 30.5 -270.5 -245.5 -215.5 -205.5 -190.5 ...
    -170.5 -150.5 -135.5 -105.5 -90.5];

%% Load model output
default.temp = nc_varget(infiles.default,'temp',[39 0 0 0],[1 inf inf inf]);
default.salt = nc_varget(infiles.default,'salt',[39 0 0 0],[1 inf inf inf]);
default.depth = cumsum( ...
    nc_varget(infiles.default,'h',[39 0 0 0],[1 inf inf inf]));

%%
ic.temp = nc_varget(infiles.ic,'Temp');
ic.salt = nc_varget(infiles.ic,'Salt');
ic.depth = cumsum( nc_varget(infiles.ic,'h'));
%%

griffies.temp = nc_varget(infiles.griffies,'temp');
griffies.salt = nc_varget(infiles.griffies,'salt');
griffies.depth = cumsum( nc_varget(infiles.griffies,'h' ));
%% Set plot commands

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
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = ic.depth;
    plottemp = ic.temp;
    eval(templotcmds);
    title('Initial Conditions')
    
    subplot(3,2,3)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = default.depth;
    plottemp = default.temp;
    eval(templotcmds);
    title('Default GOLD (40yr)')
    
    subplot(3,2,5)
    plotdepth = griffies.depth;
    plottemp = griffies.temp;
    eval(templotcmds);
    title('From Griffies CORE-2')
    
    % Salt
    subplot(3,2,2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = ic.depth;
    plotsalt = ic.salt;
    eval(saltplotcmds);
    title('Initial Conditions')
    
    subplot(3,2,4)
    latgrid = repmat(metrics.lath.data',[63 1]);
    plotdepth = default.depth;
    plotsalt= default.salt;
    eval(saltplotcmds);
    title('Default GOLD (40yr)')   
    
    subplot(3,2,6)
    plotdepth = griffies.depth;
    plotsalt = griffies.salt;
    eval(saltplotcmds);
    title('Initial Conditions')
    
    pause
%     saveas(gcf,...
%         sprintf('/ltraid3/ashao/GOLD/gold_aabw_diagnose/figs/esm2g_gold/%02d_%s.eps',lineno,wocenames{lineno}),'epsc');
    
end