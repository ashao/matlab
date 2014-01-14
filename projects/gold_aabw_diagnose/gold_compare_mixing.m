%% Set some user-modifiable variables

% Where to save the transects
outpath = '/ltraid3/ashao/uw-apl/figs/gold_aabw_diagnose/mask_srestore_under_ice_100yr/';
mkdir(outpath)
yidx = 79; % Year of run to output

% Paths to the model output (usually ocean_month.nc)
% Default gold run (gold_sis)
infiles.control = '/ltraid4/ashao/gold/default/default.newgold.00.40.nc';
% Modified run where all mixing parameters multiplied by 2
infiles.allconst = '/ltraid4/ashao/gold/MASK_SRESTORE_UNDER_ICE/ocean_month.nc';
% GOLD Initial Conditions
infiles.ic = '/ltraid4/ashao/save/GOLD/gold_aabw_diagnose/newgold/GOLD_IC.nc';

%% Extract the data

start4d = [39 0 0 0];
count4d = [1 inf inf inf];
temp.control = nc_varget(infiles.control,'temp',start4d,count4d);
salt.control = nc_varget(infiles.control,'salt',start4d,count4d);
h.control = nc_varget(infiles.control,'h',start4d,count4d);
depth.control = cumsum(h.control);
%%
start4d = [79 0 0 0];
temp.allconst= nc_varget(infiles.allconst,'temp',start4d,count4d);
salt.allconst = nc_varget(infiles.allconst,'salt',start4d,count4d);
h.allconst = nc_varget(infiles.allconst,'h',start4d,count4d);
depth.allconst = cumsum(h.allconst);
%%
temp.ic = nc_varget(infiles.ic,'Temp');
salt.ic = nc_varget(infiles.ic,'Salt');
depth.ic = cumsum( nc_varget(infiles.ic,'h'));
%% Set the plotting commands
load metrics
wocelabels = {'S1 68W','A23 30W','S2 0E','I6 30E','I8 270W','I9 245W', ...
    'S3 215W','P11 205W','P14 190W','P15 170W','P16 150W','P17 135W', ...
    'P18 105W','P19 90W'};
wocenames = {'S1','A23','S2','I6','I8','I9','S3','P11','P14','P15', ...
    'P16','P17','P18','P19'};
wocelons = [-68.5 -30.5 0.5 30.5 -270.5 -245.5 -215.5 -205.5 -190.5 ...
    -170.5 -150.5 -135.5 -105.5 -90.5];

templotcmds = ['plotdepth = squeeze(plotdepth(:,:,s3idx));' ...
    'plottemp= double(squeeze(plottemp(:,:,s3idx)));' ...    
    'contourf(latgrid,plotdepth,plottemp,-2:.4:10);' ...        
    'xlim([-75. -45]);ylim([0 2500]);' ...
    'set(gca,''ydir'',''reverse'');' ...
    'cax=colorbar;' ... 
    'ylabel(cax,''Temp'');' ...
    'caxis([-2 2]);' ...    
    'xlabel(wocelabels{lineno});'];

saltplotcmds = ['plotdepth = squeeze(plotdepth(:,:,s3idx));' ...
    'plotsalt= double(squeeze(plotsalt(:,:,s3idx)));' ...    
    'contourf(latgrid,plotdepth,plotsalt,33.5:0.05:35);' ...        
    'xlim([-75. -45]);ylim([0 2500]);' ...
    'set(gca,''ydir'',''reverse'');' ...
    'cax=colorbar;' ... 
    'ylabel(cax,''Salt'');' ...
    'caxis([34 35]);' ...    
    'xlabel(wocelabels{lineno});'];
%%

for lineno = 1:length(wocelons);
  figno = figure; % Comment out this if running Matlab without display
%       figno=figure('WindowStyle','Normal','Visible','Off') % Uncomment if running Matlab without display
    s3idx = find(abs(metrics.lonh.data- wocelons(lineno))<0.2)
    latgrid = repmat(metrics.lath.data',[63 1]);
    colormap(othercolor('BuDRd_12'))
    counter = 0;

    % Temperature
    subplot(3,2,1)
    plotdepth = depth.ic;
    plottemp = temp.ic;
    eval(templotcmds);
    xlabel('Initial Conditions')
    title(wocelabels{lineno})
    subplot(3,2,3)
    plotdepth = depth.control;
    plottemp = temp.control;
    eval(templotcmds);
        xlabel('Control (40yr)')
    
    subplot(3,2,5)
    plotdepth = depth.allconst;
    plottemp = temp.allconst;
    eval(templotcmds);
        xlabel('No restoring under ice (100yr)')
    
    % Salt
    subplot(3,2,2)
    plotdepth = depth.ic;    
    plotsalt = salt.ic;
    eval(saltplotcmds);
    xlabel('Initial Conditions')
    title(wocelabels{lineno})
    subplot(3,2,4)
    plotdepth = depth.control;
    plotsalt = salt.control;
    eval(saltplotcmds);
    xlabel('Control (20yr)')
    
    subplot(3,2,6)
    plotdepth = depth.allconst;
    plotsalt= salt.allconst;
    eval(saltplotcmds);
    xlabel('No restoring under ice (100yr)')        
    
    saveas(figno,[outpath sprintf('%02d_%s.eps',lineno,wocenames{lineno})],'epsc');
    
end
