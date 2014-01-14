m_proj('Robinson','lon',[0 360],'lat',[-80 80]);
pathsla = '/ltraid2/ashao/uw-apl/data/aviso/msla/';
files = dir([pathsla '*.nc']);
nfiles = length(files);
ndelim = 8;
outdir = [pathsla filesep 'image' filesep];
mkdir(outdir)
%%

matlabpool(6)
%%
for ifile = 1:nfiles
    
    figure('visible','off')
    filename = files(ifile).name;
    tempstr = textscan(filename,'%s','Delimiter','_');
    timestr = char(tempstr{1}(8));
    date = datenum(str2num(timestr(1:4)),str2num(timestr(5:6)),str2num(timestr(7:8)));
    lon = nc_varget([pathsla filename],'NbLongitudes');
    lat = nc_varget([pathsla filename],'NbLatitudes');
    sla = nc_varget([pathsla filename],'Grid_0001');
    m_pcolor(lon,lat,sla');
    shading flat;
    caxis([-25 25])
    colormap(othercolor('BuDRd_12'));
    cax=colorbar('SouthOutside');
    xlabel(cax,'SLA (cm)')
    m_coast('patch',[0.5 0.5 0.5]);
    m_grid('box','fancy')
    title(datestr(date))
    exportfig(gcf,[outdir datestr(date,29) '.png'],'Format','png','color','cmyk','resolution',150)
    close(gcf)
    disp(datestr(date))
end