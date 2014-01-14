function [ clim ] = him_make_clim( himfiles, outpath )
% Makes climatology from HIM output by averaging over monthly values
% himfiles:
%   Structure containing the following fields which are paths to files
%   ocmon_snap: Snapshot of data at beginning of month
%   ocean_month: Monthly averages
%   ice_month:  Ice model output
%   flux_month: Atmospheric fluxing...requires interpolation onto HIM grid

if ~exist(outpath,'dir')
    mkdir(outpath)
end

debug=false;

%% Define dimensions of HIM run
nmon=12;
if debug
    nmon=1;
end
nlay=49;
nlat=210;
nlon=360;
sepfiles=true;

%% Initialize arrays
fill4d=zeros(nmon,nlay,nlat,nlon);
fill3d=zeros(nmon,nlat,nlon);

fice_clim=fill3d;
atmp_clim=fill3d;

temp_clim=fill4d;
salt_clim=fill4d;
h_clim=fill4d;
uh_clim=fill4d;
vh_clim=fill4d;
wd_clim=zeros(nmon,nlay+1,nlat,nlon);

count_4d=-[ 1 1 1 1 ];
count_3d=-[ 1 1 1 ];

stride_4d=[12 1 1 1];
stride_3d=[12 1 1];

%% Extract data and average over every month
for mon=1:nmon
    
    start_4d=[ mon-1 0 0 0 ];
    start_3d=[ mon-1 0 0 ];
    
    
    disp(sprintf('Month %d/%d',mon,12))
    % Ocean Month variables
    temp_mon=nc_varget(himfiles.ocean_month,'temp',start_4d,count_4d,...
        stride_4d);
    salt_mon=nc_varget(himfiles.ocean_month,'salt',start_4d,count_4d,...
        stride_4d);
    uh_mon=nc_varget(himfiles.ocean_month,'uh',start_4d,count_4d,...
        stride_4d);
    vh_mon=nc_varget(himfiles.ocean_month,'vh',start_4d,count_4d,...
        stride_4d);
    wd_mon=nc_varget(himfiles.ocean_month,'wd',start_4d,count_4d,...
        stride_4d);
    h_mon=nc_varget(himfiles.ocmon_snap,'h',start_4d,count_4d,stride_4d);
            
    temp_clim(mon,:,:,:)=squeeze(mean(temp_mon,1));
    salt_clim(mon,:,:,:)=squeeze(mean(salt_mon,1));
    uh_clim(mon,:,:,:)=squeeze(mean(uh_mon,1));
    vh_clim(mon,:,:,:)=squeeze(mean(vh_mon,1));
    wd_clim(mon,:,:,:)=squeeze(mean(wd_mon,1));
    h_clim(mon,:,:,:)=squeeze(mean(h_mon,1));
end

%% Extract dimensions (time,depth,lat,lon) for use in creating output files
time_snap=nc_varget(himfiles.ocmon_snap,'Time',0,12);
time_mean=nc_varget(himfiles.ocean_month,'Time',0,12);
lath=nc_varget(himfiles.ocean_month,'xh');
lonh=nc_varget(himfiles.ocean_month,'yh');
zl=nc_varget(himfiles.ocean_month,'zl');
zi=nc_varget(himfiles.ocean_month,'zi');

geninfo.Datatype='single';
geninfo.Dimension={ 'Time', 'zl', 'lonh', 'lath' };

timeinfo.Name='Time';
timeinfo.Datatype='single';
timeinfo.Dimension={ 'Time' };

latinfo.Name='lath';
latinfo.Datatype='single';
latinfo.Dimension={ 'lath' };

loninfo.Name='lonh';
loninfo.Datatype='single';
loninfo.Dimension={ 'lonh' };

zlinfo.Name='zl';
zlinfo.Datatype='single';
zlinfo.Dimension={ 'zl' };

ziinfo.Name='zi';
ziinfo.Datatype='single';
ziinfo.Dimension={ 'zi' };

if sepfiles==true
    
    wdoutpath=[outpath filesep 'WD-clim.nc'];
    tsoutpath=[outpath filesep 'ts.nc'];
    houtpath=[outpath filesep 'H-clim.nc'];
    uhoutpath=[outpath filesep 'UH-clim.nc'];
    vhoutpath=[outpath filesep 'VH-clim.nc'];
    
    %% Temperature/Salinity Climatology
    
    make_new_netcdf(tsoutpath)
    
    tempinfo=geninfo;
    tempinfo.Name='TEMPCLIM';
    saltinfo=geninfo;
    saltinfo.Name='SALTCLIM';
    
    nc_addvar(tsoutpath,tempinfo);
    nc_addvar(tsoutpath,saltinfo);
    
    nc_varput(tsoutpath,'TEMPCLIM',single(temp_clim))
    nc_varput(tsoutpath,'SALTCLIM',single(salt_clim))
    
    %% Snapshot of Thickness Fields Climatology
    
    make_new_netcdf(houtpath);
    hinfo=geninfo;
    hinfo.Name='HCLIM';
    
    nc_addvar(houtpath,hinfo);
    nc_varput(houtpath,'HCLIM',single(h_clim))
    
    %% Mass Transport Climatology
    uhinfo=geninfo;
    uhinfo.Name='UHCLIM';
    
    
    make_new_netcdf(uhoutpath);
    nc_addvar(uhoutpath,uhinfo);
    nc_varput(uhoutpath,'UHCLIM',single(uh_clim));
    
    vhinfo=geninfo;
    vhinfo.Name='VHCLIM';
    
    make_new_netcdf(vhoutpath);
    nc_addvar(vhoutpath,vhinfo);
    nc_varput(vhoutpath,'VHCLIM',single(vh_clim));
    
    wdinfo=geninfo;
    wdinfo.Name='WDCLIM';
    
    make_new_netcdf(wdoutpath);
    wdinfo.Dimension{2}='zi';
    nc_addvar(wdoutpath,wdinfo);
    nc_varput(wdoutpath,'WDCLIM',single(wd_clim));

end

clim.H=h_clim;
clim.T=temp_clim;
clim.S=salt_clim;
clim.UH=uh_clim;
clim.VH=vh_clim;
clim.WD=wd_clim;


    function [ ] = make_new_netcdf(filepath,yeslay)
        
        if nargin<2
            yeslay=true;
        end
        
        nc_create_empty( filepath, 'clobber');
        nc_adddim(filepath,'Time',length(time_mean));
        if strcmp(wdoutpath,filepath)
            nc_adddim(filepath,'zi',length(zi));
            nc_addvar(filepath,ziinfo);
            nc_varput(filepath,'zi',zi);
        elseif yeslay
            nc_adddim(filepath,'zl',length(zl));
            nc_addvar(filepath,zlinfo);
            nc_varput(filepath,'zl',zl);
        end
        nc_adddim(filepath,'lath',length(lath));
        nc_adddim(filepath,'lonh',length(lonh));
        
        nc_addvar(filepath,timeinfo);
        
        nc_addvar(filepath,latinfo);
        nc_addvar(filepath,loninfo);
        
        nc_varput(filepath,'Time',time_mean);
        nc_varput(filepath,'lath',lath);
        nc_varput(filepath,'lonh',lonh);
        
    end
end