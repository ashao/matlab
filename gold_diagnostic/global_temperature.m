% Reduced Albedo Run
goldpath = '/ltraid3/ashao/gold_aabw_diagnose/albedo/50yr/';
goldfile.albedo.z = [goldpath 'ocean_month_z.nc'];
goldfile.albedo.isopycnal = [goldpath 'ocean_month.nc'];
goldfile.albedo.ice= [goldpath 'ice_month.nc'];
goldfile.albedo.static = [goldpath 'ocean_geometry.nc'];
goldfile.albedo.vertical = [goldpath 'Vertical_coordinate.nc'];
gold.albedo.time = nc_varget(goldfile.albedo.isopycnal,'time');
gold.albedo.layer = nc_varget(goldfile.albedo.vertical,'Layer');
gold.albedo.geolon = nc_varget(goldfile.albedo.static,'geolon');
gold.albedo.geolat = nc_varget(goldfile.albedo.static,'geolat');
gold.albedo.Ah = nc_varget(goldfile.albedo.static,'Ah');
gold.albedo.wet = logical(nc_varget(goldfile.albedo.static,'wet'));
dim.albedo.ntime = length(gold.albedo.time);
dim.albedo.nlayer =length(gold.albedo.layer);
[dim.albedo.nlon, dim.control.nlat]=size(gold.albedo.geolon);

% Control Run
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
dim.control.ntime = length(gold.control.time);
dim.control.nlayer =length(gold.control.layer);
[dim.control.nlon, dim.control.nlat]=size(gold.control.geolon);


%%
gold.albedo.avgtemp=zeros(dim.albedo.ntime);
gold.control.avgtemp=zeros(dim.albedo.ntime);

clear Ah3d
Ah3d(1,:,:)=gold.albedo.Ah;
Ah3d=repmat(Ah3d,[63 1 1]);

%%
for year = 1:dim.albedo.ntime
    
    fprintf('Year %d...',year);
    % Albedo
    fprintf('Albedo...')
    start4d = [year-1 0 0 0];
    count4d = [1 -1 -1 -1];
    temp = nc_varget(goldfile.albedo.isopycnal,'temp',start4d,count4d);
    depth = nc_varget(goldfile.albedo.isopycnal,'h',start4d,count4d);    
    wts = Ah3d.*depth;
    wts = wts/sum(makevec(wts));
    gold.albedo.avgtemp(year)=nansum(makevec(temp.*wts));
    fprintf('Done...');
    
    % Control        
    fprintf('Control...')
    start4d = [(year-1)*12 0 0 0];
    count4d = [12 -1 -1 -1];
    temp = squeeze(mean(...
        nc_varget(goldfile.control.isopycnal,'temp',start4d,count4d)));
    depth = squeeze(mean(...
        nc_varget(goldfile.control.isopycnal,'h',start4d,count4d)));
    wts = Ah3d.*depth;
    wts = wts/sum(makevec(wts));
    gold.control.avgtemp(year)=nansum(makevec(temp.*wts));
    fprintf('Done!\n')
    
end