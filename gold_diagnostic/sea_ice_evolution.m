goldpathsurf.cfc11.soldt = '/ltraid3/ashao/gold_aabw_diagnose/albedo/50yr/';

goldfile.z = [goldpath 'ocean_month_z.nc'];
goldfile.isopycnal = [goldpath 'ocean_month.nc'];
goldfile.ice= [goldpath 'ice_month.nc'];
goldfile.static = [goldpath 'ocean_geometry.nc'];
goldfile.vertical = [goldpath 'Vertical_coordinate.nc'];

gold.time = nc_varget(goldfile.ice,'time');
gold.layer = nc_varget(goldfile.vertical,'Layer');
gold.geolon = nc_varget(goldfile.static,'geolon');
gold.geolat = nc_varget(goldfile.static,'geolat');
gold.wet = logical(nc_varget(goldfile.static,'wet'));
%%
midx = 7;
load metrics
himpath = '/ltraid1/ashao/HIM/hyak_store/COMBINE/month/';
him.alb = squeeze(mean(nc_varget([himpath 'ice_month.nc'],...
    'ALB',[0 0 0 ],[120 -1 -1])));
% him.alb = squeeze(max(him.cn));
him.mincn = squeeze(min(him.cn));
him.meancn = squeeze(mean(him.cn));


%%
dim.ntime = length(gold.time);
dim.nlayer =length(gold.layer);
[dim.nlon, dim.nlat]=size(gold.geolon);

%% 
m_proj('Orthographic','lat',-90,'radius',50);
colormap(othercolor('BuDRd_12'));
for t=1:(dim.ntime)
    clf
    hold on;
    gold.icecn = squeeze((nc_varget(goldfile.ice,'ALB',...
        [t-1 0 0],[1 -1 -1]))); 
    m_pcolor(gold.geolon,gold.geolat,gold.icecn-him.alb);
%     m_contour(gold.geolon,gold.geolat,him.meancn,[0.5 0.5],'linecolor','k','LineWidth',2)
    
    shading flat;
    title(sprintf('Model Year %d',round(gold.time(t)/365)));
    m_grid;
    m_coast('patch',[0.5 0.5 0.5]);    
    caxis([-0.25 0.25])
    colorbar
    drawnow;
    
end