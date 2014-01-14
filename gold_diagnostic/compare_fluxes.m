goldpath = '/ltraid2/darr/GOLDruns_hyak/700yr-nz63/';
goldpath = '/ltraid3/ashao/gold_aabw_diagnose/albedo/50yr/';
goldfile.z = [goldpath 'ocean_month_z.nc'];
goldfile.isopycnal = [goldpath 'ocean_month.nc'];
goldfile.ice= [goldpath 'ice_month.nc'];
goldfile.static = [goldpath 'ocean_geometry.nc'];
goldfile.vertical = [goldpath 'Vertical_coordinate.nc'];

gold.time = nc_varget(goldfile.ice,'time');
% gold.layer = nc_varget(goldfile.vertical,'Layer');
gold.geolon = nc_varget(goldfile.static,'geolon');
gold.geolat = nc_varget(goldfile.static,'geolat');
gold.wet = logical(nc_varget(goldfile.static,'wet'));

fields = {'sw','lw','sh','lh'};
ncfields = upper(fields);

%%
load metrics
himpath = '/ltraid1/ashao/HIM/hyak_store/COMBINE/month/';
him.sw = squeeze(mean(nc_varget([himpath 'ice_month.nc'],...
    'SW',[0 0 0 ],[-1 -1 -1 ])));
him.lw = squeeze(mean(nc_varget([himpath 'ice_month.nc'],...
    'LW',[0 0 0 ],[-1 -1 -1 ])));
him.sh = squeeze(mean(nc_varget([himpath 'ice_month.nc'],...
    'SH',[0 0 0 ],[-1 -1 -1 ])));
him.lh = squeeze(mean(nc_varget([himpath 'ice_month.nc'],...
    'LH',[0 0 0 ],[-1 -1 -1 ])));
him.qnet = him.sw+him.lw+him.sh+him.lh;

%%
dim.ntime = length(gold.time);
dim.nlayer =length(gold.layer);
[dim.nlon, dim.nlat]=size(gold.geolon);
dim.nfields = length(fields);
%%
plotcmds = [ ['shading flat;'] ...
    ['title(sprintf(''Model Year %d'',round(gold.time(yidx)/365)));']...
    ['m_grid; m_coast(''patch'',[0.5 0.5 0.5]);']...
    ['cax=colorbar(''Location'',''SouthOutside'');']];
%%
m_proj('Orthographic','lat',-90,'lon',0,'radius',50);
colormap(othercolor('BuDRd_12'));
% for t=1:(dim.ntime/12)
for yidx=1:dim.ntime
%     yidx = (t-1)*12;    
    gold.sw = nc_varget(goldfile.ice,...
        'SW',[yidx 0 0 ],[1 -1 -1 ]);
    gold.lw = nc_varget(goldfile.ice,...
        'LW',[yidx 0 0 ],[1 -1 -1 ]);
    gold.sh = nc_varget(goldfile.ice,...
        'SH',[yidx 0 0 ],[1 -1 -1 ]);
    gold.lh = nc_varget(goldfile.ice,...
        'LH',[yidx 0 0 ],[1 -1 -1 ]);
    gold.qnet = gold.sw+gold.lw+gold.sh+gold.lh;
%     for i = 1:dim.nfields
%         subplot(2,2,i)
        m_pcolor(gold.geolon,gold.geolat,gold.qnet-him.qnet);
        eval(plotcmds)
        caxis([-100 100])
        xlabel(cax,'Q_{net}')
        drawnow;
%     end
    
   
end
