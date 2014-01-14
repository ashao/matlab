infiles.full = '/ltraid3/ashao/uw-apl/projects/stanley/saturation.coldstart.fullstanley.y90.100.nc';
infiles.fge = '/ltraid3/ashao/uw-apl/projects/stanley/saturation.coldstart.fge_only.y90.100.nc';
infiles.fge_fc = '/ltraid3/ashao/uw-apl/projects/stanley/saturation.coldstart.fge_fc.y90.y100.nc';
infiles.aux = '/ltraid3/ashao/uw-apl/data/offtrac/input/normalyear/auxiliary_input.clim.nc';

%%
[ar n2 o2] = ar_n2_o2_properties;

start4d = [0 0 0 0];
count4d = [12 1 -1 -1];

output.full.ar = nc_varget(infiles.full,'mn_ar',start4d,count4d);
output.full.n2 = nc_varget(infiles.full,'mn_n2',start4d,count4d);
output.full.o2 = nc_varget(infiles.full,'mn_o2',start4d,count4d);

output.fge.ar = nc_varget(infiles.fge,'mn_ar',start4d,count4d);
output.fge.n2 = nc_varget(infiles.fge,'mn_n2',start4d,count4d);
output.fge.o2 = nc_varget(infiles.fge,'mn_o2',start4d,count4d);


output.fge_fc.ar = nc_varget(infiles.fge_fc,'mn_ar',start4d,count4d);
output.fge_fc.n2 = nc_varget(infiles.fge_fc,'mn_n2',start4d,count4d);
output.fge_fc.o2 = nc_varget(infiles.fge_fc,'mn_o2',start4d,count4d);

forcing.temp = nc_varget(infiles.aux,'temp',start4d,count4d);
forcing.salt = nc_varget(infiles.aux,'salt',start4d,count4d);
% forcing.temp = forcing.temp([2:12 1],:,:);
% forcing.salt = forcing.salt([2:12 1],:,:);
%%
satconc.ar = ar.atmconc*trac_calcsol(forcing.temp+273.15,forcing.salt,ar.F_sol_coeffs.vol);
satconc.n2 = n2.atmconc*trac_calcsol(forcing.temp+273.15,forcing.salt,n2.F_sol_coeffs.vol);
satconc.o2 = o2.atmconc*trac_calcsol(forcing.temp+273.15,forcing.salt,o2.F_sol_coeffs.vol);

surfsat.full.ar = output.full.ar./satconc.ar;
surfsat.full.n2 = output.full.n2./satconc.n2;
surfsat.full.o2 = output.full.o2./satconc.o2;

surfsat.fge.ar = output.fge.ar./satconc.ar;
surfsat.fge.n2 = output.fge.n2./satconc.n2;
surfsat.fge.o2 = output.fge.o2./satconc.o2;

surfsat.fge_fc.ar = output.fge_fc.ar./satconc.ar;
surfsat.fge_fc.n2 = output.fge_fc.n2./satconc.n2;
surfsat.fge_fc.o2 = output.fge_fc.o2./satconc.o2;
%% Surface maps of saturation
load metrics
figure(1)
m_proj('Mercator','lon',[-240 -120],'lat',[20 65]);
colormap(othercolor('BuDRd_12'))

    subplot(3,2,1);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.full.ar(2,:,:),0.90:0.01:1.10);   
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Argon Saturation February')
    colorbar;
    
    subplot(3,2,2);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.full.ar(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Argon Saturation June')
   colorbar; 
    
   
    subplot(3,2,3);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge_fc.ar(2,:,:),0.90:0.01:1.10);   
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge Fc')
    colorbar;
    
    subplot(3,2,4);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge_fc.ar(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge Fc')
   colorbar; 
   
    subplot(3,2,5);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge.ar(2,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge')
colorbar;
    subplot(3,2,6);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge.ar(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge')
    colorbar;
    
    saveas(gcf,'figs/argon_surfsat.eps','epsc')
    %% Surface maps of saturation
load metrics
figure(1)
m_proj('Mercator','lon',[-240 -120],'lat',[20 65]);
colormap(othercolor('BuDRd_12'))

    subplot(3,2,1);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.full.n2(2,:,:),0.90:0.01:1.10);   
    clabel(cs,v);
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Nitrogen Saturation February')
    colorbar;
    
    subplot(3,2,2);    
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.full.n2(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v);
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Nitrogen Saturation June')
   colorbar; 
    subplot(3,2,3);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge_fc.n2(2,:,:),0.90:0.01:1.10);   
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge Fc')
    colorbar;
    
    subplot(3,2,4);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge_fc.n2(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge Fc')
   colorbar; 
    subplot(3,2,5);    
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge.n2(2,:,:),0.90:0.01:1.10);    
    clabel(cs,v);
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge')
colorbar;
    subplot(3,2,6);    
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge.n2(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v);
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge')
    colorbar;
    saveas(gcf,'figs/n2_surfsat.eps','epsc')
        %% Surface maps of saturation
load metrics
figure(1)
m_proj('Mercator','lon',[-240 -120],'lat',[20 65]);
colormap(othercolor('BuDRd_12'))

    subplot(3,2,1);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.full.o2(2,:,:),0.90:0.01:1.10);   
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('OxygenSaturation February')
    colorbar;
    
    subplot(3,2,2);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.full.o2(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Oxygen Saturation June')
   colorbar; 
    
       subplot(3,2,3);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge_fc.ar(2,:,:),0.90:0.01:1.10);   
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge Fc')
    colorbar;
    
    subplot(3,2,4);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge_fc.ar(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge Fc')
   colorbar; 
   
    subplot(3,2,5);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge.o2(2,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge')
colorbar;
    subplot(3,2,6);
    [cs v]=m_contourf(metrics.geolon.data,metrics.geolat.data,surfsat.fge.o2(6,:,:),0.90:0.01:1.10);    
    clabel(cs,v)
    caxis([0.90 1.10])
    m_coast('patch',[0 0 0]);
    m_grid;
    title('Fge')
    colorbar;
    saveas(gcf,'figs/o2_surfsat.eps','epsc')
    