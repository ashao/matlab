
infiles.u10 = '/ltraid4/ecmwf/era-interim/4xdaily/u10.2002.2013.nc';
infiles.v10 = '/ltraid4/ecmwf/era-interim/4xdaily/v10.2002.2013.nc';

ecmwf.time = nc_varget(infiles.u10,'time');
ecmwf.lon = nc_varget(infiles.u10,'longitude');
ecmwf.lat = nc_varget(infiles.u10,'latitude');

%% Plot windspeed magnitude
ntime = length(ecmwf.time);
colormap(othercolor('Blues9'))

for t = 1:12:(ntime - 1)
    clf
    start3d = [t-1 0 0];
    count3d = [1 inf inf];
    ecmwf.u10 = nc_varget(infiles.u10,'u10',start3d,count3d);
    ecmwf.v10 = nc_varget(infiles.v10,'v10',start3d,count3d);
    fprintf('Done Loading %d\n',t)
    
%     subplot(2,2,4)
    m_proj('Stereographic','lon',240,'lat',77,'radius',20)
%     m_contourf(ecmwf.lon,ecmwf.lat,sqrt(ecmwf.u10.^2 + ecmwf.v10.^2),0:1:25,'LineColor','None');
    m_pcolor(ecmwf.lon,ecmwf.lat,sqrt(ecmwf.u10.^2 + ecmwf.v10.^2));
    shading flat
    m_grid('xtick',5,'ytick',5,'xaxislocation','top');
    m_coast('line','Color','Black');
    caxis([0 20])
    cax = colorbar;
    ylabel(cax,'Wind Speed (m/s)')
    time = double(ecmwf.time(t))/24 + datenum(1900,1,1);
    title(datestr(time,0))    
    drawnow; pause(0.1)
%     subplot(2,2,[1 2])
%     m_proj('Robinson','lon',[0 360],'lat',[-80 80])
%     m_contourf(ecmwf.lon,ecmwf.lat,sqrt(ecmwf.u10.^2 + ecmwf.v10.^2),0:1:25,'LineColor','None');
%     %     shading flat
%     m_grid;
%     m_coast('line','Color','Black');
%     caxis([0 20])
%     %     cax = colorbar;
%     %     ylabel(cax,'Wind Speed (m/s)')
%     time = double(ecmwf.time(t))/24 + datenum(1900,1,1);
%     title(datestr(time,0))
%     
%     subplot(2,2,3)
%     m_proj('Stereographic','lon',0,'lat',-90,'radius',40)
%     m_contourf(ecmwf.lon,ecmwf.lat,sqrt(ecmwf.u10.^2 + ecmwf.v10.^2),0:1:25,'LineColor','None');
%     %     shading flat
%     m_grid('xtick',5,'ytick',5);
%     ;
%     m_coast('line','Color','Black');
%     caxis([0 20])
%     %     cax = colorbar;
%     %     ylabel(cax,'Wind Speed (m/s)')
%     time = double(ecmwf.time(t))/24 + datenum(1900,1,1);

%     pause
end
