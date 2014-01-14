
trackpath = 'C:\Users\ashao\data\aviso\';
trackfiles = dir([trackpath '*.mat']);

figure(1)
clf
m_proj('Mercator','lon',[0 360],'lat',[-70 -30]);
hold on
for i=1:length(trackfiles)
    fprintf('%d...',i)
    load([trackpath trackfiles(i).name]);
    m_plot(track.lon,track.lat)
    
end
%%
m_grid;
m_coast('patch',[0 0 0]);
