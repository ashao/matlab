datapath = {'/ltraid4/aviso/alongtrack/tp_cf/'; ...
    '/ltraid4/aviso/alongtrack/j1_cf/' ; ...
    '/ltraid4/aviso/alongtrack/j2_cf/' };
outpath = '/ltraid4/aviso/alongtrack/tracks/';

for trackidx = 1:254 % Extract all 254 tracks
   
    track = aviso_track_timeseries(datapath,'vxxc',trackidx);
    outfile = [outpath sprintf('t%03d.mat',trackidx)];
    save(outfile,'track');
    
end

%%
