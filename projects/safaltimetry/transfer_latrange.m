old_track_path = 'C:\Users\ashao\Data\altimetry\longterm\';
new_track_path = 'C:\Users\ashao\Data\vxxc_matlab\';

oldfiles = dir([old_track_path 't*.mat']);
newfiles = dir([new_track_path 't*.mat']);

ntracks = length(oldfiles);
%%
for tidx = 1:ntracks
    
    fprintf('Track %03d\n',tidx)
    old_track = load([old_track_path oldfiles(tidx).name]);
    load([new_track_path newfiles(tidx).name])
    
    track.latrange = old_track.track.latrange;    
    track.number = tidx;    
    save([new_track_path newfiles(tidx).name],'track')    
    
end