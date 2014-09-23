inpath = '/ltraid3/ashao/uw-apl/data/safaltimetry/vxxc_matlab/';
load('/ltraid3/ashao/uw-apl/matlab/projects/safaltimetry/sokolov/sokolov.mat');
files = dir([inpath 't*.mat']);
nfiles = length(files)

pfm = sokolov.pfm';
safm = sokolov.safm';

for t = 1:nfiles
    
    load([inpath files(t).name])
    clf;
    worldmap([-90 -30],[0 360])
    
    % Transform into map projection to account for 0-360 loop
    [pfmx pfmy] = mfwdtran(sokolov.pfm(:,2),sokolov.pfm(:,1));
    [safmx safmy] = mfwdtran(sokolov.safm(:,2),sokolov.safm(:,1));        
    [trackx tracky] = mfwdtran(track.lat,track.lon);
    
    % Intersection with Polar Front
    p=InterX([pfmx' ; pfmy'],[trackx' ; tracky']);
    [intersectlat intersectlon] = minvtran(p(1,:),p(2,:));
    [southbound sidx] = min(intersectlat);
    southlon = intersectlon(sidx);
    
    % Intersection with SAF
    p=InterX([safmx' ; safmy'],[trackx' ; tracky']);
    [intersectlat intersectlon] = minvtran(p(1,:),p(2,:));
    [northbound nidx] = max(intersectlat);
    northlon = intersectlon(nidx);
    
    % Verify before saving
    plotm(track.lat,track.lon)
    plotm(sokolov.pfm(:,2),sokolov.pfm(:,1))
    plotm(sokolov.safm(:,2),sokolov.safm(:,1))
    scatterm(southbound,southlon);
    scatterm(northbound,northlon);
    track.latrange(1) = southbound;
    track.latrange(2) = northbound;      
    
    drawnow
    save([inpath files(t).name],'track')
    
    
    
    
end