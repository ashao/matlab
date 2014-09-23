inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/curvefit/processed2/';
trackpath = '/ltraid4/aviso/alongtrack/sla/vxxc_matlab/'

files = dir([inpath 't*.mat']);
nfiles = length(files);

for fileidx = 1:nfiles
    
    load([inpath files(fileidx).name]);
    load([trackpath files(fileidx).name]);
    ntrans = length(opt_track.skewness);
    
    for transidx = 1:ntrans
                
        opt_track.processed.lat(transidx) = opt_track.optpar{transidx}(2);
        opt_track.processed.lon(transidx) = ...
            interp1(track.lat,track.lon,opt_track.processed.lat(transidx));
        opt_track.processed.width(transidx) = ...
            sqrt( ...
            opt_track.optpar{transidx}(1)^2 + ... 
            opt_track.optpar{transidx}(3)^2 );
    end
    
    save([inpath files(fileidx).name],'opt_track')
    
end