inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/timevary/toprocess/';
outpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/timevary/processed/curvefit/';
if ~exist(outpath)
    mkdir(outpath)
end
files = dir([inpath 't*.mat']);
nfiles = length(files)

for i=1:nfiles
    
    load([inpath files(i).name]);
    nwindows = length(opt_track);
    for winidx = 1:nwindows
        
        ntrans = length(opt_track(winidx).skewness);
        for tidx = 1:ntrans
            datay = opt_track(winidx).lat{tidx};
            dataskew = opt_track(winidx).skewness{tidx};
            SNR = opt_track(winidx).std{tidx};
            x0 = [0.5 mean(datay) 0.5];
            LB = [0.01 x0(2)-10 0.01];
            UB = [10 x0(2)+10 10];
            [opt_track(winidx).optpar{tidx} opt_track(winidx).R2{tidx}] = ...
                jet_model_estimate(datay',dataskew,SNR,LB,UB,x0);
            fprintf('Track %d Window: %d Transition: %d Parameters: %e %e %e\n', ...
                opt_track(winidx).tracknum,winidx,tidx,opt_track(winidx).optpar{i})
            
        end
    end
    
    save([outpath files(i).name],'opt_track')
    
end
