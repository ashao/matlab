function [ ] = timevary_jetmodel_estimate_even
inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/timevary/annual/toprocess/';
outpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/timevary/annual/processed/curvefit/';
if ~exist(outpath)
    mkdir(outpath)
end
files = dir([inpath 't*.mat']);
nfiles = length(files)

for fileidx=2:2:nfiles
    load([inpath files(fileidx).name]);
    nwindows = length(opt_track);
    for winidx = 1:nwindows
        
        ntrans = length(opt_track(winidx).skewness);
        for tidx = 1:ntrans
            
            try
                datay = opt_track(winidx).lat{tidx};
                dataskew = opt_track(winidx).skewness{tidx};
                meansla = nanmean(opt_track(winidx).mean{tidx});
                SNR = meansla./opt_track(winidx).std{tidx}*1.5;
                x0 = [0.5 mean(datay) 0.5 0.5];
                LB = [0.01 x0(2)-5 0.01 0.01];
                UB = [5 x0(2)+5 5 1];
                
                [opt_track(winidx).optpar{tidx} opt_track(winidx).R2{tidx}] = ...
                    jet_model_estimate_curvefit(datay',dataskew,LB,UB,x0);
                %             opt_track(winidx).R2{tidx} = 1-opt_track(winidx).fval{tidx}./nanvar(dataskew)/length(dataskew);
                fprintf('Track %d Window: %d Transition: %d R2: %f Parameters: %e %e %e %e\n', ...
                    opt_track(winidx).tracknum,winidx,tidx,opt_track(winidx).R2{tidx},opt_track(winidx).optpar{tidx})
            end
        end
    end
    
    save([outpath files(fileidx).name],'opt_track')
    
end
