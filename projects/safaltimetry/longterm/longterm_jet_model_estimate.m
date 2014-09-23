function longterm_jet_model_estimate()
inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/range/toprocess/';
outpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/curvefit/processed/';
if ~exist(outpath)
	mkdir(outpath)
end
files = dir([inpath 't*.mat']);
nfiles = length(files);
options = optimset('Display','iter','TolFun',1e-12);

for i=1:nfiles
   
    load([inpath files(i).name]);
    ntrans = length(opt_track.skewness);
    opt_track.tracknum = i;
    for tidx = 1:ntrans
        datay = opt_track.lat{tidx};
        dataskew = opt_track.skewness{tidx};
        meansla = nanmean(opt_track.mean{tidx});
        SNR = meansla./opt_track.std{tidx};
%         x0 = [0.5 mean(datay) 0.5];
%         LB = [0.01 x0(2)-5 0.01];
%         UB = [5 x0(2)+5 5];
%         [opt_track.optpar{tidx} opt_track.fval{tidx} opt_track.simskew{tidx}] = ...
%             jet_model_estimate(datay',dataskew,SNR,LB,UB,x0);
%         opt_track.R2{tidx} = 1-opt_track.fval{tidx}./nanvar(dataskew)/length(dataskew);

% Do this with the curvefit to the Gaussian derivative instead
         x0 = [mean(datay) 2 10];
         LB = [x0(1)-5 0.01  0.01];
         UB = [x0(1)+5 20 20];
        if length(datay)>3

		try
      		  [opt_track.optpar{tidx} opt_track.fval{tidx}] = lsqcurvefit(@gaussian_derivative,x0,datay,dataskew',LB,UB,options);
        


      			  fprintf('Track %d Transition: %d Parameters: %e %e %e\n', ...
           		 opt_track.tracknum,tidx,opt_track.optpar{tidx})        
		catch err
			fprintf('Error in Track %d:',i)
		end
	end
%         clf; hold on;
%          plot(datay,dataskew)
%          plot(datay,gaussian_derivative(opt_track.optpar{tidx},datay),'-k')
%          pause
%          plot(datay,dataskew,'k')
%          drawnow
        
    end
    save([outpath files(i).name],'opt_track')
    
    
end
