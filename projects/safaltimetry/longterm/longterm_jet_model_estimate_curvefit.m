% function longterm_jet_model_estimate_curvefit()
inpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/toprocess/';
outpath = '/ltraid3/ashao/uw-apl/projects/saf_altimetry/longterm/curvefit/processed2/';
if ~exist(outpath)
    mkdir(outpath)
end
files = dir([inpath 't*.mat']);
nfiles = length(files);
options = optimset('Display','Off','TolFun',1e-12);

for i=1:nfiles
    
    load([inpath files(i).name]);
    ntrans = length(opt_track.skewness);
    opt_track.tracknum = i;
    for tidx = 1:ntrans
        datay = opt_track.lat{tidx}';
        dataskew = opt_track.skewness{tidx};
        notnan = ~isnan(dataskew);
        datay = datay(notnan);
        dataskew = dataskew(notnan);
        meansla = nanmean(opt_track.mean{tidx});
%         SNR = abs(meansla./opt_track.std{tidx});
    SNR = 0.7;
        %         x0 = [0.5 mean(datay) 0.5];
        %         LB = [0.01 x0(2)-5 0.01];
        %         UB = [5 x0(2)+5 5];
        %         [opt_track.optpar{tidx} opt_track.fval{tidx} opt_track.simskew{tidx}] = ...
        %             jet_model_estimate(datay',dataskew,SNR,LB,UB,x0);
        %         opt_track.R2{tidx} = 1-opt_track.fval{tidx}./nanvar(dataskew)/length(dataskew);
        
        % Do this with the curvefit to the Gaussian derivative instead
        
        if length(datay)>3            
%             try
                x0 = [1 mean(datay) 1 0.7];
                LB = [0.01 x0(2)-10 0.01 0.01];
                UB = [10 x0(2)+10 10 1];
                [opt_track.optpar{tidx} opt_track.R2{tidx}] = ...
                    jet_model_estimate_curvefit(datay,dataskew,LB,UB,x0);
                fprintf('Track %d Transition: %d SNR: %f Parameters: %e %e %e %f R^2: %f\n', ...
                    opt_track.tracknum,tidx,SNR,opt_track.optpar{tidx},opt_track.R2{tidx})
%                 ylim([-1 1])
%                title(sprintf('R^2 = %f',opt_track.R2{tidx})); axis equal;
%                drawnow;
%                pause(1)
%                 pause(1)
%                 drawnow
%                 pause(1)
%                 pause
%             catch
%                 fprintf('Error in Track %d:',i)
%             end
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
