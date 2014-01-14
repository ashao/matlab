function [ ] = run_lags( inpath , outpath, lags )
matlabpool(4)	
	nlags=length(lags);
	for ii=1:nlags
		fullinpath=[inpath filesep sprintf('lag_%02d',lags(ii)) filesep ];
		fulloutpath=[outpath filesep sprintf('lag_%02d',lags(ii)) filesep ];
		nfiles=length(dir([fullinpath filesep '*.mat']));
		nfiles
        parfor jj=1:nfiles
            optcruncher_2011_diff(fullinpath, fulloutpath,jj,jj);
        end
    end
    matlabpool close
end
