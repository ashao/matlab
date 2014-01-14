function [ ] = postprocess_all_lags(inpath)

    lags=-60:5:60;
    for ii=1:length(lags)
       
        fullpath=[inpath sprintf('lag_%02d',lags(ii)) filesep];
        disp(sprintf('Lag %d',lags(ii)))
        postprocess_diff(fullpath);
        
    end

end