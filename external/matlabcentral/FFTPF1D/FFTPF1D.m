% function [g] = FFTPF1D (X,binsize, f, P)
% Discrete Fourier Transform Low/High Pass Filter.
% This is a simply implement of such a filter for a given 1-D data.
% X: the array of you data, each data point is a bin of signal 
% binsize: the bin size in your data data
% f: the cutoff of wave length, 
% P: 1, or 0, true or false : true, low pass, which eliminating higher
% than f frequence signals, false, high pass, which eliminating lower
% than f frequence signals. 
% It could be easily modified into 2-D version, or translated into 
% R/S-language.
% Example:
% a=1:2:1000;
% b=sin(a) + sin(2.*a) + sin(0.1 .* a) + ...
% c=FFTPF1D(b, 2, 10, 1)
% Author: Zhihua Zhang. invokey@gmail.com
% Feel free to use and distribute this code for training or 
% academic propose.

function [g] = FFTPF1D(X, binsize, f, P)

    M=length(X);
    % get the fft position of cutoff frequence (reference
    % frequency)
    f = f / binsize;
    xidx=1:1:M;
    fftref = abs(fft(sin(xidx .* 2 .* pi ./f)));
    bounder = find( (max(fftref) - fftref) < (max(fftref) ./ 1000 ));
    % do the filter
    fftx = fft(X)
    if ( P )
        fftx(bounder(1):1:bounder(2))=0
    else
        fftx(1:1:bounder(1)) = 0;
        fftx(bounder(2):1:M) = 0;
    end

    g=real(ifft(fftx));    
end
