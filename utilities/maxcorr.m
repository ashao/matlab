function [corr lag] = maxcorr(A,B,varargin)


[C lags] = xcov(A,B,varargin{:});

[null maxidx] = max(abs(C));

corr = C(maxidx);
lag = lags(maxidx);