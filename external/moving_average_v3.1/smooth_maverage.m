function [Y,Nsum,IND] = smooth_maverage(X,Fr,Fc,IND);
%SMOOTH_MAVERAGE   Smooths elements by moving average, ignoring NaN's.
%
%   Syntax:
%     [Y,Nsum,IND] = smooth_maverage(X,Fr,Fc,IND);
%
%   Input:
%     X   - Matrix with finite elements with/without Nan's.
%     Fr  - Window semi-length in the rows. A positive scalar (default 0).
%     Fc  - Window semi-length in the columns. A positive scalar (default
%           0). 
%     IND - Indicates de linear index of the elements to be smoothed.
%           By default it smooths the NaN's elements.  
%
%   Output:
%     Y    - X with the IND elements smoothed.
%     Nsum - Number of not NaN's elements that fixed on the moving window.
%            Provided to get a sum instead of a mean: Y(IND).*Nsum. Is a
%            vector of length equal as IND.
%     IND - Indicates de linear index of the elements that were smoothed.
%
%   Description: 
%      This program interpolates the elements defined by IND or the NaN's
%      ones, by averaging it along with the surrounding elements that fit
%      on the little matrix of size (2Fr+1)x(2Fc+1) centered on it and
%      ignoring NaN's. It smooths also in the edges. If Fc is 0 or empty,
%      the smoothing will be done columnwise; rowwise with Fr is 0 or
%      empty. If Fc is not specified and X is a row vector, it will be
%      smoothed by Fr.
%
%   Example:
%      x = round(rand(5)*10)
%      IND = [1 13 23 24]; 
%      x([1 13 17 23]) = NaN
%      smooth_maverage(x,1)
%      smooth_maverage(x,1,0,IND)
%      smooth_maverage(x,2,2,IND)
%      smooth_maverage(x,1,1)
%      smooth_maverage(x,[],1)
%
%
%   See also NANMEAN on the Statistical Toolbox and NANMOVING_AVERAGE,
%   NANMOVING_AVERAGE2 by Carlos Vargas.

% Copyright 2008  Carlos Vargas, nubeobscura@hotmail.com
%	$Revision: 1.0 $  $Date: 2008/03/04 11:00:00 $

%   Written by
%   M. in S. Carlos Adrián Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico,  march 2008
%
%   nubeobscura@hotmail.com
%
%   Download from:
%   http://www.mathworks.com/matlabcentral/fileexchange/loadAuthor.do?objec
%   tType=author&objectId=1093874

%% Errors checking
if ~nargin
 error('Interp2nanmovingaverage:Inputs','There are not inputs.')
end
if ndims(X)~=2
 error('Interp2nanmoving_average:Inputs','Entry must be a matrix.')
end
[M,N] = size(X);
if nargin<2 || isempty(Fr)
 Fr = 0;
end
if nargin<3 || isempty(Fc)
 Fc = 0;
 if M==1  % row vector?
  Fc = Fr;
  Fr = 0;
 end
end
if nargin<4 || isempty(IND)
 IND = 1:M*N; IND(~isnan(X(:))) = [];
 if isempty(IND)
  Y = X;
  return
 end
end

%% MAIN
Y = X;
Nind = length(IND);
ynans = repmat(NaN,2*Fr+1,2*Fc+1);
Ny = (2*Fr+1)*(2*Fc+1);
Nsum = repmat(NaN,Nind,1);
for k = 1:Nind
 [i,j] = ind2sub([M N],IND(k));
 rows = i-Fr:i+Fr; rowsi = ((rows>0)+(rows<M+1))>1;
 cols = j-Fc:j+Fc; colsi = ((cols>0)+(cols<N+1))>1;
 y = ynans;
 y(rowsi,colsi) = X(rows(rowsi),cols(colsi));
 nnan = ~isnan(y(:));
 Nsum(k)   = sum(  nnan); 
 Y(IND(k)) = sum(y(nnan))/Nsum(k);
end
% Resize in order to get the sum by Y(IND).*Nsum:
Nsum = reshape(Nsum,size(Y(IND)));

% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com