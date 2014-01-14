function [gh11,nh11]=variof1(x1,icode);
%
% function [gh11,nh11]=variof1(x1,icode);
%
% function to compute variograms or covariograms, in 1D or 2D
% the data are on a (possibly incomplete) regular grid.
% missing values are entered as NaN into the regular grid matrix
% the program computes variograms in the frequency domain by
% using 2D-FFT.
%
% input: x1:    data matrix. Missing values are indicated by NaN
%        icode: a code to indicate which function to compute
%               =1 : variogram
%               =2 : covariogram
%
% output: gh11: variogram or covariogram depending on icode.
%         nh11: number of pairs available
%
% this program uses the functions FFT2, IFFT2, FFT2SHIFT and CONJ which are
% standard MATLAB functions.
% author: D. Marcotte, dmarcotte@mail.polymtl.ca

[n,p]=size(x1);				% dimensions of data matrix
nrows=2*n-1;
ncols=2*p-1;

% find the closest multiple of 8 to obtain a good compromise between
% speed (a power of 2) and memory required

nr2=ceil(nrows/8)*8;
nc2=ceil(ncols/8)*8;

% form an indicator  matrix: 1's for all data values
%                            0's for missing values
% in data matrix, replace missing values by 0;

x1id=~isnan(x1);			% 1 for a data value; 0 for missing
x1(~x1id)=zeros(sum(sum(~x1id)),1);	% missing replaced by 0

fx1=fft2(x1,nr2,nc2);			% fourier transform of x1

if icode==1
  fx1_x1=fft2(x1.*x1,nr2,nc2);		% fourier transform of x1*x1
end
clear x1;

fx1id=fft2(x1id,nr2,nc2);		% fourier transform of the indicator matrix
clear x1id

% compute number of pairs at all lags

nh11=round(real(ifft2(conj(fx1id).*fx1id)));

% compute the different structural functions according to icode

if icode==1;				% variogram is computed
  gh11=real(ifft2(conj(fx1id).*fx1_x1+conj(fx1_x1).*fx1id-2*conj(fx1).*fx1));
  gh11=gh11./max(nh11,1)/2;

else					% covariogram is computed

  m1=real(ifft2(conj(fx1).*fx1id))./max(nh11,1);	% compute tail mean
  m2=real(ifft2(conj(fx1id).*fx1))./max(nh11,1);	% compute head mean
  clear fx1id
  gh11=real(ifft2(conj(fx1).*fx1));
  gh11=gh11./max(nh11,1)-m1.*m2;
end

clear fx1 fx1id fx1_fx1

% reduce matrix to required size and shift so that the 0 lag appears at the center of each matrix

nh11=[nh11(1:n,1:p) nh11(1:n,nc2-p+2:nc2);nh11(nr2-n+2:nr2,1:p) nh11(nr2-n+2:nr2,nc2-p+2:nc2)];
gh11=[gh11(1:n,1:p) gh11(1:n,nc2-p+2:nc2);gh11(nr2-n+2:nr2,1:p) gh11(nr2-n+2:nr2,nc2-p+2:nc2)];

gh11=fftshift(gh11);
nh11=fftshift(nh11);
