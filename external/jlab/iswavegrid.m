function[bool,dk]=iswavegrid(k)
%ISWAVEGRID  Tests whether a matrix is in WAVEGRID format
%
%   ISWAVEGRID(K) returns true if K is a matrix in the WAVEGRID
%   format, which implies a simple relation exists between a complex-
%   valued member of K and its (ii,jj) position in the matrix.  
%  
%   [BOOL,DK]=ISWAVEGRID(K) optionally returns the grid spacing DK.
%  
%   See also SUB2K, K2SUB
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003, 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if strcmp(k, '--t')
  iswavegrid_test,return
end
  
tol=1e-10;
bool=0;

[M,N]=size(k);
if M==N
  dk=abs(k(2)-k(1));
  [ii,jj]=k2sub(k,dk,N);
  [ii2,jj2]=ind2sub([N N],(1:length(k(:)))');
  bool=all(ii2(:)==ii(:)) && all(jj2(:)==jj(:));
end

if ~bool,dk=nan;end

function[]=iswavegrid_test
  
dk=1;N=5;
reporttest('ISWAVEGRID for 5 x 5 dk=1 grid', iswavegrid(wavegrid(dk,N)))

dk=pi;N=5;
reporttest('ISWAVEGRID for 5 x 5 dk=pi grid', iswavegrid(wavegrid(dk,N)))

