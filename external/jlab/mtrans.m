function[varargout]=mtrans(varargin)
% MTRANS  Multitaper "eigentransform" computation.
% 
%   [F,W]=MTRANS(X,PSI) returns the multitaper "eigentransform" matrix for 
%   use in multitaper spectral estimates or eigenspectral SVD analysis.                 
%
%       X  --  M x N matrix containing N length M time series
%     PSI  --  M x K matrix of K data tapers
%       W  --  [M/2] x K x  N eigentransfrom matrix (for real X)
%              [M/2] x K x 2N eigentransfrom matrix (for complex X)
%
%   In the above, [M/2] means M/2 if M is even, and (M-1)/2 is M is odd.
%
%   [F,W1,W2,...,WN]=MTRANS(X1,X2,...,XN,PSI) also works, where X1,...,XN
%   are all M x N matrices.
%
%   See also: SLEPTAP, HERMFUN, MSPEC, MSVD.
%
%   Usage:  [f,w]=mtrans(x,psi);  
%           [f,wx,wy]=mtrans(x,y,psi);    
%           [f,wx,wy,wz]=mtrans(x,y,z,psi);    
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2011 J.M. Lilly --- type 'help jlab_license' for details        
  
%Sort out input arguments
psi=varargin{end};
x=varargin(1:end-1);


%Size check
for i=1:length(x)
    sizex(i,:)=size(x{i});
end
if ~allall(sizex(:,1)==sizex(1,1))||~allall(sizex(:,2)==sizex(1,2))
   error('All input arguments should be the same size')
end

sizex=sizex(1,:);
psimat=vrep(psi,sizex(2),3);


f=fourier(sizex(1));

index=1:length(f);
varargout{1}=f;

for i=1:size(x,2)
    mmat=mtrans1(x{i},psimat);
    varargout{i+1}=mmat(index,:,:);
end

function[mmat]=mtrans1(x,psimat)

Nnans=length(find(isnan(x)));
if Nnans>0
    disp(['MTRANS swapping ' int2str(Nnans) ' NANs for zeros.'])
    vswap(x,nan,0);
end
x=permute(x,[1 3 2]);
xmat=vrep(x,size(psimat,2),2);
mmat=fft(psimat.*xmat,[],1);


