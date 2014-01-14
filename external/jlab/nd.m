function[varargout]=nd(varargin)
% ND  Number of the last nonsingleton dimension.
%  
%   N=ND(X) returns the number of the last nonsingleton dimension of X.
%
%   This provides a useful definition of the dimensionality of X.
%   Unlike Matlab's NDIMS, which thinks that a column vector and a
%   scalar both have dimension 2, ND defines the dimension of a scalar
%   to be zero and that of a column vector to be one.  A row vector,
%   however, has ND equal to 2.
%  
%   Note that ND(X) changes if X is permuted. See NNSD for a different
%   definition of dimensionality that is independent of permutation.
%  
%   [N1,N2,...NM]=ND(X1,X2,...XM) returns the dimensions of multiple
%   input arguments.  If zero or one output arguments are given, a
%   single row array [N1 N2 ...NM] is output.
%
%   See also NNSD.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2010 J.M. Lilly --- type 'help jlab_license' for details        
  
if strcmp(varargin{1}, '--t')
  nd_test,return
end


for i=1:nargin
  x=varargin{i};
  sx=size(x);
  if isempty(x)
    nd1=nan;
  elseif isscalar(x)
    nd1=0;
  else
   nd1=max(find(sx~=1));
  end
  
%  lsx=ndims(x);
%   if jiscol(x)
%     lsx=1;
%   elseif isscalar(x)
%     lsx=0;
%   end

  nd(i)=nd1;
  varargout{i}=nd(i);
end

if nargout==0 || nargout==1
   varargout{1}=nd;
end

function[]=nd_test

x=nd([],1,(1:10),[ (1:10)' (1:10)']);
bool(1)=aresame(x,[nan 0 2 2]);
%disp('Should be NAN 0 2 2')

x=(1:10)';
z(:,:,3)=x;
q(:,:,:,1)=z;
q(:,:,:,2)=z;
c(1,1,:)=x;  

x=nd(x,z,q,c,permute(c,[3 2 1]));
bool(2)=aresame(x,[1 3 4 3 1]);
%disp('Should be 1 3 4 3 1')
reporttest('ND',all(bool))

