function[varargout]=nnsd(varargin)
%NNSD  Number of nonsingleton dimensions
%  
%   N=NNSD(X) returns the number of nonsingleton dimensions of X.
%   Unlike Matlab's NDIMS, which thinks that an array and a scalar
%   both have dimension 2, NNSD defines the dimension of a scalar to
%   be zero and that of an array to be one.  Singleton dimensions are
%   never counted.
% 
%   [N1,N2,...NM]=NNSD(X1,X2,...XM) returns the dimensions of multiple 
%   input arguments.  If zero or one output arguments are given, a
%   single row array [N1 N2 ...NM] is output.  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002,2004 J.M. Lilly --- type 'help jlab_license' for details        
  
  
if strcmp(varargin{1}, '--t')
  nnsd_test,return
end


for i=1:nargin
  x=varargin{i};
  sx=size(x);
  nd(i)=length(sx)-length(find(sx==1));
  varargout{i}=nd(i);
end

if nargout==0 || nargout==1
   varargout{1}=nd;
end

function[]=nnsd_test

x=nnsd(1,(1:10),[ (1:10)' (1:10)']);
bool(1)=aresame(x,[0 1 2]);
%disp('Should be 0 1 2')

x=(1:10)';
z(:,:,3)=x;
q(:,:,:,1)=z;
q(:,:,:,2)=z;

x=nnsd(x,z,q);
bool(2)=aresame(x,[1 2 3]);
%disp('Should be 1 2 3')
reporttest('NNSD',all(bool))
