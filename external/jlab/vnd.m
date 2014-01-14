function[varargout]=vnd(varargin)
% VND  Number of dimensions
%  
%   N=VND(X) returns the number of dimensions of X.  Unlike matlab,
%   which thinks that an array and a scalar both have dimension 2, VND
%   defines the dimension of a scalar to be zero and that of an array to
%   be one.  Singleton dimensions are never counted.  
% 
%   [N1,N2,...NM]=VND(X1,X2,...XM) returns the dimensions of multiple 
%   input arguments.  If zero or one output arguments are given, a
%   single row array [N1 N2 ...NM] is output.  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
    
if strcmp(varargin{1}, '--t')
  vnd_test,return
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

function[]=vnd_test
reporttest('VND', all(vnd(1,(1:10),[ (1:10)' (1:10)'])==[0 1 2]))

