function[varargout]=vpower(varargin)
%VPOWER  Raises array to the specified power.
%
%   Y=VPOWER(X,N) is equivalent to Y=X.^N;
%
%   [Y1,Y2,...YN]=VPOWER(X1,X2,...XN,POWER,N) also works.
%
%   VPOWER(X1,X2,...XN,POWER,N); with no output arguments overwrites
%   the original input variables.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmp(varargin{1}, '--t')
  vpower_test,return
end

n=varargin{end};

for i=1:length(varargin)-1
  varargout{i}=varargin{i}.^n;
end

eval(to_overwrite(nargin-1))

    
function[]=vpower_test
x1=[1 2; 3 4];
x2=x1;
ans1=[1 4 ; 9 16];

vpower(x1,x2,2);
reporttest('VPOWER', aresame(x1,ans1) && aresame(x2,ans1))
