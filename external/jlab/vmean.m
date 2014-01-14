function[varargout] = vmean(varargin)
%VMEAN  Mean over finite elements along a specified dimension.
%
%   Y=VMEAN(X,DIM) takes the mean of all finite elements of X along      
%   dimension DIM. 
%                                                                         
%   [Y,NUM]=VMEAN(X,DIM) also outputs the number of good data points NUM, 
%   which has the same dimension as X.                              
%
%   [Y1,Y2,...YN]=VMEAN(X1,X2,...XN,DIM) also works.
%
%   VMEAN(X1,X2,...XN,DIM);  with no output arguments overwrites the 
%   original input variables.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(varargin{1}, '--t')
  vmean_test,return
end

dim=varargin{end};

for i=1:length(varargin)-1
  [varargout{i},numi{i}]=vsum(varargin{i},dim);

  numi{i}(numi{i}==0)=nan;
  varargout{i}=varargout{i}./numi{i};
end

for i=length(varargin):nargout
  varargout{i}=numi{i-length(varargin)+1};
end

eval(to_overwrite(nargin-1))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]=vmean_test
x1=[1 2 ; nan 4];
x2=[inf 6; nan 5];
ans1=[3/2 4]';
ans2=[6 5]';

vmean(x1,x2,2);
reporttest('VMEAN output overwrite', aresame(x1,ans1) && aresame(x2,ans2))

x1=[1 2 ; nan 4];
ans1=[3/2 4]';
ans2=[2 1]';

[y1,y2]=vmean(x1,2);
reporttest('VMEAN mean & num', aresame(y1,ans1) && aresame(y2,ans2))










