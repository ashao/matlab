function[varargout]=vempty(varargin)
%VEMPTY   Initializes multiple variables to empty sets.
%
%   [X1,X2, ... XN]=VEMPTY  is equivalent to
%		
%      X1=[]; X2=[]; .... XN=[];
%
%   thus initializing all the output variables to empty sets.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information (C)
%   2004 J.M. Lilly --- type 'help jlab_license' for details

  
if nargin~=0
  if strcmp(varargin{1}, '--t')
   vempty_test,return
  end
end


for i=1:nargout
  varargout{i}=[];
end

function[]=vempty_test

z=[];
[x,y]=vempty;
reporttest('VEMPTY', all(x==z&y==z))
  
