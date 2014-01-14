function[varargout]=vdiff(varargin)
%VDIFF	Length-preserving first central difference.
%
%   DX=VDIFF(X,DIM) differentiates X along dimension DIM using the first 
%   central difference; DX is the same size as X.                                 
%                                                                        
%   [D1,D2,...,DN]=VDIFF(X1,X2,...,XN,DIM) for multiple input variables 
%   also works. 
%
%   VDIFF(X1,X2,...,DIM); with no output arguments overwrites the
%   original input variables.
%
%   DXDT=VDIFF(DT,...) optionally uses scalar timestep DT to approximate
%   a time derivative, i.e. DXDT equals DX divided by DT.
%   _____________________________________________________________________
%
%   First and last points
%
%   The first and last points must be treated differently, as the central 
%   difference is not defined there.  Three different methods can be used.
%
%   VDIFF(...,STR) specifies which method to use.
%
%        'endpoint'  uses the first forwards / first backwards difference
%                    at the first and last point, respectively.  
%        'periodic'  treats the array as being periodic along dimension DIM,
%                    so that the central difference is defined at endpoints.
%        'nans'      fills in the first and last values with NANs.
%
%   The default behavior is 'endpoint'.
%   _____________________________________________________________________
%
%   'vdiff --t' runs some tests.
%
%   Usage:  x=vdiff(x,dim);
%           x=vdiff(dt,x,dim);
%           x=vdiff(dt,x,dim,'periodic');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2011 J.M. Lilly --- type 'help jlab_license' for details    
  
%   I am so irritated by diff

 
if strcmp(varargin{1}, '--t')
  vdiff_test,return
end


if length(varargin{1})==1
    dt=varargin{1};
    varargin=varargin(2:end);
else
    dt=1;
end

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='endpoints';
end
    
%if length(varargin{end})==1
   n=varargin{end};
   varargin=varargin(1:end-1);
%else
 %  n=1;
%end

for i=1:length(varargin)
  varargout{i}=vdiff1(varargin{i},n,str)./dt;
end

if nargin>1
  eval(to_overwrite(length(varargin)))
end

function[y]=vdiff1(x,n,str)
  	  
if ~isempty(x)
    y=vshift(x,1,n)./2-vshift(x,-1,n)./2;
    %y=vshift(x,1,n)-x;
        
    if strcmp(str(1:3),'end')
        y=vindexinto(y,vindex(x,2,n)-vindex(x,1,n),1,n);
        y=vindexinto(y,vindex(x,size(x,n),n)-vindex(x,size(x,n)-1,n),size(x,n),n);
        %y(1,:)=x(2,:)-x(1,:);
        %y(end,:)=x(end,:)-x(end-1,:);
    elseif strcmp(str(1:3),'nan')
        y=vnan(y,1,n);
        y=vnan(y,size(y,n),n);
    elseif strcmp(str(1:3),'per')
        %Do nothing
    end
else
    y=[];
end
    
function[]=vdiff_test

y1=(1:4)';
y2=2*(1:4)';
[x1,x2]=vdiff(y1,y2,1);
bool=aresame(x1,[1 1 1 1]').*aresame(x2,2*[1 1 1 1]');
reporttest('VDIFF', bool)
vdiff(y1,y2,1);
bool=aresame(y1,[1 1 1 1]').*aresame(y2,2*[1 1 1 1]');
reporttest('VDIFF output overwrite', bool)

dt=pi;
y1=(1:4)';
y2=2*(1:4)';
[x1,x2]=vdiff(pi,y1,y2,1);
bool=aresame(x1,[1 1 1 1]'./dt).*aresame(x2,2*[1 1 1 1]'./dt,1e-10);
reporttest('VDIFF with non-unit time step', bool)

