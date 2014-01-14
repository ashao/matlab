function[h]=fastcontour(varargin)
%FASTCONTOUR  Lightning-fast "fake" contouring for large matrices.
%
%   FASTCONTOUR(X,Y,Z,C) contours matrix Z at levels C, where X and Y
%   are the horizontal (row) and vertical (column) coordinates for Z.
%
%   FASTCONTOUR draws a "fake" contour by drawing many small points
%   instead of a closed curve.  For large matrices with complicated 
%   contours, this approach can be vastly faster than a true contour.
%
%   FASTCONTOUR(X,Y,Z), FASTCONTOUR(Z,C), and FASTCONTOUR(Z) also work.
% 
%   Unlike CONTOUR, if the contour levels C are not chosen, the single
%   level C=0 is used.
%
%   FASTCONTOUR(..., STR) specifies the style to use.  STY may be a 
%   string, e.g. STR='b.' (the default),  or a cell array of strings, 
%   e.g. STR{1}='b.', STR{2}='r*'.  The styles are cycled through if
%   more contour levels than sytles are specified.
%
%   For small matrices, avoid '.' as this marker will not show up.  
%
%   H=FASTCONTOUR(...) also returns the line handles to the "contours".
%
%   FASTCONTOUR is particularly handy when drawing coastlines.
%
%   Usage:  h=fastcontour(z);
%           h=fastcontour(x,y,z);
%           h=fastcontour(z,c);
%           h=fastcontour(x,y,z,c);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details

%[x,y,z]=peaks;
%fastcontouf(x,y,z,'bo')

c=0;
str{1}='b.';

na=nargin;
if iscell(varargin{end}) || ischar(varargin{end})
    if iscell(varargin{end})
        str=varargin{end};
    elseif ischar(varargin{end})
        str{1}=varargin{end};
    end
    varargin=varargin(1:end-1);
    na=na-1;
end

if na==1
    z=varargin{1};
elseif na==2
    z=varargin{1};
    c=varargin{2};
elseif na==3
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
elseif na==4
    x=varargin{1};
    y=varargin{2};
    z=varargin{3};
    c=varargin{4};
end


if na==1 || na==2
    x=(1:size(z,2));
    y=(1:size(z,1));
end


if (~ismat(x) && ismat(y)) || (ismat(x) && ~ismat(y))
    error('X and Y must either be both matrices or both vectors');
end

if ~ismat(x)
   [x,y]=meshgrid(x,y);
end

%Check hold state
if ishold
    washold=1;
else 
    cla
    washold=0;
end

hold on
for i=1:length(c)
    index=fastcontour1(x,y,z-c(i));
    if ~isempty(index)     
        nstr=mod(i-1,length(str))+1;  %Cycle through strings
        h(i)=plot(x(index),y(index),str{nstr},'markersize',5);
    else
        h(i)=nan;
    end
end

%Return to hold state
if ~washold 
   hold off
end

if nargout==0
    clear h
end

function[index]=fastcontour1(x,y,z)
bool1=sign(vshift(z,1,1).*vshift(z,-1,1))<0;
bool2=sign(vshift(z,1,2).*vshift(z,-1,2))<0;

%Remove border points
bool=bool1|bool2;
bool([1 end],:)=0;
bool(:,[1 end])=0;
%spy(bool)
index=find(bool);



