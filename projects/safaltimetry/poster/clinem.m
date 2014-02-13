% This function plots a 3D line (x,y,z) encoded with scalar color data (c)
% using the specified colormap (default=jet);
%
% SYNTAX: h=cline(x,y,z,c,colormap);
%
% DBE 09/03/02

function h=clinem(x,y,c,cmap,varargin);

if nargin==0  % Generate sample data...
    x=linspace(-10,10,101);
    y=2*x.^2+3;
    z=sin(0.1*pi*x);
    c=exp(z);
    w=z-min(z)+1;
    cmap='jet';
elseif nargin<3
    fprintf('Insufficient input arguments\n');
    return;
elseif nargin==3
    cmap='jet';
elseif nargin < 5
    varargin = {};
end

cmap=colormap(cmap);                      % Set colormap
yy=linspace(min(c),max(c),size(cmap,1));  % Generate range of color indices that map to cmap
cm = spline(yy,cmap',c);                  % Find interpolated colorvalue
cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;

h = zeros(length(y)-1);
% Lot line segment with appropriate color for each data pair...
for i=1:length(y)-1
    if ~isnan(cm(:,i))
        h(i)=linem([x(i) x(i+1)],[y(i) y(i+1)],'color',[cm(:,i)],varargin{:});
    end
end
    return