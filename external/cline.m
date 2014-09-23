% This function plots a 3D line (x,y,z) encoded with scalar color data (c)
% using the specified colormap (default=jet);
%
% SYNTAX: h=cline(x,y,z,c,colormap);
%
% DBE 09/03/02

function h=cline(x,y,c,cmap,cmin,cmax);

if nargin==0  % Generate sample data...
  x=linspace(-10,10,101);
  y=2*x.^2+3;
  c=exp(z);
  w=z-min(z)+1;
  cmap='jet';
elseif nargin<3
  fprintf('Insufficient input arguments\n');
  return;
elseif nargin==3
  cmap='jet';
end

if ischar(colormap);
    cmap=colormap(cmap); 
else
    cmap = colormap;
end% Set colormap

if nargin < 5    
    yy=linspace(min(c),max(c),size(cmap,1));
else
    yy=linspace(cmin,cmax,size(cmap,1));
end% Generate range of color indices that map to cmap
cm = spline(yy,cmap',c);                  % Find interpolated colorvalue
cm(cm>1)=1;                               % Sometimes iterpolation gives values that are out of [0,1] range...
cm(cm<0)=0;
nanidx = sum(isnan(cm))>0;

% Plot line segment with appropriate color for each data pair...
  for i=(length(x)-1):-1:1
    if sum( isnan(cm(:,i)) ) == 0
        h(i)=linem([x(i) x(i+1)],[y(i) y(i+1)],'color',[cm(:,i)],'LineWidth',2);
    end
  end

return