function[]=line_color(x,y,c,clim,w,cmap,interval)
% line_color(x,y,c,color_limits,linewidth,cmap)
% 
% Default:  cmap = jet(256)
%           linewidth = 2
%           color_limits = [min(c), max(c)]
%
% Luc Rainville 6/5/08

h=gca;

if nargin<6
  cmap = jet(256);
end

if nargin<5
  w=2;
end

if nargin<4
  clim(1) = min(c);
  clim(2) = max(c);
end

if nargin < 7
    interval = 1;
end
h=gca;
%cmap=colormap;
colormap(cmap);

N=size(cmap,1);
C = round((0.5*(c(1:end-1)+c(2:end))-clim(1))/diff(clim)*(N-1))+1;
C(C<1)=1;
C(C>N)=N;

flag=0;
for n=1:interval:length(x)-1
  if isfinite(C(n))
    plotm(x(n:n+1), y(n:n+1), 'color',cmap(C(n),:),'linewidth',w)
    flag=1;
  end
  if flag==1
      hold on
  end
end
caxis(clim);




