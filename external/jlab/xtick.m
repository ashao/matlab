function[]=xtick(i1,i2)
%XTICK  Sets locations of x-axis tick marks.
%
%   XTICK(DX) or XTICK DX where DX is a number uses DX as the interval
%   between successive x-tick marks, with the endpoints being the
%   axes limits, and applies this to the current axis.
%
%   XTICK(X) where X is an array sets the 'xtick' property of the
%   current axis to X.
%
%   XTICK(H,DX) applies the change to axes H. 
% 
%   See also YTICK, ZTICK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details  

if strcmp(i1,'--t')
    return
end


if(nargin==1)
   h=gca;
   dx=i1;
else
   h=i1;
   dx=i2;
end


if ischar(dx)
  dx=str2double(dx);
end

if length(dx)==1
  x=get(h,'xlim');
  set(gca,'xtick',(x(1):dx:x(2)))
else
  set(gca,'xtick',dx);
end
