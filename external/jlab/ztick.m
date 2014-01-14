function[]=ztick(i1,i2)
%ZTICK  Sets locations of z-axis tick marks.
%
%   ZTICK(DZ) or ZTICK DX where DX is a number uses DZ as the interval
%   between successive z-tick marks, with the endpoints being the
%   axes limits, and applies this to the current axis.
%
%   ZTICK(Z) where Z is an array sets the 'ztick' property of the
%   current axis to Z.
%
%   ZTICK(H,DZ) applies the change to axes H. 
% 
%   See also XTICK, YTICK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details  

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
  x=get(h,'zlim');
  set(gca,'ztick',(x(1):dx:x(2)))
else
  set(gca,'ztick',dx);
end
