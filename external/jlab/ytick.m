function[]=ytick(i1,i2)
%YTICK  Sets locations of y-axis tick marks.
%
%   YTICK(DY) or YTICK DY where DY is a number uses DY as the interval
%   between successive y-tick marks, with the endpoints being the
%   axes limits, and applies this to the current axis.
%
%   YTICK(Y) where Y is an array sets the 'ytick' property of the
%   current ayis to Y.
%
%   YTICK(H,DY) applies the change to axes H.
% 
%   See also XTICK, ZTICK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details   
    

if strcmp(i1,'--t')
    return
end

if(nargin==1)
   h=gca;
   dy=i1;
else
   h=i1;
   dy=i2;
end

if ischar(dy)
  dy=str2double(dy);
end

if length(dy)==1
  y=get(h,'ylim');
  set(gca,'ytick',(y(1):dy:y(2)))
else
  set(gca,'ytick',dy);
end
