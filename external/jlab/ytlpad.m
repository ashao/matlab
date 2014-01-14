function[str]=ytlpad(h,n)
%YTLPAD  Pads the ytick labels with a leading space.
%
%   YTLPAD shifts the ytick labels of the current axis over slightly
%   to the right by inserting a column of spaces on the left-hand
%   side.
%
%   STR=YTLPAD returns the new string array.  
%
%   YTLPAD(H) applies itself to axes with handle H.   
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details 
   
if nargin==0
  h=gca;
end

yt=get(h,'yticklabel');
sp=char(real(' ')*ones(size(yt(:,1))));
yt=[sp yt];
set(h,'yticklabel',yt);
yt=get(h,'yticklabel');

if nargout==1
  str=yt;
end
