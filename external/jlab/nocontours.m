function[h]=nocontours(h)
%NOCONTOURS  Removes contours from a CONTOURF plot.
%
%   NOCONTOURS(H) where H is an array of handles to lines, output by
%   CONTOURF, sets the linestyle for all the lines to 'none'.  This
%   removes the contours from the CONTOURF plot.
%
%   NOCONTOURS with no input arguments sets the linestyle of patch
%   objects in the current axes to 'none'.
%
%   H=NOCONTOURS also returns the handles to the patch objects.  
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details        

  
if nargin==0
  [h,bline,bpatch]=linehandles(gca);
  h=h(bpatch);
  %h=patchhandles(gca);
end

set(h,'linestyle','none')
  
if nargout ==0
  clear h
end

