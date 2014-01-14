function[]=jimage(x,y,z)
%JIMAGE  Image plot with scaled cdatamapping and normal ydir.
%  
%   JIMAGE(Z) makes an IMAGE plot of matrix Z with 'cdatamapping' set
%   to 'scaled' and 'ydir' set to 'normal'.  Thus a matrix will appear
%   right-side up, and with each entry being colored according to its
%   position between the minimum and the maximum values of Z. 
%  
%   Unlike PCOLOR, which either throws away the last row and column of
%   your matrix ('faceted' and 'flat' shading) or interpolates between
%   elements ('interp' shading), IMAGE actually shows you the relative
%   values of each matrix element.  
%
%   JIMAGE(X,Y,Z) with vector X and Y uses these as the X- and Y-axes,
%   respectively.  Unlike IMAGE, JIMAGE is picky that the lengths of X
%   and Y match the size of Z.
%
%   For visibility, JIMAGE also sets 'grid off' and 'tickdir' to 'out'.
%
%   Example:  jimage(peaks)  
%  
%   Usage: jimage(z)
%          jimage(x,y,z)
%  
%   See also IMAGE.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(x,'--t')
    return
end

if nargin==1
  image(x,'cdatamapping','scaled')
else
  if length(x)~=size(z,2)
    error('Length of X must match number of columns of Z')
  end
  if length(y)~=size(z,1)
    error('Length of Y must match number of rows of Z')
  end
  image(x,y,z,'cdatamapping','scaled')
end

set(gca,'ydir','normal')
outticks
grid off


