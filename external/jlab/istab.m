function[b]=istab(x)
% ISTAB  Tests whether elements of a string are tab markers.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        
 
  
x=real(x);
blk=9;
b=(x==blk);
