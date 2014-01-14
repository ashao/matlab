function[bool]=iseven(x)
%ISEVEN Tests whether the elements of an array are even
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details
  
bool=floor(x/2)==x/2;
