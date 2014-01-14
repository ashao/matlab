function[bool]=isodd(x)
%ISODD Tests whether the elements of an array are odd
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details
bool=iseven(x+1);
