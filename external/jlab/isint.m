function[b,n]=isint(x)
%ISINT Tests if an array is composed of integers.
%                                 
%   B=ISINT(X) returns a boolean array B of the same size as X, with 
%   each element equal to one if the corresponding element of X is an
%   integer, and zero otherwise.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details

resid=x-floor(real(x));
b=abs(resid)==0;



