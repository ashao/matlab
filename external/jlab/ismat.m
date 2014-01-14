function[b]=ismat(x)
%ISMAT   Test whether the argument is a 2-D matrix; false for scalars.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details
b=(size(x,1)>1).*(size(x,2)>1).*(ndims(x)==2);
