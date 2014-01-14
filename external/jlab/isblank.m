function[b]=isblank(x)
%ISBLANK  Tests whether elements of a string are blanks.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
 
x=real(x);
blk=real(' ');
b=(x==blk);
