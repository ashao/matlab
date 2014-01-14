function[x]=packstrs(x)
%PACKSTRS Removes empty entries from a cell array of strings
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2008 J.M. Lilly --- type 'help jlab_license' for details        
  
x=deblankstrs(x);
n=cellength(x);
x=x(n>=1);
