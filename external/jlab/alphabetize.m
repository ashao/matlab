function[str]=alphabetize(str)
% ALPHABETIZE Sorts a string matrix by its first column
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details
firstcol=real(str(:,1));
[firstcol,index]=sort(firstcol);
str=str(index,:);
