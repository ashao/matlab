function[]=usefulstrings
% USEFULSTRINGS Note on numeric representations of useful strings
%
%   9  <==> tab.
%   10  <==> newline with carriage return.
%   13  <==> return to beginning of line
%   138 <==> newline without carriage return.
%   32  <==> space
%
%   Use real(X) to find the numeric representation of a string X.
%   Use char(X) to convect numeric array X into a string.
%
%   e.g. char([116 101 120 116]) ='text';
%        real('text') = [116 101 120 116];
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details    
  
