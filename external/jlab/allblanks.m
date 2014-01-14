function[bool]=allblanks(x)
%ALLBLANKS  Equals one if string argument is all blanks or empty, else 0 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details
  
if isempty(x)
  bool=1;
else  
  bool=all(real(x)=='32');
end
