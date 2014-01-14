function[b]=issquare(x)
%ISSQUARE   Tests whether the argument is a square matrix
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        

b=0;
if nd(x)==2
  if size(x,1)==size(x,2)
     b=1;
  end
end
  
