function[x]=deblankstrs(x)
%DEBLANKSTRS Deblanks all elements in a cell array of strings
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
for i=1:length(x)
  x{i}=deblank(x{i});
end
