function[str]=strscat(x)
%STRSCAT  Concatenates a cell array of strings into one long string
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details      
  
str=[];  
for i=1:length(x)
  str=[str x{i}];
end
