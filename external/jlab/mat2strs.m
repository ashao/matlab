function[y]=mat2strs(x)
%MAT2STRS  Converts a string matrix into a cell array of strings.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 1999, 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
for i=1:size(x,1)
  y{i}=deblank(x(i,:));
end

y=packstrs(y);
