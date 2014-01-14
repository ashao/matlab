function[y]=sray2strs(x)
%SRAY2STRS  Converts a string array w/ returns into a cell array of strings.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details      

y=[];

index=find(real(x)==10);
  
if ~isempty(index)
  x(index)=',';
  y=list2strs(deblank(x));
end

y=packstrs(y);


