function[b]=ismname(x)
%ISMNAME  Tests whether a string is an m-file name; cells OK  
%
%   ISMNAME STR tests whether a string STR has the ending ".m".
%
%   ISMNAME(X) where X is a cell array of strings returns a boolean
%   array testing whether each string has the ending ".m".
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details  
  
if iscell(x)
  for i=1:length(x)
   b(i)=ismname1(x{i});
  end
else
   b=ismname1(x);
end

function[b]=ismname1(x)
 
x=deblank(x);
if length(x)>=2
    b=strcmp(x(end-1:end),'.m');
else 
    b=false;
end

if strcmp(x(1),'#') || strcmp(x(1),'.')  %Exclude funny names
   b=false;
end

