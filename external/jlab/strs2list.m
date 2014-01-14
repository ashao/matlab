function[row]=strs2list(x)
%STRS2LIST  Converts a cell array of strings into a comma-delimited list.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003--2006 J.M. Lilly --- type 'help jlab_license' for details      
  
if ~iscell(x)
  xtemp=x;
  clear x
  x{1}=xtemp;
end

%x=packstrs(x);
  
  
M=length(x);
for i=1:M
    n(i)=length(x{i});
end
N=max(n);

row=[];

for i=1:M
   row=[row,x{i},',']; 
end

row=row(1:end-1);
