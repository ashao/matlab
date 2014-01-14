function[]=unitplot(x,str)
%UNITPLOT  Plot with automatic axis  
%
%   UNITPLOT
%
%   'unitplot --t' runs a test.
%
%   Usage: []=unitplot();
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(x, '--t')
    unitplot_test,return
end
 
t=(0:size(x,1)-1);
t=t./(size(x,1)-1);
%t=t-mean(t);

if nargin ==2
    plot(t,x,str);
else
    plot(t,x);
end


function[]=unitplot_test
 
%reporttest('UNITPLOT',aresame())
