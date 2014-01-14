
function[]=monthlyticks
%MONTHLYTICKS   Set x-axis appropriate for months
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        

set(gca,'xtick',(1:12),'xlim',[0.5 12.5],'xticklabel',('JFMAMJJASOND')')
