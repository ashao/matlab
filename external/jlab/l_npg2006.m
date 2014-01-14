%
%L_NPG2006  Script for loading float used in Lilly and Gascard (2006).
%
%   Try 'type l_npg2006'.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details


%The "pomme" dataset is proprietary; it is not distributed here.
load pomme
use pomme

jj=11; 
vindex(id,num,p,t,lat,lon,cv,cx,p,jj,2);

t_insitu=t;
t=sw_ptmp(35+0*t_insitu,t_insitu,p,0*p);

[latc,lonc]=vmean(lat,lon,1);
cx=latlon2xy(lat,lon,latc,lonc);

ii=(1:1120);
vindex(id,num,p,t,lat,lon,cv,cx,ii,1);

num = num-datenum(2001,1,1)+1; 

dt=1/6;  %Time step in days
 
description='NPG2006 float dataset for Lilly and Gascard (2006)';
link='http://www.nonlin-processes-geophys.net/13/467/2006/npg-13-467-2006.html';
creator='Created by the L_NPG2006 script by J. M. Lilly';
timestamp=datestr(now);

matsave -v6 npg2006 description link creator timestamp num dt lat lon p t cx cv 


