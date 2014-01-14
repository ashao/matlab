function[mf]=monthfrac(num)
%MONTHFRAC  Convert date from 'datenum' format to 'month.fraction'.
%  
%   MF=MONTHFRAC(NUM) where NUM is an array of dates in Matlab's 'datenum'
%   format, returns the fraction of the month at each date.
%
%   That is, for a date at noon on the second day of April, a month having 
%   thirty days, MONTHFRAC will return 4+ 1.5/30 = 4.0500.
%
%   The actual number of days in each month is used, that is, accounting 
%   for the fact that a month could have 28, 29, 30, or 31 days.
%  
%   See also YEARFRAC, DATENUM, DATEVEC
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2009 J.M. Lilly --- type 'help jlab_license' for details        

bool=find(isnan(num));
num(bool)=0;

[y,mo,d,h,mi,s] = datevec(num);

%Number of days in this month
nd=datenum(y,mo+1,1)-datenum(y,mo,1,0,0,0);

%Day we're at in this month (Full days completed, one less than today's date)
na=datenum(y,mo,d,h,mi,s)-datenum(y,mo,1,0,0,0);

mf=mo+na./nd;
mf(bool)=nan;


  
  
