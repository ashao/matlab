function[num]=yf2num(yf)
%YF2NUM  Convert date in 'year.fraction' format to 'datenum' format.
%
%   YF2NUM(YF) where YF is an array of dates in 'year.fraction' format
%   returns the array in Matlab's 'datenum' format.
%
%   See also YEARFRAC, DATENUM, DATEVEC.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2009 J.M. Lilly --- type 'help jlab_license' for details  
    
if strcmp(yf, '--t')
  yf2num_test,return
end

y=floor(yf);

%Number of days in this year?
d0=datenum(y-1,12,31);
d1=datenum(y,12,31);
nd=d1-d0;
num=d0+nd.*(yf-y)+1;

function[]=yf2num_test

yearf=(1850:(1/360):2004)';
num=yf2num(yearf);
yearf2=yearfrac(num);
bool=maxmax(abs(yearf-yearf2))<1e-10;
reporttest('YF2NUM and YEARFRAC, daily resolution, 1e-10 cutoff',bool)

yearf=(1990:(1/360/24):2004)';
%yearf=(1850:(1/360/24):2004)';
num=yf2num(yearf);
yearf2=yearfrac(num);
bool=maxmax(abs(yearf-yearf2))<1e-10;
reporttest('YF2NUM and YEARFRAC, hourly resolution, 1e-10 cutoff',bool)
