function[ds1,ds2,ds3]=aquasal(sst,n)
%AQUASAL  Aquarius salinity change with brightness temperature.
%
%   [DS1,DS2,DS3]=AQUASAL(SST) returns the salinity change DS for a 
%   differential brightness temperature change for each of the three 
%   radiometer beams, given a sea surface temperature array SST.
%
%   DS=AQUASAL(SST,N) returns the salinty change for beam number N.
%
%   See also AQUAPRINT, AQUAPLOT.
%
%   'aquasal --f' makes a sample figure.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2007 J.M. Lilly --- type 'help jlab_license' for details    

% if strcmp(sst, '--t')
%   aquasal_test,return
% end
if strcmp(sst, '--f')
  aquasal_fig,return
end
    
ds1=frac(1,0.0176*sst+0.2325);
ds2=frac(1,0.0186*sst+0.2442);
ds3=frac(1,0.0198*sst+0.2569);

if nargin==2
    if n==2
        ds1=ds2;
    elseif n==3
        ds1=ds3;
    end
    clear ds2 ds3
end

%function[]=aquasal_test


function[]=aquasal_fig
figure
t=(0:.2:30)';
[ds1,ds2,ds3]=aquasal(t);
plot(t,[1./ds1 1./ds2 1./ds3]),hold on,linestyle g r b
title('Brightness temperature change per salinity change')
xlabel('Sea surface temperature (^\circ C)')
plot(15.79,0.510,'g*')
plot(15.79,0.538,'r*')
plot(15.79,0.569,'b*')
