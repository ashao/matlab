function[index]=periodindex(varargin)
%PERIODINDEX  Returns time index in increments of instantaneous period.
%
%   INDEX=PERIODINDEX(DT,OMEGA,N) returns a time index INDEX that skips 
%   every N times the instantaneous period 2*pi/OMEGA.
%
%   OMEGA is a column vector of instantaneous frequency in radians per unit 
%   time as computed by INSTFREQ. DT is the sample time, a scalar.  The
%   units of 2*pi/OMEGA should be the same as the units of DT.
%
%   Then PERIODINDEX constructs an index into a vector of the same length
%   as OMEGA that skips every N times the instantaneous period 2*pi/OMEGA.
%   Thus N=1 returns an index with one sample per period, etc.
%
%   The first and last N periods are omitted from the index.  Specifically,
%   the first value of INDEX is N times first value of 2*pi/OMEGA, and 
%   the last value of INDEX is not larger than the length of OMEGA miums
%   N times the last value of 2*pi/OMEGA.
%
%   OMEGA can have interior NANS, in which case the first and last N 
%   periods of each contiguous block of non-NAN data are omitted as above.
%
%   INDEX=PERIODINDEX(OMEGA,N) also works, in which case DT is taken to be
%   as unity.  Thus OMEGA must have units of radians per sample interval.
%
%   PERIODINDEX is useful with ELLIPSEPLOT for plotting ellipses a
%   specified number of periods apart.
%
%   See also ELLIPSEPLOT.
%
%   Usage: index=periodindex(dt,omega,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 


if nargin==3
    dt=varargin{1};
    varargin=varargin(2:3);
else
    dt=1;
end

om=varargin{1};
N=varargin{2};
index=[];

if ~isempty(om)
  if anyany(~isnan(om))
    skip=round(N*frac(2*pi,om.*dt));
     
    [L,ia,ib]=blocklen(~isnan(skip));
    %alternating ones and zeros
    % ia,ib,skip(ia),skip(ib)
    for i=1:length(ia)
        if ~isnan(skip(ia(i)))
            if ia(i)+skip(ia(i))<= length(skip)
                first=skip(ia(i)+skip(ia(i)));
                if ~isnan(first)
                    index(end+1)=ia(i)+skip(ia(i));  
                    while index(end)+skip(index(end))<ib(i)
                        index(end+1)=index(end)+skip(index(end));
                    end
                    if index(end)>ib(i)-skip(ib(i))
                        index=index(1:end-1);
                    end
                end
            end
        end
    end
 end
end
 
%reporttest('PERIODINDEX',aresame())
