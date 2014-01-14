function[y]=relog(x)
%RELOG   RELOG(X)=REAL(LOG(X))
%
%   RELOG returns the real part of the natural log of its argument, 
%   as in RELOG(a exp(i phi)) = a. 
%
%   See also IMLOG.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2006--2007 J.M. Lilly --- type 'help jlab_license' for details  

warning('off','MATLAB:log:logOfZero')
y=real(log(x));
warning('on','MATLAB:log:logOfZero')