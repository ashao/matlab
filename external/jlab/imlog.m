function[y]=imlog(x)
%IMLOG   IMLOG(X)=UNWRAP(IMAG(LOG(X)))
%
%   IMLOG returns the imaginary part of the natural log of its 
%   argument, as in IMLOG(a exp(i phi)) = phi.  The output is 
%   unwrapped using UNWRAP so that it varies continuously.
%
%   See also RELOG.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2006--2007 J.M. Lilly --- type 'help jlab_license' for details  

warning('off','MATLAB:log:logOfZero')
y=unwrap(imag(log(x)));
warning('on','MATLAB:log:logOfZero')
