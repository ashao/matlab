function[y]=frac(x1,x2)
%FRAC   FRAC(A,B)=A./B;
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details  
warning('off','MATLAB:divideByZero')
y=x1./x2;
warning('on','MATLAB:divideByZero')
