function[rp,rm]=quadform(a,b,c)
%QUADFORM  Implements the quadratic formula.
%
%   [RP,RM]=QUADFORM(A,B,C) implements the quadratic formula to
%   solve the equation A*X^2+B*X+C=0   
%
%       RP = 1/2 [- B + SQRT (B^2 - 4*A*C)]
%       RM = 1/2 [- B - SQRT (B^2 - 4*A*C)]
%
%   Only real-valued roots are returned; complex-valued roots are 
%   converted into NANs. RP and RM are thus both real-valued.
%
%   A, B, and C may be either arrays of the same size, or scalars,
%   or a combination.  RP and RM have the same size as the input.
%
%   'quadform --t' runs a test 
%   _______________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details    


if strcmp(a, '--t')
  quadform_test,return
end

index=find(b.^2<4*a.*c);

rp=frac(1,2*a).*(-b+sqrt(b.^2-4*a.*c));
rm=frac(1,2*a).*(-b-sqrt(b.^2-4*a.*c));

if ~isempty(index)
    rp(index)=nan;
    rm(index)=nan;
    rp=real(rp);
    rm=real(rm);
end
%a.*rp.^2+b.*rp+c
%a.*rm.^2+b.*rm+c

function[]=quadform_test
N=1000;
a=rand(N,1);
b=randn(N,1);
c=randn(N,1);
index=find(b.^2-4*a.*c>0);
vindex(a,b,c,index,1);

[rp,rm]=quadform(a,b,c);
b1=a.*rp.^2 + b.*rp+c;
b2=a.*rm.^2 + b.*rm+c;

tol=1e-10;
reporttest('QUADFORM verify plus and minus roots',allall(abs(b1.*b2)<=tol))




