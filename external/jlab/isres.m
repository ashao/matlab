function[b,m]=isres(s1,s2,k1,k2,tol)
%ISRES  Test whether input wavenumbers form a resonant triad.
%
%   ISRES(S1,S2,K1,K2) where S1 and S2 are plus or minus one and K1
%   and K2 are arrays of identical size, returns a matrix of size K1
%   will elements equal to one if S1*OM(K1)+S2*OM(K2)=OM(S1*K1+S2*K2)
%   for the corresponding elements of K1 and K2, and zero otherwise.
%  
%   ISRES returns zero for elements for at which either K1 or K2, or
%   both, are zero.  ISRES returns NAN for elements for at which 
%   either K1 or K2, or both, are NAN.
%   
%   ISRES(S1,S2,K1,K2,TOL) uses numerical tolerance TOL, which has a 
%   default value of 10^-6.
%  
%   [B,M]=ISRES(S1,S2,K1,K2) also returns matrix M, the absolute
%   value of the difference between the two sides of the resonance
%   condition.
%
%   See also TRIADRES, VTRIADRES.
%
%   Usage: isres(s1,s2,k1,k2);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

if nargin<5  
   tol=1e-6;  
end
om1=om(k1);
om2=om(k2);
om3=om(s1*k1+s2*k2);
m=abs(s1*om1+s2*om2-om3);
b=(m<tol).*(om1~=0).*(om2~=0);

index=find(isnan(m));
if ~isempty(index)
  b(index)=nan;
end

  
  
