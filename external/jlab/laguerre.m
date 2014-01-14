function[y]=laguerre(x,k,c)
%LAGUERRE Generalized Laguerre polynomials
%
%   Y=LAGUERRE(X,K,C) where X is a column vector returns the
%   generalized Laguerre polynomials specified by parameters K and C.
%  
%   LAGUERRE is used in the computation of the generalized Morse
%   wavelets and uses the expression given by Olhede and Walden (2002),
%  "Generalized Morse Wavelets", Section III D. 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2009 F. Rekibi and J. M. Lilly 
%                         --- type 'help jlab_license' for details  


y=0*x;
for m=0:k
   %Log of gamma function much better ... trick from Maltab's ``beta'' 
   %fact=gamma(k+c+1)./(gamma(c+m+1).*gamma(k-m+1));  
   fact=exp(gammaln(k+c+1)-gammaln(c+m+1)-gammaln(k-m+1));
   y=y+((-1).^m).*fact.*(x.^m)./gamma(m+1);
end
