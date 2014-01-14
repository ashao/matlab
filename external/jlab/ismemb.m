function[bool]=ismemb(a,b,tol)
%   ISMEMB Tests whether the elements of an array are members of a set.
%
%   ISMEMB(A,B) for scalar A is true if A is a member of the array B,
%   and false otherwise.  For array A, ISMEMB tests whether each
%   element of A is a member of B, and returns an array of size
%   SIZE(A).
%
%   Equality is judged to within a small numerical tolerance specified
%   by ISMEMB(A,B,TOL), with a default of TOL=1e-10;
%
%   See also SUBSET, LOOKUP
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2003--2008 J.M. Lilly --- type 'help jlab_license' for details
  
if nargin==2  
  tol=1e-10;  
end

if any(~isfinite(a)) || any(~isfinite(b))
  error('Elements of A and B should be finite (non-NAN, non-INF).')
end

bool=true(size(a));

a=a(:);
b=b(:);

%bool1=(osum(a,-b)==0);
bool1=(abs(osum(a,-b))<=tol);
bool(sum(bool1,2)==0)=0;  %this is equivalent to "any"

%This convienently gets around explicitly needing to know the 
%dimensionality of A, and also reshaping.  
  
  
