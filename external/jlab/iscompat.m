function[bool]=iscompat(a,b)
%ISCOMPAT  Tests whether an array's size is "compatible" with another's.
%
%   ISCOMPAT(A,B) returns true if the sizes of the dimensions of A and
%   B is the same, up to the dimensionality of B.  The size of A is
%   then said to be "compatible" with the size of B. This permits, for
%   example, B to be a matrix of spatial positions, and A to contain
%   some time-varying  parameter at those positions.
%
%   By 'dimensionality' of B is meant ND(B), not NDIMS(B).
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2003, 2004 J.M. Lilly --- type 'help jlab_license' for details  
    
if strcmp(a, '--t')
  iscompat_test,return
end

sizea=size(a);
sizeb=size(b);
ndb=nd(b);

if ndb>0
  bool=all(sizeb(1:ndb)==sizea(1:ndb));
else
  %B is scalar; special case
  bool=1;
end

function[]=iscompat_test
x=[1 2; 3 4];y(:,:,2)=x;
b(1)=iscompat(y,x);
b(2)=iscompat(y,1);
b(3)=iscompat(y,1:2);
b(4)=iscompat(y,(1:2)');
b(5)=iscompat(x,(1:2)');  %This would fail with the usual NDIMS

bool=aresame(b,[1 1 0 1 1]);
reporttest('ISCOMPAT',bool)

