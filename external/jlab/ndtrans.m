function[y]=ndtrans(x,n)
%NDTRANS  Generalized transpose of a potentially multidimensional vector.
%
%   Y=NDTRANS(X,N), where X is an array having exactly one
%   non-singleton dimension, 'transposes' X so that the
%   non-singleton dimension is oriented along dimension N.
%	
%   When X is complex-valued, the sign of the imaginary terms in Y
%   will be the same as in X.  NDTRANS is a therefore a conjugate
%   transpose operator, generalized to multiple dimensions.
%
%   Y=NDTRANS(X) returns a column vector for any orientation
%   of the input array X.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details        


L=length(x);
ndim=ndims(x);

if nargin==1
   n=1;
end

if ~issing(x)
   %if more than one dimension is non-singleton, error
   error('Input array must have only one non-singleton dimension.')
else
   str='y=reshape(x,';
   for i=1:n-1
       str=[str,'1,'];
   end
   if n==1
      %Reshape requires at least two size elements.
      str=[str,'L,1);'];
   else
      str=[str,'L);'];
   end   
   eval(str)
end







