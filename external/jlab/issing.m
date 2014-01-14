function[b,n]=issing(x)
%ISSING Tests whether the argument is a singleton array.
%
%   ISSING(X) is true (=1) if exactly one dimension of array X has
%   greater than one element, and is zero otherwise.  Such an array X
%   is called a singleton array or a vector.
%
%   ISSING(X) is defined to be false if X is a scalar, that is if all
%   dimensions of X have only one element.
%
%   [B,N]=ISSING(X) returns the result of the singleton array test in
%   boolean array B and also returns the number N of the dimension of
%   X having the greatest length.  Thus N=1 for a column vector, N=2
%   for a row vector, etc.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information (C)
%   2000, 2004 J.M. Lilly --- type 'help jlab_license' for details

szx=size(x);
b=(length(find(szx>1))==1);
if nargout ==2
   [maxlen,n]=max(szx);
end


