function[b]=isset(x)
%ISSET  Tests whether an input array is a set.
%  
%   ISSET(X) returns true if no elements of X are repeated, or 
%   if X is empty or a scalar, and zero otherwise.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2003, 2004 J.M. Lilly --- type 'help jlab_license' for details
  
x=x(:);  
if ~isempty(x)
    
  x=row2col(x);
  x2=osum(x,-x);
  
  %Outer sum should have all zeros along the diagonal and
  %nowhere else.  
  
  b=all(diag(x2)==0) && (length(find(x2==0))==length(x));

else
  %True for the empty set.
  b=1;
end
