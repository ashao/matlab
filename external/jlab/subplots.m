function[h]=subplots(m,n)
%SUBPLOTS  Return handles to subplots of the current figure.
%
%   H=SUBPLOTS(M,N) returns a vector H containing all handles to the
%   subplots of the current figure, where M is the number of rows and
%   N is the number of columns.
%
%   AXES(H(i)) is then equivalent to SUBPLOT(M,N,i).  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details      

  
for i=1:n
  for j=1:m
     index=(j-1)*n+i;
     h(index)=subplot(m,n,index);
  end
end  
