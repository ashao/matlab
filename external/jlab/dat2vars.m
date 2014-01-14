function[varargout]=dat2vars(data)
%DAT2VARS  Put the columns of a matrix into named vectors.
%
%   [O1,O2,O3,...]=DAT2VARS(DATA) puts the nth column of DATA into the
%   Nth output argument.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000--2007 J.M. Lilly --- type 'help jlab_license' for details
  

for i=1:size(data,2) 
       varargout{i}=data(:,i);
end




