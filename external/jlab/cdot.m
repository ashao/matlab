function[x]=cdot(a,b)
% CDOT  X=CDOT(A,B) <==> X=REAL(A.*CONJ(B));
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details
  
x=real(a.*conj(b));
