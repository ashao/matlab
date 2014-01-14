function[y]=dk(T,s,bool)
%DK   The Dirichlet kernel
% 
%   DK(T,OM) returns the real-valued Dirichlet kernel.
%
%   DK(T,OM,N) returns either the real-valued Dirichlet kernel (N==0)
%   or the complex-valued (rotated) Dirichlet kernel (N==1).
%
%   T may either be a scalar or an array of the same size as OM.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==2
  bool=0;
end

y=frac(1,2*pi).*frac(sin(s.*T./2),s./2);

index=find(s==0);

if ~isempty(index)
  if length(T)>1
     y(index)=frac(T(index),2*pi);
  else
     y(index)=frac(T,2*pi);
  end
end

if bool
  y=y.*exp(-i.*s.*T./2);
end

%length(find(isnan(y)))
