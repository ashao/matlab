function[y]=fk(T,s)
% FK  The Fejer kernel FK(T,OM)
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        

y=frac(1,2*pi).*frac(sin2(s.*T./2),squared(s/2).*T); 

index=find(s==0);
if ~isempty(index)
  if length(T)>1
     y(index)=frac(T(index),2*pi);
  else
     y(index)=frac(T,2*pi);
  end
end


function[y]=sin2(x)
y=sin(x).^2;


