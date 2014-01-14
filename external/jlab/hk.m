function[y]=hk(T,s,so)
%HK The "H" kernel (product of two Dirichlet kernels)
%
%    HK(T,OM,OM_R)=
%         4*PI*REAL(DK(T,OM+OM_R,1).*CONJ(DK(T,OM-OM_R,1)));
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        

  

if length(T)>1
  if length(s)>1
    if jiscol(T) && jisrow(s)
      T1=T;
      s1=s;
      T=T1*ones(size(s1));
      s=ones(size(T1))*s1;
    end
  end
end

y=4*pi.*real(dk(T,s+so,1).*conj(dk(T,s-so,1)));




