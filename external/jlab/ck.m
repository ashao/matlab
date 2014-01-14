function[y]=ck(T,s,so)
%CK The "C" kernel (difference of two Dirichlet kernels)
%
%   CK(T,OM,OM_R) is the "C" kernel defined in 
%
%            Lilly and Lettvin (2004), Signal Processing
%
%   If T is a scalar, the output is of size SIZE(OM).  If OM
%   is a scalar, the output is of size SIZE(T).  Otherwise the output
%   is a matrix of size LENGTH(T(:)) x LENGTH(S(:)).  
%  
%   OM_R must always be of size SIZE(OM) or a scalar.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002, 2004 J.M. Lilly --- type 'help jlab_license' for details        



if length(T)>1 && length(s)>1
      s=s(:);
      T=T(:);
      T1=T;
      s1=s;
      T=osum(T1,0*s1);
      s=osum(0*T1,s1);
end

y=sqrt(-1).*frac(pi,so).*(dk(T,s+so,1)-dk(T,s-so,1));
y(so==0)=0;
%length(find(isnan(y)))


