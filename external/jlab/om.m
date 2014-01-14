function[s,ds]=om(k)
%OM  Dispersion relation for gravity-capillary waves.
%
%   S=OM(K) where K is an array of wavenumbers, returns the 
%   frequency S  for a gravity-capillary wave.  The units of K
%   are [rad cm^-1] and the units of S are  [rad s^-1].
%
%   S=OM with no arguments returns the frequency at which the
%   phase speed is minimized, equal to (4*G^3./T)^(1/4).
%
%   [S,DS]=OM(K) also returns DS/DK, the derivative of the 
%   frequency with respect to wavenumber.
%
%   The wavenumber array may be complex-valued, K = Kx + i*Ky.
%  
%   See also GC_PARAMS, KMIN.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    


%/*******************************************************
[g,T,km,sm]=gc_params;  %cgs units

if nargin==0
    s=(4*g.^3./T).^(1/4);
else
    k=abs(k);
    s=sqrt(g.*k+T.*k.*k.*k);    
    if nargout>1
      ds=0.5.*(1./s).*(g+3*T.*k.^2);
    end
end

%\*******************************************************
