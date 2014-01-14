function[x]=kmat(theta)
%KMAT 2x2 phase shift matrix through specified angle.
%
%   K=KMAT(PHI) creates a phase shift matrix
%
%       K=[e.^(i*PHI)     0;
%          0     e.^(-i*PHI)]
%
%   such that K*X phase-shifts the first element of the column vector X
%   forward by PHI radians, and phase-shfits the second the second
%   element backwards by PHI radians.
%
%   If LENGTH(PHI)>1, then K will have dimension 2 x 2 x LENGTH(PHI).
%
%   See also JPOLY, VECTMULT, JMAT, IMAT, and TMAT.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2010 J.M. Lilly --- type 'help jlab_license' for details        

if length(theta)>1;
  theta=theta(:);
end

x=zeros(2,2,length(theta));
x(1,1,:)=rot(theta);
x(2,2,:)=rot(-theta);


if size(theta,1)~=length(theta(:));
    x=reshape(x,[2,2,size(theta)]);
end