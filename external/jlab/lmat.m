function[x]=lmat(theta)
%LMAT 2x2 phase twist matrix through specified angle.
%
%   L=LMAT(PHI) creates a phase shift matrix
%
%       L=[  COS(PHI)    i SIN(PHI);
%          i SIN(PHI)      COS(PHI)]
%
%   such that L*X phase-shifts the first element of the column vector X
%   forward by PHI radians, and phase-shfits the second the second
%   element backwards by PHI radians.
%
%   If LENGTH(PHI)>1, then L will have dimension 2 x 2 x LENGTH(PHI).
%
%   See also JPOLY, VECTMULT, JMAT3, IMAT, JMAT, KMAT, and TMAT.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details        

if length(theta)>1;
  theta=theta(:);
end

x=zeros(2,2,length(theta));
x(1,1,:)=cos(theta);
x(1,2,:)=sqrt(-1)*sin(theta);
x(2,1,:)=sqrt(-1)*sin(theta);
x(2,2,:)=cos(theta);


if size(theta,1)~=length(theta(:));
    x=reshape(x,[2,2,size(theta)]);
end