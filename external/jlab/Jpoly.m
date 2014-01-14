%
%  JPOLY  Special polynomials, matrices, and functions.
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, *JPOLY*, JGRAPH, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
% 
%  Polynomials
%   laguerre   - Generalized Laguerre polynomials.  [co-authored with F. Rekibi]
%   hermpoly   - Hermite polynomials.  [co-authored with F. Rekibi]
%   bellpoly   - Complete Bell polynomials.
%
%  Special functions
%   dawson      - The Dawson function. [By P.J. Acklam]
%   dawsonderiv - Derivatives of the Dawson function.
%
%  Special 2 x 2 matrices
%   imat       - 2x2 identify matrix.
%   jmat       - 2x2 rotation matrix through specified angle.
%   kmat       - 2x2 phase shift matrix through specified angle.
%   lmat       - 2x2 Hilbert reflections matrix through specified angle.
%   tmat       - 2x2 complex grouping matrix.
%
%  Special 3 x 3 matrices
%   jmat3      - 3x3 rotation matrix through specified angle.
%
%  Special N x N matrices
%   amat       - Mean-taking matrix.
%   fmat       - Unitary Fourier transform matrix.
%   hmat       - Hilbert transform matrix.
%
%  Signal processing kernel functions [see Lilly and Lettvin (2004)]
%   ck         - The "C" kernel (difference of two Dirichlet kernels).
%   dk         - The Dirichlet kernel.
%   fk         - The Fejer kernel.
%   hk         - The "H" kernel (product of two Dirichlet kernels).
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details   

help Jpoly

if 0  
             %Polynomials
             laguerre   %- Generalized Laguerre polynomials.  [co%-authored with F. Rekibi]
             hermpoly   %- Hermite polynomials.  [co%-authored with F. Rekibi]
             bellpoly   %- Complete Bell polynomials.
          
            %Special functions
             dawson      %- The Dawson function. [By P.J. Acklam]
             dawsonderiv %- Derivatives of the Dawson function.
          
            %Special 2 x 2 matrices
             imat       %- 2x2 identify matrix.
             jmat       %- 2x2 rotation matrix through specified angle.
             kmat       %- 2x2 phase shift matrix through specified angle.
             lmat       %- 2x2 Hilbert reflection matrix through specified angle.
             tmat       %- 2x2 complex grouping matrix.
          
            %Special 3 x 3 matrices
             jmat3      %- 3x3 rotation matrix through specified angle.
          
            %Special N x N matrices
             amat       %- Mean%-taking matrix.
             fmat       %- Unitary Fourier transform matrix.
             hmat       %- Hilbert transform matrix.
          
            %Signal processing kernel functions [see Lilly and Lettvin (2004)]
             ck         %- The "C" kernel (difference of two Dirichlet kernels).
             dk         %- The Dirichlet kernel.
             fk         %- The Fejer kernel.
             hk         %- The "H" kernel (product of two Dirichlet kernels).
  end