function A = morsearea(C,ga,be)
%MORSEAREA  Time-frequency concentration area of Morse wavelets. [with F. Rekibi]
%
%   A=MORSEAREA(C,GAMMA,BETA) calculates the area of time/frequency
%   concentration region for the generalized Morse wavelets specified
%   by parameters C, GAMMA, and BETA. 
%  
%   The input parameters must either be matrices of the same size, or
%   some may be matrices and the others scalars.  
% 
%   MORSEAREA uses the area formula of Olhede and Walden (2002),
%   "Generalized Morse Wavelets", at the bottom right of page 2664.
%
%   Usage: A = morsearea(C,ga,be);
%
%   'morsearea --f' generates a sample figure
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 F. Rekibi and J. M. Lilly 
%                         --- type 'help jlab_license' for details  

if strcmp(C,'--f')
  morsearea_fig;return
end

r=((2*be)+1)./ga;
A=2*pi*(C-1).*gamma(r+1-(1./ga)).*gamma(r+(1./ga))./(ga.*gamma(r).^2);

function[]= morsearea_fig
C1=logspace(0.1,3,100)';
be1=linspace(1,10,101)';
C=osum(C1,0*be1);
be=osum(0*C1,be1);
ga=2;


figure
A = morsearea(C,ga(1),be);
contourf(C1,be1,log10(A'),20)
xlog
xlabel('Parameter C')
ylabel('Parameter \beta')
title('Log10 Area with \gamma =2')
colorbar
