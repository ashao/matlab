function[K]=sub2k(m,n,dk,N)
%SUB2K Convert an I,J subscript pair into a wavenumber
%  
%   SUB2K(I,J,DK,N), where I and J are indices into an NxN matrix, 
%   returns the complex-valued wavenumber
%       K=(J+iI)*DK - (1+i)*(N+1)/2*DK
%   which corresponds to this location on a regular wavenumber grid.
%   Note the row index "I" determines a Y-wavenumber.
%
%   SUB2K is inverted by K2SUB.  See also WAVEGRID.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details      
  
if strcmp(m, '--t')
  sub2k_test,return
end

K=(n*dk-(N+1)/2*dk)+sqrt(-1)*(m*dk-(N+1)/2*dk);


function[]=sub2k_test

x=sub2k((1:5),0*(1:5)+3,1,5);
bool(1)=aresame(x,sqrt(-1)*[-2 -1 0 1 2]);
%disp('Should be sqrt(-1) x -2 -1 0 1 2')

x=sub2k((1:4),0*(1:4)+2,1,4);
bool(2)=aresame(x,-0.5 + sqrt(-1)*[-1.5 -0.5 0.5 1.5]);
%disp('Should be -0.5 + sqrt(-1) x -1.5 -0.5 0.5 1.5')

x=sub2k((1:4),0*(1:4)+2,pi,4)./pi;
bool(3)=aresame(x,-.5 + sqrt(-1)*[-1.5 -0.5 0.5 1.5]);
%disp('Should be -0.5 + sqrt(-1) x -1.5 -0.5 0.5 1.5')

x=sub2k(0*(1:5)+2+1,(1:5),1,5);
bool(4)=aresame(x,[-2 -1 0 1 2]);
%disp('Should be -2 -1 0 1 2')

x=sub2k(0*(1:4)+2,(1:4),1,4);
bool(5)=aresame(x,[-1.5 -.5 .5 1.5] -sqrt(-1)*0.5);
%disp('Should be -1.5 -0.5 0.5 1.5 - sqrt(-1) x 0.5')
reporttest('SUB2K',all(bool))
  
  
