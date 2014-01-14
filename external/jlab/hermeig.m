function l = hermeig(R,j,str)
%HERMEIG  Eigenvalues of orthonormal Hermite functions.
%
%   L=HERMEIG(R,N) returns the N+1 lowest-order eigenvalues of the
%   orthogonal Hermite windowing functions for and area R.  
%  
%   L=HERMEIG(R,N,STR) returns specificies the algorithm to use. 
%   STR='Bayram', the default, uses that of Bayram and Baraniuk 2001.
%   STR='Simons', uses the (incorrect) version of Simons et al. 2003.
%
%   R may be an array but N must be a scalar.
%
%   L has LENGTH(R(:)) rows and N+1 colums.  Note that L(1) is the 
%   eigenvalue of the 'zeroth' order Hermite function, etc. 
%  
%   See also HERMFUN, HERMPOLY.
%
%   'hermeig --f' generates a test figure, Fig. 2c of Simons et
%   al. 2003, using both algorithms.
%
%   Usage:  lambda=hermeig(r,n);
%           lambda=hermeig(r,n,'bay');
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004--2007 F. Rekibi and J. M. Lilly 
%                          --- type 'help jlab_license' for details
   
% 05.08.07  JML fixed bug to include N+1 columns

if strcmp(R,'--f')
  hermeig_fig;return
end

if size(R,1)>1
	R=R(:);
end

l=zeros(j+1,length(R));

if nargin==2
	str='Bayram';
end

if strcmp(str,'Simons')
	for i=0:j
		l(i+1,:)=gammainc(R.^2/2,i+1)./factorial(i);
	end
elseif strcmp(str,'Bayram')
	S=0;
	for i=0:j
		S=S+2^(-i).*R.^(2*i)./factorial(i);
		l(i+1,:)=1-exp(-R.^2/2).*S; %lambda(i)
	end
end

l=l';

function[]=hermeig_fig

J=15;
R=[2 3 4];
L1= hermeig(R,J);
L2= hermeig(R,J,'Simons');

%start numbering at zero
t=(0:J);
figure,
subplot(121),plot(t,L1)
title('Bayram algorithm')
subplot(122),plot(t,L2)
title('Simons algorithm')
