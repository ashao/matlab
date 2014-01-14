function[bi,ci,ai]=ag2bc(ga,ar,be1,c1)
%AG2BC  Convert Morse wavelet parameters A and Gamma to Beta and C.   
%
%   AG2BC is used for determining parameter values of the Generalized
%   Morse Wavelets of Olhede and Walden (2002). 
%  
%   [BI,CI]=AG2BC(GA,AR,B,C) where GA is a "Gamma" parameter value and
%   "AR" is an "Area" parameter value, and B and C are arrays of "Beta"
%   and "C" parameter values, returns arrays BI and CI such that
%   MORSEAREA(CI,GA,BI) is approximately equal to AR.  BI and CI are of
%   size LENGTH(B) x 1.  
%
%   In other words, AG2BC solves for a curve C as a function of Beta
%   for fixed Gamma such that these three parameters create a set of
%   wavelets concentrated over a region having the specified area AR.
%  
%   [BI,CI,AI]=AG2BC(GA,AR,B,C) optionally returns the exact area AI
%   associated with BI, CI, and GA.
%
%   B and C are optional and have default values equal to 
%       C=LOGSPACE(0,3,1000)';    B=LINSPACE(1,20,1000)';
%
%   'ag2bc --t' runs a test.
%
%   Usage:  [bi,ci]=ag2bc(ga,ar,b,c);
%           [bi,ci]=ag2bc(ga,ar);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        

if strcmp(ga,'--t')
  ag2bc_test;return
end
  
if nargin==2
  c1=logspace(0,3,1000)';
  be1=linspace(1,20,1000)';
end
  
c=osum(c1,0*be1);
be=osum(0*c1,be1);
ga=1;
A = morsearea(c,ga(1),be);

[ci,bi,ai]=vzeros(size(A,1),length(ar));
for i=1:length(ar)
  [m,ii]=min(abs(A-ar(i)));
  index=sub2ind(size(A),ii,1:size(ii,2));
  ci(:,i)=c(index)';  
  bi(:,i)=be(index)';
  ai(:,i) = morsearea(ci(:,i),ga(1),bi(:,i));
end

function[]=ag2bc_test
%ar=[10 150];
ar=150;
[bi,ci,ai]=ag2bc(1,ar);

reporttest('AG2BC exact vs. approximate area', maxmax(abs(ai-ar)<2))

if 0
figure,
plot(bi,ai),hlines(ar)
title('Exact vs. approximate area')
xlabel('Beta')
ylabel('Area')
end
