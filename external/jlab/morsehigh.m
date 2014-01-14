function[f]=morsehigh(varargin)
%MORSEHIGH  High-frequency cutoff of the generalized Morse wavelets.
%
%   FALPHA=MORSEHIGH(GAMMA,BETA,ALPHA) returns the high frequency cutoff
%   FALPHA of the generalized Morse wavelet specified by GAMMA and BETA, 
%   with cutoff level FALPHA.
%
%   Specifically, if PSI is the wavelet and PSIMAX is its maximum value, 
%   then FALPHA is the highest *radian* frequency at which 
%
%      PSI(FALPHA)/PSIMAX > ALPHA.
%
%   This gives a way to choose the high-frequency cutoff in the wavelet
%   transform.  See Lilly and Olhede (2009d) for details. 
%  
%   The input parameters may either all be scalars, or GAMMA and BETA
%   may be scalars of the same size with scalar ALPHA.
%   ___________________________________________________________________
%
%   Precision vs. speed
%
%   MORSEHIGH(..., N) uses 1/N times the peak frequency MORSEFREQ as the
%   numerical interval.  N=100 is the default; choose a smaller value
%   for faster speed but diminished precision. 
%   ___________________________________________________________________
%  
%   See also MORSEFREQ, MORSEWAVE.
%
%   'morsehigh --t' runs a test.
%   'morsehigh --f' generates some sample figures.
%
%   Usage: falpha=morsehigh(ga,be,alpha);
%          falpha=morsehigh(ga,be,alpha,N);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2009 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--f')
    morsehigh_figure,return
elseif strcmp(varargin{1}, '--t')
    morsehigh_test,return
end

gamma=varargin{1};
beta=varargin{2};
alpha=varargin{3};
if nargin>3
    N=varargin{4};
end
    
if nargin~=3&&nargin~=4;
    error('MORSEHIGH takes either three or four input arguments.')
end

ompeak=morsefreq(gamma,beta);
N=100;

dom=vrep(ompeak/N,10*N,3);
dom(:,:,1)=ompeak;

ommat=cumsum(dom,3);

amat=vrep(morsea(gamma,beta),10*N,3);
betamat=vrep(beta,10*N,3);
gammamat=vrep(gamma,10*N,3);

morse=frac(1,2)*amat.*(ommat.^betamat).*exp(-ommat.^gammamat);

kk=vsum(morse>alpha,3);
ii=vrep((1:size(gamma,1))',size(gamma,2),2);
jj=vrep(1:size(gamma,2),size(gamma,1),1);

index=sub2ind(size(morse),ii,jj,kk);

%f=frac(ommat(index),ompeak);
f=ommat(index);

function[]=morsehigh_figure
morsehigh_figure1;
morsehigh_figure2;

function[]=morsehigh_figure1

ga1=(1:.1:11);
be1=(1:.1:10);

[ga,be]=meshgrid(ga1,be1);

ompeak=morsefreq(ga,be);
ommax=morsehigh(ga,be,0.1,10);

figure
%p=morseprops(ga,be);
contourf(ga1,be1,(ompeak./ommax),20),colorbar,nocontours

title('Morse Wavelet  f_{peak} / f_{0.1} )')
xlabel('Gamma Parameter')
ylabel('Beta Parameter')

disp('This is the ratio of the peak frequency to the ALPHA=0.1 frequency,')
disp('that is, the frequency at which the wavelet has decayed to 0.1 of its')
disp('peak value.  Large values of this ratio mean long high-frequency tails.')

function[]=morsehigh_test

ga=3;
be=4;
alpha=0.1;
ompeak=morsefreq(ga,be);
N=100;

dom=ompeak/N+zeros(4*N,1);
om=cumsum(dom,1);

morse=frac(1,2)*morsea(ga,be).*(om.^be).*exp(-om.^ga);

fhigh=morsehigh(ga,be,0.1);
reporttest('MORSEHIGH',aresame(fhigh,om(max(find(morse>alpha))),1e-4))


function[]=morsehigh_figure2

ga=3;
be=4;
alpha=0.1;
ompeak=morsefreq(ga,be);
N=100;

dom=ompeak/N+zeros(4*N,1);
om=cumsum(dom,1);

morse=frac(1,2)*morsea(ga,be).*(om.^be).*exp(-om.^ga);
fhigh=morsehigh(ga,be,alpha);
figure,
plot(om,morse),xlabel('Radian frequency')
title('MORSEHIGH frequency with \alpha=0.05')
vlines(fhigh),hlines(.1)


