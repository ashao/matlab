function h=hermfun(t,j)
%HERMFUN  Orthonormal Hermite functions. [with F. Rekibi]
%
%   H=HERMFUN(T,N) generates the orthonormal Hermite functions
%   [H0,...HN] on a time axis specfied by the column vector T.
%
%   HERMFUN uses the expression of Simons et al. 2003.
%  
%   Note that H(:,1) is the 'zeroth' order Hermite function, etc.
%   H thus has N+1 columns.
%  
%   See also HERMEIG, HERMPOLY.
%
%   'hermfun --f' generates a sample figure; compare with the
%   Hermite function figure at
%
%     http://en.wikipedia.org/wiki/Hermite_polynomials#Definition
%
%   Usage:  h=hermfun(t,n);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2004--2007 F. Rekibi and J. M. Lilly
%                         --- type 'help jlab_license' for details

% 05.08.07  JML fixed bug to include N+1 columns

%   'hermfun --f' generates a sample figure; compare with Figure 2 
%   of Simons, van der Hilst, and Zuber (2003), JGR 108 B5.

if strcmp(t,'--f')
  hermfun_fig;return
end

if size(t,1)==1
    t=t';
end

H=hermpoly(t,j);

E=exp(-t.^2/2)*ones(1,j+1);
HE=H.*E;

h=zeros(length(t),j);
for k=1:j+1
	h(:,k)=HE(:,k)*frac(1,pi^(1/4))*frac(1,sqrt(2^(k-1)*factorial(k-1)));
end


% Normalisation

function[]=hermfun_fig
t=(-5:0.1:5.1);
tnorm=(t-t(1))/(t(length(t))-t(1));


h=hermfun(t,5);
figure,plot(t,h),linestyle k r b y g m
title('Hermite functions H_0 -- H_5')

function[]=hermfun_fig2

t=(-5:0.1:5.1);
R=(2:4);

tnorm=(t-t(1))/(t(length(t))-t(1));

j=5;
h=hermfun(t,j);
figure,plot(t,h)
for i=1:size(h,2)
	h(:,i)=h(:,i)/max(abs(h(:,i)));
end


% Trace des 5 premieres fonctions de Hermite ortho
figure
subplot 211
plot(tnorm,h);  
linestyle r b g m-- y--
legend('j=0','j=1','j=2','j=3','j=4',-1);
xlabel('Normalized window length')
ylabel('Normalized value')

% Trace des valeurs propres de la fonction de Hermite

j=15;
l = hermeig(R,j);

subplot 223
T=(0:j);
plot(T,l)
legend('R=2','R=3','R=4',0);
xlabel('Number');
ylabel('lambda R(j)')

% Trace de l'energie

subplot 224

U=zeros(length(t),length(R));
for i=1:length(R)
  htemp=hermfun(t,R(i).^2).^2;
  eigtemp=hermeig(R(i),R(i).^2);
  U(:,i)=frac(1,R(i).^2).*(htemp*eigtemp');
end
plot(tnorm,U);  
linestyle r b g m-- y--
xlabel('Normalized window length')
ylabel('Energy');

% subplot 325
% j=6;
% plot(h(:,size(h,2)));
% axis tight
% title('6 eme fct de Hermite');
% 
% subplot 326
% j=7;
% h=hermfun(t,j);
% W = WignerDist(h(:,end));

