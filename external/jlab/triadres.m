function[th,k2]=triadres(k1,k2,s)
%TRIADRES  Solves the triad resonance condition given two waves.
%
%   TH=TRIADRES(K1,K2) where K1 and K2 are two wavenumbers with
%   units [rad cm^-1] returns the angle TH between wave one and 
%   wave two for which the resonance condition, McGoldrick's 
%   (1965) equation 2.3, is satisfied.
%
%   K1 and K2 may each either be a scalar or an array.  TH is
%   a matrix of size LENGTH(K2) by LENGTH(K1).
%
%   [TH,CK2]=TRIADRES(K1,K2,S) solves for the sum interaction if 
%   S==1 and the difference interaction if S==-1, and returns
%   the complex-valued matrix CK2 of SIZE(TH) containing the 
%   complex-valued K2 vectors satisfying the resonace condition
%   when K1 is purely real and positive.  
%
%   See also VTRIADRES, TRIADEVOLVE.
%
%   Usage: th=triadres(k1,k2);
%
%   'triadres --t' runs a test.
%   'triadres --f' makes a sample figure.  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(k1,'--t')
    triadres_test;return
end
if strcmp(k1,'--f')
    triadres_fig;return
end


%Hack to prevent failure at zero 
k1(k1==0)=1e-10;


k1=row2col(k1);
k2=col2row(k2);
[k2,k1]=ndgrid(k2,k1);
[N,M]=size(k1);

a=k1./kmin;
b=k2./kmin;

if nargin==2
    s=1;
end

a2=a.^2;
a3=a.^3;
a4=a.^4;
b2=b.^2;
b3=b.^3;
b4=b.^4;

D=s*8.*a2.*b2;
C=4*(3*a3.*b+2*a.*b+3*a.*b3);
B=s*2*(1+4*a2+4*b2+3*a4+3*b4+6*a2.*b2);
A=-6+4*a.*b-6*a2.*b2+3*a3.*b+3*a.*b3-6*a2-6*b2 ...
    -1*s*4*sqrt((a+a3).*(b+b3)).*((a+b)./a./b+(a3+b3)./a./b);

th=nan*ones(N,M);
for i=1:N
    for j=1:M
	cthi=roots([D(i,j),C(i,j),B(i,j),A(i,j)]);
        %note can't use isreal since it returns a scalar
        index=find(imag(cthi)==0&abs(cthi)<1);
        if length(index)>=1
           th(i,j)=acos(cthi(index(1)));
        end
    end    
end

th(th==0)=nan;
k2=k2.*exp(sqrt(-1)*th);


function[]=triadres_fig

k1=[3 5/2 2 3/2 1 2/3 1/2 2/5 1/3].*kmin;
k2=(.0001:.05:4)'.*kmin;


disp(char(10))
disp('Recreates McGoldrick''s version of the resonant condition hodograph, his')
disp('Figure 1, using identical panels.  This is shown in the left panel.')  
disp('This is for sum-type interactions, defined by McGoldricks (2.2).  The')
disp('right hand panel is for difference-type interactions, obtained by changing')
disp('the sign of a few terms in McG (2.3).') 


figure
subplot(121)
[th,ck2]=triadres(k1,k2);
polar(0,4),hold on,plot([flipud(conj(ck2));ck2]/kmin),axis equal
title('McGolrick Figure 1; Sum interactions')
text(2,5,'k1=[3 5/2 2 3/2 1 2/3 1/2 2/5 1/3]')

subplot(122)
[th,ck2]=triadres(k1,k2,-1);
polar(0,4),hold on,plot([flipud(conj(ck2));ck2]./kmin),axis equal
title('Difference Interactions')

%[th,ck2]=triadres(k1,k2);
%k3=k1+k2;

function[]=triadres_test

k1=[3 5/2 2 3/2 1 2/3 1/2 2/5 1/3].*kmin;
k2=(.0001:.05:4)'.*kmin;

ck1=ones(length(k2),1)*k1;

tol=1e-8;
[th,ck2]=triadres(k1,k2);
bool=aresame(om(ck1)+om(ck2),om(ck1+ck2),tol);
reporttest('TRIADRES sum triads satisfy resonance',bool)

[th,ck2]=triadres(k1,k2,-1);
bool=aresame(abs(om(ck2)-om(ck1)),om(ck1-ck2),tol);
reporttest('TRIADRES difference triads satisfy resonance',bool)
