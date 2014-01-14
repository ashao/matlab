function[a]=triadevolve(ahat,k,t)
%TRIADEVOLVE McGoldrick's evolution functions for a resonant wave triad.
%
%   A=TRIADEVOLVE(AHAT,K,T) implements McGoldrick's (1965) 
%   equation 4.6 for the evolution of a resonant wave triad, 
%   neglecting viscosity and assuming one member of the triad 
%   initially has zero amplitude.
%
%   Here AHAT and K are arrays of length 2 specifying the 
%   initial amplitudes and wavenumbers, respectively, of two 
%   waves.  The units of K are radians cm^-1.
%
%   Note that my definition of AHAT differs from McGoldrick's
%   by a factor of two; see Lilly and Lettvin (2006).
%
%   See also RESCOEFF.
%
%   Usage:  a=triadevolve(ahat,k,t);
%    
%   'triadevolve --f' makes a sample figure.  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001--2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(ahat, '--t'),return,end
if strcmp(ahat,'--f')
  triadevolve_fig;return
end
  
k(3)=k(1)+k(2);

[g3,g2,g1]=rescoeff(k(1),k(2),k(3),'mcg');
B12=g3/2;
B13=g2/2;
B23=g1/2;

m=B23.*ahat(2).^2./(B13.*ahat(1).^2);

if m<1
	Theta=2*real(ahat(1).*sqrt(B12.*B13).*t);
	[sn,cn,dn] = ellipj(Theta,m);
	a(:,1)=ahat(1).*dn;
	a(:,2)=ahat(2).*cn;
	a(:,3)=ahat(2).*sn*sqrt(B12./B13);
else
	Theta=real(ahat(2).*sqrt(B12.*B23).*t);
	[sn,cn,dn] = ellipj(Theta,1./m);
	a(:,1)=ahat(1).*cn;%note the switch here
	a(:,2)=ahat(2).*dn;
	a(:,3)=ahat(1).*sn*sqrt(B12./B23);
end



function[]=triadevolve_fig
%/********************************************************
k=[3.667*exp(sqrt(-1)*47*pi/360);3.667*exp(-sqrt(-1)*47*pi/360)];
k(3)=k(1)+k(2);
ahat=[.687 .487 0]'/10/2;%convert to cm
mu=[0 0 -pi/2]';  
dmu=mu(1)+mu(2)-mu(3); %dmu is +pi/2

K=[-flipud(k);k];

[g3,g2,g1]=rescoeff(k(1),k(2),k(3),'mcg');
B12=g3/2;B13=g2/2;B23=g1/2;

p=2*pi/om(k(3));

t=p*(0:.1:32)';
a=triadevolve(ahat(1:2),k(1:2),t);
X=conj([a(:,1).*rot(mu(1)) a(:,2).*rot(mu(2)) a(:,3).*rot(mu(3))])';
X=[conj(flipud(X));X];

[S,B,B1]=dmspec(1,K,X);

clear sig asym
sig=sqrt(vsum(abs(X.^2),1));
asym=dmasym(1,K,B1,sig,rot(0));

if 0
tf=p*(0:.5:10)';
psi=2*[0*ahat ;ahat.^2];
[Bf,B1f]=psi2b(tf,K,psi,'direct'); 
asymf=dmasym(1,K,B1f,sig(1),rot(0));
end
%\********************************************************
%figure,plot(t,imag(B1)),hold on,plot(tf,imag(B1f))
figure,
subplot(211)
plot(t./p,a./abs(ahat(1)));
linestyle three 
axis([0 32 -1.1 1.1])
title('Resonant wave triad')
grid off
set(gca,'xtick',(0:4:32),'box','on')
hlines(0,'k:')
fixlabels([0 -1])
ylabel('Wave amplitude' )
set(gca,'xticklabel',[])
vlines(3.6,'k:')

subplot(212)
plot(t./p,asym,'k')
set(gca,'ytick',(-.8:.2:.8))
axis([0 32 -0.75 0.75])
%axis([0 32 -2 2])
grid off
set(gca,'xtick',(0:4:32),'box','on')
hlines(0,'k:')
vlines(3.6,'k:')
h=packrows(2,1);

fixlabels([0 -1])

letterlabels(3)
ylabel('Asymmetry \gamma_a(i)' )
xlabel('Nondimensional time (periods of wave three)')


ka=2*kfun(1,1,k(1),k(2));
m=a(1,1)*a(1,2)*ka;
axes(h(1)),plot(t./p,m.*t./abs(ahat(1)),'k-.')
m=3*(2*ahat(1).^2)*(2*ahat(2).^2)*ka./(sig(1,1).^3);
axes(h(2)),plot(t./p,m.*t,'k-.')

t1=min(find(a(:,2)<0));
tnew=t;tnew(t<t(t1))=nan;
ka=2*kfun(-1,1,k(1),k(3));  
m=a(t1,1)*a(t1,3)*ka;
axes(h(1)),plot(tnew./p,-m*(tnew-t(t1))./abs(ahat(1)),'k-.')
m=-3*(2*a(t1,1).^2)*(2*a(t1,3).^2)*ka./(sig(t1).^3);
axes(h(2)),plot(tnew./p,m.*(tnew-t(t1)),'k-.')


fontsize jpofigure
set(gcf,'paperposition',[2 2 3.5 5])
%print -dpsc mcgfig5.ps 
%!gv mcgfig5.ps &

