

function[y]=dawsonderiv(varargin)
%DAWSONDERIV Derivatives of the Dawson function.
%
%   DCELL=DAWSONDERIV(X,N) returns the first N derivatives of the 
%   Dawson function at elements X. 
%
%   All elements of DCELL have the same size as X.  DCELL{1} is the 
%   first derivative, and so forth.
%
%   A form for the derivatives of the DAWSON function is given by 
%   Lilly and Olhede (2008b), and involves Hermite polynomials.
%
%   See also DAWSON, HERMPOLY.
%
%   'dawsonderiv --t' runs a test.
%   'dawsonderiv --f' generates a sample figure.
%
%   Usage: dn=dawsonderiv(x,n);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    dawsonderiv_test,return
end
if strcmp(varargin{1}, '--f')
    dawsonderiv_fig,return
end
 
x=varargin{1};
n=varargin{2};

herm=hermpoly(x(:),n);
iherm=hermpoly(sqrt(-1)*x(:),n);

for k=1:n+1
    hermcell{k}=reshape(herm(:,k),size(x));
    ihermcell{k}=reshape(iherm(:,k),size(x));
    %real(sqrt(-1).^(k-1).*reshape(iherm(:,k),size(x)));
end
y=zeros(size(x));
for k=1:n
     y=y+choose(n,k).*hermcell{n-k+1}.*((sqrt(-1).^(k-1)).*ihermcell{k});
     %                      H_(n-k}                               H{k-1}
end
y=((-1).^n).*(hermcell{n+1}.*dawson(x)-y);

function[]=dawsonderiv_test
dt=0.01;
t=(-15:dt:15)';
dk1=dawson(t);
dk1([1 end],:)=nan;  %Avoid differentiation errors at edges
for k=1:5
    dk2=dawsonderiv(t,k);
    dk1=vdiff(dk1,1)./dt;
    err=vsum(abs(dk1-dk2).^2,1)./vsum(abs(dk1).^2,1);
    reporttest(['DAWSONDERIV matches numerical differentiation for n=' int2str(k)],err<1e-3)
end

function[]=dawsonderiv_fig
dt=0.01;
t=(-10:dt:10)';
dk(:,1)=dawson(t);
dk(:,1)=dk(:,1)./maxmax(dk(:,1));
for k=1:5
    dk(:,k+1)=dawsonderiv(t,k);
    dk(:,k+1)=dk(:,k+1)./maxmax(dk(:,k+1));
end
figure,
plot(t,dk),title('The Dawson function and first 4 (normalized) derivatives'),ylim([-2 1.1])

%plot(t,[psi1 psi2])
