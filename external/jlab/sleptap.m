function[v,lambda]=sleptap(varargin)
%SLEPTAP  Calculate Slepian tapers.
%
%   [PSI,LAMBDA]=SLEPTAP(N,P,K) calculates the K lowest-order Slepian
%   tapers PSI of length N and time-bandwidth product P, together with
%   their eigenvalues LAMBDA. PSI is N x K and LAMBDA is K x 1.
%
%   K is optional and defaults to 2P-1.  
%   P is optional and defaults to 4.
%   
%   For N<256, SLEPTAP uses the tridiagonal method described in Percival 
%   and Walden (1993).  For N>256, it first computes tapers for N=256 and 
%   then spline-interpolates.
%   
%   Note that N may also be an array of lengths.  In this case PSI is a 
%   cell array of matrices, with PSI{1} being N(1) x K, PSI{2} being 
%   N(2) x K, etc.  Then LAMBDA is K x LENGTH(N).
%   _____________________________________________________________________
%
%   Normalization
%
%   By default, the tapers are set to have unit energy. Alternatively
%   SLEPTAP(...,'bandpass') uses the "bandpass" normalization in which the
%   tapers are rescaled so that the maximum value of the Fourier transform
%   of the first taper is set to two. 
%
%   See WAVETRANS for details on bandpass normalization.  
%   _____________________________________________________________________
%
%   See also MSPEC and MSVD.
%   
%   'sleptap --t' runs some tests.  
%
%   Usage:  [psi,lambda]=sleptap(n); 
%           [psi,lambda]=sleptap(n,p,k); 
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2011 J.M. Lilly --- type 'help jlab_license' for details        


if strcmp(varargin{1}, '--t')
  sleptap_test,return
end
n=varargin{1};
str='energy';
if isstr(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end
na=length(varargin);

%Nmax=256;
Nmax=512;
if na==1
	p=4;
else
    p=varargin{2};
end

if na<=2
	k=2*p-1;
else
    k=varargin{3};
end


%if n>Nmax, interpolate
if anyany(n>=Nmax)
    [v1,d]=sleptap_one(Nmax,p,k);
end

for j=1:length(n)
    if n(j)==Nmax
        v{j}=v1;
    elseif n(j)<=Nmax
        disp(['SLEPTAP calculating tapers of length ',int2str(n(j)),'.'])
        v{j}=sleptap_one(n(j),p,k);
    elseif n(j)>Nmax
        disp(['SLEPTAP interpolating to length ',int2str(n(j)),'.'])
        for i=1:k
            v{j}(:,i)=sinterp(v1(:,i),n(j));
        end
    end
end


%normalize
for j=1:length(n)
    for i=1:k
        v{j}(:,i)=v{j}(:,i)/sqrt(v{j}(:,i)'*v{j}(:,i));
        if v{j}(round(end/2),i)<0
           v{j}(:,i)=-v{j}(:,i);
        end
    end
end

if nargout==2
    lambda=zeros(k,length(n));
    for j=1:length(n)
        if n(j)>=Nmax
            lambda(:,j)=sleptap_lambda_one(Nmax,p,k,v1);
        else
            lambda(:,j)=sleptap_lambda_one(n(j),p,k,v{j});
        end
    end
end

if length(v)==1
    v=v{1};
end

%Bandpass normalization
if strfind(str,'ban')
    Vmax=max(abs(fft(v(:,1))));
    v=v*frac(2,Vmax);
end

function[v,d]=sleptap_one(n,p,k)
w=p./n;

%taper calculation using tridiagonal matrix
tic 
mat=zeros(n,n);
index=(1:n+1:n*n);
mat(index)=((n-1-2*(0:n-1))./2).^2.*cos(2*pi*w);
index2=index(1:length(index)-1)+1;
index3=index(2:length(index))-1;
mat(index2)=(1:n-1).*(n-(1:n-1))./2;
mat(index3)=(1:n-1).*(n-(1:n-1))./2;
toc

%tic
%t=[0:n-1]';
%tmat=osum(t,-t);
%mat=frac(sin(2*pi*w*tmat),pi*tmat);
%vswap(mat,nan,2*w);
%toc

OPTIONS.disp=0;
OPTIONS.maxit=2000;
%OPTIONS.tol=1e-8;

[v,d]=eigs(mat,k,'LA',OPTIONS);
%[v,d]=eigs(double(mat),max(2*p-1,k),'LM',OPTIONS);
%v=v(:,1:k);d=d(1:k);


function[lambda]=sleptap_lambda_one(n,p,k,v)
w=p./n;
i=(0:n-1)'*ones(1,n);
j=i';
A=pi*(i-j);
index=find(A==0);
A(index)=1;
A=sin(2*pi*w*(i-j))./A;
A(index)=2*w;
%find eigenvalues
for i=1:k
    %note unexplained matlab quirk: dividing
    %two column vectors gives you a column vector,
    %dividing two row vectors gives you a scalar
    lambda(i,1)=(A*v(:,i))'/v(:,i)';
end

%/***************************************************
%here's some garbage
if 0
%see how much spline-interpolated ones vary from others
xx=v(:,1);xx=xx/sqrt(xx'*xx);figure,plot(xx)
xx=diff(xx);xx=xx/sqrt(xx'*xx);hold on,plot(xx,'g')
xx=diff(xx);xx=xx/sqrt(xx'*xx);plot(xx,'r')
xx=diff(xx);xx=xx/sqrt(xx'*xx);plot(xx,'c');

for i=1:4
v(:,i)=v(:,i)/sqrt(v(:,i)'*v(:,i));
end

figure,plot(v)


l1=lambda;
v1=v;


k=4;
n=256;
w=4/n;
mat=zeros(n,n);
index=(1:n+1:n*n);
mat(index)=((n-1-2*(0:n-1))./2).^2.*cos(2*pi*w);
index2=index(1:length(index)-1)+1;
index3=index(2:length(index))-1;
mat(index2)=(1:n-1).*(n-(1:n-1))./2;
mat(index3)=(1:n-1).*(n-(1:n-1))./2;
[v,d]=eigs(mat,k);
%[d,index]=sort(diag(d));
%v=v(:,index);
%d=flipud(d);
%v=fliplr(v);
for i=1:4
v(:,i)=v(:,i)/sqrt(v(:,i)'*v(:,i));
end

A=calcdefmat(n,w);

for i=1:k
	lambda(i,1)=mean((A*v(:,i))\v(:,i));
end
v=v(:,1:k);
 
for i=1:4
	v1i(:,i)=interp1((1:100)/100',v1(:,i),(1:256)'/256,'cubic');
end
for i=1:4
v1i(:,i)=v1i(:,i)/sqrt(v1i(:,i)'*v1i(:,i));
end
end
%end garbage
%\***************************************************



function[a]=calcdefmat(n,w)
i=(0:n-1)'*ones(1,n);
j=i';
a=sin(2*pi*w*(i-j))./(pi*(i-j));
a(isnan(a))=2*w;


function[]=sleptap_test

[psi,lambda]=sleptap(200);

tol=1e-6;
bool=false(1,size(psi,2));
for j=1:size(psi,2)       
        bool(1,j)=aresame(vsum(psi(:,j).^2,1),1,tol);
end
reporttest('SLEPTAP unit energy',allall(bool))

[psi,lambda]=sleptap([200 512 1024]);

tol=1e-6;
bool=false(length(psi),size(psi{1},2));
for i=1:length(psi)
    for j=1:size(psi{1},2)       
        bool(i,j)=aresame(vsum(psi{i}(:,j).^2,1),1,tol);
    end
end

reporttest('SLEPTAP unit energy with interpolation & cell output',allall(bool))

[psi,lambda]=sleptap(200,4,1,'bandpass');

reporttest('SLEPTAP bandpass normalization',aresame(maxmax(abs(fft(psi))),2,1e-10))
