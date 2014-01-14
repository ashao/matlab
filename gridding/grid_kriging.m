function [ Zz Zzstd] = grid_kriging( x,y,z, Xx, Yy,beta )
%KRIGINGGRID Grids data in 2D by Kriging Interpolation
%   Input:
%       x, y, z: Input fields (vectors)
%       Xx, Yy: Interpolaton points (vectors)
%   Output:
%       Zz: Interpolated values.

% Calculate distances between data points and interpolation points
%

%   Create vector of data points
nanidx=find(isnan(z));
x(nanidx)=[];
y(nanidx)=[];
z(nanidx)=[];
Z=[z ; 0];

%   Estimate variogram:
[X1 X2]=meshgrid(x);
[Y1 Y2]=meshgrid(y);
[Z1 Z2]=meshgrid(z);
D=sqrt((X1-X2).^2+(Y1-Y2).^2);
G=0.5*(Z1-Z2).^2;
idx=1:length(z);
D2=D.*(diag(x*NaN)+1);
lag=mean(min(D2));
hmd=max(D(:))/2;
max_lags=floor(hmd/lag);
LAGS=ceil(D/lag);
for i=1:max_lags;
    SEL=(LAGS==i);
    DE(i)=mean(mean(D(SEL)));
    GE(i)=mean(mean(G(SEL)));
end

% Fit to Exponential variogram model: alpha*r^beta
disp('Estimating alpha via least squares')
alpha=fitpow(DE,GE,beta);
disp(['Alpha=' num2str(alpha)])

%   Create vector of interpolation points
[Xx Yy]=meshgrid(Xx,Yy);
intpt=[Xx(:) Yy(:)];
size(intpt)
%   Kriging Interpolation

n=length(xypt);
disp('Initializing V matrix')
V=zeros(n+1,n+1);
for i=1:n
    for j=1:n
        V(i,j)=rcalc(xypt(i,:)-xypt(j,:));
    end
end

V(n+1,:)=1;
V(:,n+1)=1;
V(n+1,n+1)=0;

Vinv=inv(V);

numintpt=length(intpt);
for i=1:n
    if mod(i,100)==1
        %	disp(['Calculating Variogram Point ' num2str(i) '/' num2str(n)])
    end
    for j=1:numintpt
        Vstar(j,i)=powvar(alpha,rcalc(intpt(j,:)-xypt(i,:)),beta);
    end
end
Vstar(:,n+1)=1;
disp('Calculating interpolated Z values')
Zz=Vstar*Vinv*Z;
disp('Calculating Errors')
Zzstd=sqrt(Vstar*Vinv*Z);

end

function [ alpha ] = fitpow( r,funval, beta )
X=r.^beta;
alpha=X\funval;
figure
hold on
scatter(r,funval)
plot(r,alpha*r.^beta)
drawnow
end

function [ out ] = powvar(alpha,r,beta)
out=alpha.*r.^beta;
end

function [ r ] = rcalc( invec )

r= sqrt( invec(:,1).^2 +invec(:,2).^2);

end
