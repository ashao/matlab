function[xi,yi,ti]=gridadvect(u,v,x1,y1,t,xo,yo,dt)
%GRIDADVECT  Advect a set of particles given an evolving velocity field.
%
%   [XI,YI,TI]=GRIDADVECT(U,V,X,Y,T,XO,YO,DT) advects a set of particles
%   with initial positions (XO,YO) using the velocity field (U,V) with a
%   finer time step DT than given by the velocity snapshots.
%
%   U and V are 3-D arrays with the y-coordinate in *rows*, the x-
%   coordinate in *columns*, and time in "pages" along the third 
%   dimension.
%
%   X and Y are vectors specifying the physical dimensions of the grid,
%   in kilometers.  T is a vector specifying the time axis, in days.
%   XO and YO are also in kilometers, and U and V are in centimeters
%   per second.
%
%   The integration is carried out with time step DT, in days.
%
%   XI and YI give the evolving positions of the particles at times TI.
%   XI and YI are oriented with time in rows and particles in columns.
%
%   Note that owing to the constraints of the interpolation, TI will 
%   being at the second snapshot, time T(2), and will proceed with time
%   step DT until time T(END-2).   
%   __________________________________________________________________
%
%   Interpolation algorithm
%
%   GRIDADVECT uses a double cubic interpolation algorithm, with basin 
%   escape being explicitly prevented.
%
%   Firstly, interpolation of the velocities at the four corners of each 
%   grid box containing a particle are cubicly interpolated in time using 
%   two time steps ahead and two behind.  
%
%   Within each box, the four velocities are then interpolated to the 
%   particle position, again with a cubic algorithm.  
%
%   See CUBENITERP for details.
%   __________________________________________________________________
%
%   Usage:  [xi,yi,ti]=gridadvect(u,v,x,y,t,xo,yo,dt);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007 J.M. Lilly --- type 'help jlab_license' for details


t=t(:);
x1=x1(:);
y1=y1(:);
if length(x1)~=size(u,2)||length(x1)~=size(v,2)
    error('Length of X must be the same as the number of columns in U and V.')
end
if length(y1)~=size(u,1)||length(y1)~=size(v,1)
    error('Length of Y must be the same as the number of rows in U and V.')
end
if length(t)~=size(u,3)||length(t)~=size(v,3)
    error('Length of T must be the same as the number of pages in U and V.')
end
etao=xo+sqrt(-1)*yo;

ti=(t(2):dt:maxmax(t(1:end-2))-dt)';
N=length(ti);
etai=oprod([1;zeros(N-1,1)],etao(:));

xmax=maxmax(x1);
xmin=minmin(x1);
ymax=maxmax(y1);
ymin=minmin(y1);
dx=abs(x1(2)-x1(1));

%Grid sides as complex-valued array
x=oprod(x1(:),ones(size(etai,2),1));
y=oprod(y1(:),ones(size(etai,2),1));
eta=x+sqrt(-1)*y;

%Grid midpoints as complex-valued array
xmid=frac(1,2)*(x1(1:end-1)+x1(2:end));
ymid=frac(1,2)*(y1(1:end-1)+y1(2:end));
xmid=oprod(xmid(:),ones(size(etai,2),1));
ymid=oprod(ymid(:),ones(size(etai,2),1));
etamid=xmid+sqrt(-1)*ymid;

%Velocity as complex-valued array
zeta=(u+sqrt(-1)*v)*dt*0.864;



for i=2:N
   i./N*100;
   etai(i,:)=gridadvect1(eta,zeta,etamid,etai(i-1,:),t,ti(i-1));
   etai(i,:)=prevent_basin_escape(etai(i,:),xmin,xmax,ymin,ymax,dx);
end
  
xi=real(etai);
yi=imag(etai);

function[etai]=gridadvect1(eta,zeta,etamid,etao,t,ti) 

eta1=oprod(ones(size(etamid(:,1))),etao(:));

[mm,ii]=min(abs(imag(etamid-eta1)));
[mm,jj]=min(abs(real(etamid-eta1)));
[mm,kk]=min(abs(ti-t));

imat(1,:,1)=ii;
imat(1,:,2)=ii;
imat(2,:,1)=ii+1;
imat(2,:,2)=ii+1;

jmat(1,:,1)=jj;
jmat(1,:,2)=jj+1;
jmat(2,:,1)=jj;
jmat(2,:,2)=jj+1;

kmat=kk+0*imat;
index1=sub2ind(size(zeta),imat,jmat,kmat-1);
index2=sub2ind(size(zeta),imat,jmat,kmat);
index3=sub2ind(size(zeta),imat,jmat,kmat+1);
index4=sub2ind(size(zeta),imat,jmat,kmat+2);

%Time interpolation
zetamat=cubeinterp(t(kk-1),t(kk),t(kk+1),t(kk+2),zeta(index1),zeta(index2),zeta(index3),zeta(index4),ti); 

etamat=zeros(size(imat));
etamat(1:end)=real(eta(jmat,1))+sqrt(-1)*imag(eta(imat,1));

%Space interpolation
eta1=etamat(1,:,1);
eta2=etamat(1,:,2);
eta3=etamat(2,:,1);
eta4=etamat(2,:,2);

zeta1=zetamat(1,:,1);
zeta2=zetamat(1,:,2);
zeta3=zetamat(2,:,1);
zeta4=zetamat(2,:,2);


zeta=cubeinterp(eta1,eta2,eta3,eta4,zeta1,zeta2,zeta3,zeta4,etao);  
etai=etao+zeta;

function[eta]=prevent_basin_escape(eta,xmin,xmax,ymin,ymax,dx) 
tol=1e-3*dx;
index=find(real(eta)>xmax);
if ~isempty(index)
    eta(index)=sqrt(-1)*imag(eta(index))+xmax-tol;
end
index=find(real(eta)<xmin);
if ~isempty(index)
    eta(index)=sqrt(-1)*imag(eta(index))+xmin+tol;
end

index=find(imag(eta)>ymax);
if ~isempty(index)
    eta(index)=real(eta(index))+sqrt(-1)*(ymax-tol);
end

index=find(imag(eta)<ymin);
if ~isempty(index)
    eta(index)=real(eta(index))+sqrt(-1)*(ymin+tol);
end


%zeta=(zeta1+zeta2+zeta3+zeta4)/4;




function[]=gridadvect_test



