function [x0s,s,id,l,k,k0]=cokri2(x,x0,id,model,c,sv,itype,avg,ng)
%
% COKRI2 : This function is called from COKRI. The description for input and
%          output is given in COKRI. The only new variables are 'k0' which is
%          the right member matrix of the cokriging system and 'ng' which is
%          the total number of points for block discretization.
% Author: D. Marcotte
% Version 2.0  97/aug/14 
% External subroutines: trans, means

x0s=[];s=[];l=[];k=[];k0=[];nc=0;

% here we define the equations for the various covariograms. Any new model
% can be added here.

Gam=['h==0                                              '; %nugget
  'exp(-h)                                           '; %exponential
  'exp(-(h).^2)                                      '; %gaussian
  '1-(1.5*min(h,1)/1-.5*(min(h,1)/1).^3)             '; %spherical
  '1-h                                               '; %linear
  '1-3*min(h,1).^2+2*min(h,1).^3                     '; %modele Trochu
  '(h.^2).*log(max(h,eps))                           '];%spline palque mince

% definition of some constants

[n,t]=size(x);
[rp,p]=size(c);
r=rp/p;
[m,d]=size(x0);

% if no samples found in the search radius, return NaN

if n==0,
  x0s=NaN*ones(m/ng,p);
  s=NaN*ones(m/ng,p);
  return
end
cx=[x(:,1:d);x0];

% calculation of left covariance matrix K and right covariance matrix K0

k=zeros(n*p,(n+m)*p);
for i=1:r,
  
  % calculation of matrix of reduced rotated distances H
  
  [t]=trans(cx,model,i);
  t=t*t';
  h=sqrt(-2*t+diag(t)*ones(1,n+m)+ones(n+m,1)*diag(t)');
  h=h(1:n,:);
  ji=(i-1)*p+1; js=i*p ;
  
  % evaluation of the current basic structure
  g=eval(Gam(model(i,1),:));
  k=k+kron(g,c(ji:js,:));
end
k0=k(:,n*p+1:(n+m)*p); k=k(:,1:n*p);

% constraints are added according to cokriging type

if itype==99,
  
  % the system does not have to be solved
  return
end

if itype==1.5;
  
  % krigeage oridinaire avec 1 seule contrainte sur la var. auxiliaire
'arrêt'  
end

  
  if itype==2,
    
    % cokriging with one non-bias condition (Isaaks and Srivastava, 1990, p.410)
    
    k=[k,ones(n*p,1);ones(1,n*p),0]; k0=[k0;ones(1,m*p)]; nc=1;
  elseif itype>=3,
    
    % ordinary cokriging (Myers, Math. Geol, 1982)
    
    t=kron(ones(1,n),eye(p));
    k=[k,t';t,zeros(p,p)];
    k0=[k0;kron(ones(1,m),eye(p))]; nc=p;
    
    %   % cokriging with one non-bias condition in the z direction
    
    if itype == 3.5,
      
      t=kron(cx(1:n,d),eye(p));
      k=[k,[t;zeros(p,p)];[t',zeros(p,p+p)]];
      t=kron(cx(n+1:n+m,d)',eye(p));
      k0=[k0;t];
      nc=nc+p;
    end;
    
    if itype >=4,
      
      % universal cokriging ; linear drift constraints
      nca=p*d;
      t=kron(cx(1:n,:),eye(p));
      k=[k,[t;zeros(p,nca)];[t',zeros(nca,nc+nca)]];
      t=kron(cx(n+1:n+m,:)',eye(p));
      k0=[k0;t];
      nc=nc+nca;
    end;
    if itype==5,
      
      % universal cokriging ; quadratic drift constraints
      
      nca=p*d*(d+1)/2;
      cx2=[];
      for i=1:d,
        for j=i:d,
          cx2=[cx2,[cx(:,i).*cx(:,j)]];
        end
      end
      t=kron(cx2(1:n,:),eye(p));
      k=[k,[t;zeros(nc,nca)];[t',zeros(nca,nc+nca)]];
      t=kron(cx2(n+1:n+m,:)',eye(p));
      k0=[k0;t];
      nc=nc+nca;
    end
  end
  
  % columns of k0 are summed up (if necessary) for block cokriging
  
  m=m/ng;
  t=[];
  for i=1:m,
    for ip=1:p,
      j=ng*p*(i-1)+ip;
      t=[t,means(k0(:,j:p:i*ng*p)')'];
    end
  end
  k0=t;
  
  t=x(:,d+1:d+p);
  if itype<3,
    
    % if simple cokriging or cokriging with one non bias condition, the means
    % are substracted
    
    t=(t-ones(n,1)*avg)';
  else
    t=t';
  end
  
  % removal of lines and columns in K and K0 corresponding to missing values
  
  z=zeros(n*p,1);
  z(:)=t;
  iz=~isnan(z);
  iz2=[iz;ones(nc,1)]~=0;
  nz=sum(iz);
  
  % if no samples left, return NaN
  
  if nz==0,
    x0s=nan;
    s=nan;
    return;
  else
    k=k(iz2,iz2');
    k0=k0(iz2,:);   
    id=id(iz,:);
    
    % check there is at least one sample for each variable
    
    [n1,m1]=size(k);
    for j=n1:-1:n1-nc+1;
      if k(j,:)==zeros(1,m1)
        k=k(1:j-1,1:j-1);
        k0=k0(1:j-1,:);
        m1=m1-1;
      end
    end
    
    
    % solution of the cokriging system by gauss elimination
    
    l=inv(k)*k0;
    
    % calculation of cokriging estimates
    
    t2=l(1:nz,:)'*z(iz);
    t=zeros(p,m);
    t(:)=t2;
    
    % if simple or cokriging with one constraint, means are added back
    
    if itype<3,
      t=t'+ones(m,1)*avg;
    else
      t=t';
    end
    x0s=t;
    
    % calculation of cokriging variances
    
    s=kron(ones(m,1),sv);
    t=zeros(p,m);
    t(:)=diag(l'*k0);
    s=s-t';
  end
  
  