function [x0s,s,sv,idout,l,k,k0]=cokri(x,x0,model,c,itype,avg,block,nd,ival,nk,rad,ntok)
%
% COKRI performs point or block cokriging in D dimensions (any integer)
%       of P variables (any integer) with a combination of K basic models
%       (any integer).
%
% Syntax:
% [x0s,s,sv,id,l]=cokri(x,x0,model,c,itype,avg,block,nd,ival,nk,rad,ntok)
%
% Input description:
%   x: The n x (p+d) data matrix. This data matrix can be imported from an
%      existing ascii file. Missing values are coded 'nan' (not-a-number).
%   x0: The m x d matrix of coordinates of points to estimate.
%   model: Each row of this matrix describes a different elementary structure.
%          The first column is a code for the model type, the d following
%          columns give the ranges along the different coordinates and the
%          subsequent columns give rotation angles (a maximum of three).
%	   For more details on how to specify rotations, type help trans.
%          The codes for the current models are:
%             1: nugget effect
%             2: exponential model
%             3: gaussian model
%             4: spherical model
%             5: linear model
%          Note: a linear model is specified by arbitrary ranges and a sill
%                such that sill/range gives the desired slope in the direction
%                condidered.
%   c: The (rp x p) coefficient matrix of the coregionalization model.
%      Position (i,j) in each submatrix of size p x p give the sill of the
%      elementary component for each cross-variogram (variogram) between
%      variable i and variable j.
%   itype: Code to indicate which type of cokriging is to be performed:
%             1:  simple cokriging
%             2:  ordinary cokriging with one nonbias condition
%                 (Isaaks and Srivastava).
%             3:  ordinary cokriging with p nonbias condition.
%             4:  universal cokriging with drift of order 1.
%             5:  universal cokriging with drift of order 2.
%             99: cokriging is not performed, only sv is computed.
%   block: Vector (1 x d), giving the size of the block to estimate;
%          any values when point cokriging is required.
%   nd: Vector (1 x d), giving the discretization grid for block cokriging;
%       put every element equal to 1 for point cokriging.
%   ival: Code for cross-validation.
%             0:  no cross-validation
%             1:  cross-validation is performed  by removing one variable at a
%                 time at a given location.
%             2:  cross-validation is performed by removing all variables at a
%                 given location.
%   nk: Number of nearest neighbors in x matrix to use in the cokriging
%       (this includes locations with missing values even if all variables
%       are missing).
%   rad: Search radius for neighbors.
%   ntok: Points in x0 will be kriged by groups of ntok grid points.
%         When ntok>1, the search will find the nk nearest samples within
%         distance rad from the current ntok grid points centroid.
%컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴컴
% Output description:
%
%   For the usual application, only x0s and s are required and the other
%   output matrices may be omitted.
%
%   x0s: m x (d+p) matrix of the m points (blocks) to estimate by the
%        d coordinates and p cokriged estimates.
%   s: m x (d+p) matrix of the m points (blocks) to estimate by the
%      d coordinates and the p cokriging variances.
%   sv: 1 x p vector of variances of points (blocks) in the universe.
%   id: (nk x p) x 2 matrix giving the identifiers of the lambda weights for
%       the last cokriging system solved.
%   l: ((nk x p) + nc) x (ntok x p) matrix with lambda weights and
%      Lagrange multipliers of the last cokriging system solved.
%
% Author: D. Marcotte
% Version 2.0  97/aug/14 (matlab4 and 5 compatible)
% External subroutines: cokri2, trans, means

%
% casesen off;

x0s=[];s=[];sv=[];id=[];l=[];k=[];k0=[];

% definition of some constants.
[m,d]=size(x0);

% check for cross-validation

if ival>=1,
   ntok=1;
   x0=x(:,1:d);
   nd=ones(1,d);
   [m,d]=size(x0);
end
[rp,p]=size(c);
[n,t]=size(x);
nk=min(nk,n);
ntok=min(ntok,m);
idp=[1:p]';
ng=prod(nd);

% compute point (ng=1) or block (ng>1) variance

t2=[];
for i=1:d,
   nl=prod(nd(1:i-1));
   nr=prod(nd(i+1:d));
   t=[.5*(1/nd(i)-1):1/nd(i):.5*(1-1/nd(i))]';
   t2=[t2,kron(ones(nl,1),kron(t,ones(nr,1)))];
end
grid=t2.*(ones(ng,1)*block);
t=[grid,zeros(ng,p)];

% for block cokriging, a double grid is created by shifting slightly the
% original grid to avoid the zero distance effect (Journel and Huijbregts, p.96).

if ng>1,
   grid=grid+ones(ng,1)*block/(ng*1e6);
end
[x0s,s,id,l,k,k0]=cokri2(t,grid,[],model,c,sv,99,avg,ng);

% sv contain the variance of points or blocks in the universe

for i=1:p,
   sv=[sv,means(means(k0(i:p:ng*p,i:p:ng*p))')];
end

% start cokriging

for i=1:ntok:m,
   nnx=min(m-i+1,ntok);
   ['kriging points #',num2str(i),' to ', num2str(i+nnx-1)]

   % sort x samples in increasing distance relatively to centroid of 'ntok'
   % points to krige

   centx0=ones(n,1)*means(x0(i:i+nnx-1,:));
   tx=[x(:,1:d)-centx0].*[x(:,1:d)-centx0]*ones(d,1);
   [tx,j]=sort(tx);

   % keep samples inside search radius; create an identifier of each sample
   % and variable (id)

   id=[];
   t=[];
   ii=1;
   tx=[tx;nan];
   while ii<=nk & tx(ii)<rad*rad,
      t=[t;x(j(ii),:)];
      id=[id;[ones(p,1)*j(ii),idp]];
      ii=ii+1;
   end
   t2=x0(i:i+nnx-1,:);

   %  if block cokriging discretize the block

   t2=kron(t2,ones(ng,1))-kron(ones(nnx,1),grid);

   % check for cross-validation

   if ival>=1,
      est=zeros(1,p);
      sest=zeros(1,p);

      % each variable is cokriged in its turn

      if ival==1,
         np=1;
      else
         np=p;
      end
      for ip=1:np:p,

         % because of the sort, the closest sample is the sample to
         % cross-validate and its value is in row 1 of t; a temporary vector
         % keeps the original values before performing cokriging.

         vtemp=t(1,d+ip:d+ip+np-1);
         t(1,d+ip:d+ip+np-1)=ones(1,np)*nan;
         [x0ss,ss,idout,l,k,k0]=cokri2(t,t2,id,model,c,sv,itype,avg,ng);
         est(ip:ip+np-1)=x0ss(ip:ip+np-1);
         sest(ip:ip+np-1)=ss(ip:ip+np-1);
         t(1,d+ip:d+ip+np-1)=vtemp;
      end
      x0s=[x0s;[t2,est]];
      s=[s;[t2,sest]];
   else
      [x0ss,ss,idout,l,k,k0]=cokri2(t,t2,id,model,c,sv,itype,avg,ng);
      x0s=[x0s;[x0(i:i+nnx-1,:),x0ss]];
      s=[s;[x0(i:i+nnx-1,:),ss]];
   end
%keyboard
end
