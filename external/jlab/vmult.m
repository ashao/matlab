function[varargout]=vmult(varargin)
%VMULT   Vector-matrix multiplication without "dimensional" hassle.
%	
%   N = VMULT(M,V,DIM) replicates the vector V along dimension DIM and 
%   multiplies the resulting matrix to M.  '
%
%   VMULT is essentially a shorthand for the commands
%
%     N=VMULT(M,V,2)    <==>    N=M.*(V*ones(size(M(1,:)))) 
%     N=VMULT(M,V,1)    <==>    N=M.*(ones(size(M(:,1)))*V)
%
%   where the first multiplies the column vector V along each column of 
%   M, and the second multiplies the row vector V along each row of M.
%
%   [N1,N2,...]= VMULT(M1,M2,... V,DIM) also works.
%
%   VMULT(M1,M2,... V,DIM) where the MI are named matrices overwrites
%   these with the output matrices.
%  
%   See also OSUM, OPROD, ODOT
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(varargin{1}, '--t')
  vmult_test,return
end

dim=varargin{end};
vect=varargin{end-1};

for i=1:length(varargin)-2
  varargout{i}=vmult1(varargin{i},vect,dim);
end
  
eval(to_overwrite(nargin-2))

function[y]=vmult1(x,vect,dim)
n=size(x,dim);
vect=vrep(vect,n,dim);
vsize(x,vect);
y=x.*vect;


function[]=vmult_test
x1=[1 1; 2 2];
v1=[1 2];
ans1=x1.*(ones(size(x1(:,1)))*v1);
vmult(x1,v1,1);
reporttest('VMULT row vector', aresame(x1,ans1))

x1=[1 1; 2 2];
v1=[1 2]';
ans1=x1.*(v1*ones(size(x1(1,:))));
vmult(x1,v1,2);
reporttest('VMULT col vector', aresame(x1,ans1))

x1(:,:,1)=[1 1; 2 2];
x1(:,:,2)=[1 1; 2 2];
v1=[1 1;1 1];
ans1=x1;
vmult(x1,v1,3);
reporttest('VMULT 3-d case', aresame(x1,ans1))

