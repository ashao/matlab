function[varargout]=vadd(varargin)
%VADD   Vector-matrix addition without "dimensional" hassle.
%	
%   N = VADD(M,V,DIM) replicates the vector V along dimension DIM and 
%   adds the resulting matrix to M.  
%
%   VADD is essentially a shorthand for the commands
%
%     N=VADD(M,V,2)    <==>    N=M+V*ones(size(M(1,:))) 
%     N=VADD(M,V,1)    <==>    N=M+ones(size(M(:,1)))*V
%
%   where the first adds column vector V to each column of M, and the
%   second adds row vector V to each row of M.
%
%   [N1,N2,...]= VADD(M1,M2,... V,DIM) also works.
%
%   VADD(M1,M2,... V,DIM) where the MI are named matrices overwrites
%   these with the output matrices.
%  
%   See also OSUM, OPROD, ODOT
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details    
  
if strcmp(varargin{1}, '--t')
  vadd_test,return
end

dim=varargin{end};
vect=varargin{end-1};
for i=1:length(varargin)-2
  varargout{i}=vadd1(varargin{i},vect,dim);
end
  
eval(to_overwrite(nargin-2))

function[y]=vadd1(x,vect,dim)
n=size(x,dim);
vect=vrep(vect,n,dim);
vsize(x,vect);
y=x+vect;


function[]=vadd_test
x1=[1 1; 2 2];
v1=[1 2];
ans1=x1+ones(size(x1(:,1)))*v1;
vadd(x1,v1,1);
reporttest('VADD row vector', aresame(x1,ans1))

x1=[1 1; 2 2];
v1=[1 2]';
ans1=x1+v1*ones(size(x1(1,:)));
vadd(x1,v1,2);
reporttest('VADD col vector', aresame(x1,ans1))

x1(:,:,1)=[1 1; 2 2];
x1(:,:,2)=[1 1; 2 2];
v1=[1 1;1 1];
ans1=x1+1;
vadd(x1,v1,3);
reporttest('VADD 3-d case', aresame(x1,ans1))

