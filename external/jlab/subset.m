function[varargout]=subset(varargin)
% SUBSET Extract subset of A given B.  
%
%   SUBSET(A,B), where A and B are arrays of the same size, and B is a
%   boolean array, returns the elements of A for which the
%   corresponding element of B is true.  SUBSET(A,B) is then
%   equivalent to A(FIND(B)). The output is a column array of length
%   LENGTH(FIND(B)).
%
%   If A is an MxN matrix and B is a length M column vector, SUBSET
%   returns the rows of A for which B is true.  The output will have
%   LENGTH(FIND(B)) rows and N columns. 
%
%   This behavior generalizes to higher dimensions if the size of A is
%   compatible with the size of B, as defined in ISCOMPAT. In this
%   case SUBSET will extract those "slabs" of A for which B is true.
%   The output will be an array with size 
%     [ LENGTH(FIND(B)) SIZE(A,ND(B)+1) ....  SIZE(A,ND(A))].
%  
%   SUBSET(A1,A2,...AN,B) applies itself to each AI given B. 
%
%   See also LOOKUP, ISMEMB
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003, 2004 J.M. Lilly --- type 'help jlab_license' for details      
  
if strcmp(varargin{1}, '--t')
 subset_test,return
end


  
n=nargin;
b=varargin{n};
x=varargin(1:n-1);

index=find(b);
lfb=length(index);

y=x;

for i=1:n-1
   a=x{i};
   
   sa=size(a);
 
   la=numel(a);
   lb=numel(b);
   
   nda=nd(a);  %Matlab's NDIMS is not smart.   
   ndb=nd(b);
   
   if ~iscompat(a,b)
      error('The size of A must be compatible with the size of B.')
   end  
   
   if nda==ndb
     y{i}=a(index);
   else    
     newsizea=[lb la./lb];
     a=reshape(a,newsizea);
     a=a(index,:);
     
     newsizea=[lfb sa(ndb+1:end)];
     a=reshape(a,newsizea);
     y{i}=a;
   end  
       
   
end
varargout=y;

  

function[]=subset_test
  
x=(1:10)';  
q=[x x];     
bool=0*q;
bool([5 15:18])=1;
r(:,:,1)=q;
r(:,:,2)=q;
r(:,:,3)=q;
x=subset(q,bool)';
b(1)=aresame(x,[5 5 6 7 8]);
%disp('Should be 5 5 6 7 8')
x=subset(r,bool);
b(2)=aresame(x,osum([5 5 6 7 8]',0*[1 1 1]'));
%disp('Should be three columns of 5 5 6 7 8')
r2=reshape(r,length(bool(:)),3);
x=subset(r2,bool(:));
b(3)=aresame(x,osum([5 5 6 7 8]',0*[1 1 1]'));
%disp('Should be three columns of 5 5 6 7 8')
reporttest('SUBSET',all(b))

