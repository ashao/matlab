function[ii,jj]=lookup(a,b,tol)
% LOOKUP Locate elements of one array within another
%  
%   [II,JJ]=LOOKUP(A,B) looks up the elements of B within A and
%   returns two arrays II,JJ such that A(II)==B(JJ). Note that A must
%   be a set, i.e., it must have no repeated elements, which B can
%   have any number of elements.  II and JJ need not have the same
%   length as B(:), since some elements of B may not be members of A. 
%
%   LOOKUP uses two different algorithms.  When A is a 'wavenumber'
%   grid, i.e. it is a square matrix with regularly spaced
%   complex-valued entries such as that created by WAVEGRID, then
%   there exists a unique relationship between an entry's value and
%   its position within A.  This relationship is diagnosed from A and
%   then applied to the values of B, resulting in a fast computation
%   of II and JJ.
%  
%   If A is not a wavenumber matrix, a slower method is used.  which
%   tests for equality between A and B with a numerical tolerance
%   TOL=1e-10.  [II,JJ]=LOOKUP(A,B,TOL) specifies TOL.
%
%   See also SUBSET, ISMEMB, WAVEGRID, ISWAVEGRID
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2003, 2004 J.M. Lilly --- type 'help jlab_license' for details        

if strcmp(a, '--t')
  lookup_test,return
end

if nargin==2
  tol=1e-10;
end

ii=[];
jj=[];

bool=iswavegrid(a);

%Use a split algorithm
if bool
     %If a grid, we know [ii,jj] at once
   N=sqrt(length(a(:)));
   dk=abs(a(2)-a(1));
   [ii,jj,index1]=k2sub(b,dk,N);
       
%    %strip out of range values
%    ii=ii(:);jj=jj(:);
%    index1=find(ii>0 && ii<N+1 && jj>0 && jj<N+1);
     if ~isempty(index1)        
%        ii=ii(index1);jj=jj(index1);
         index=sub2ind([N N],ii,jj);
         jj=index1;
         ii=index;
    end 
   
elseif ~bool  %If not a grid, have to use a slow method
   if ~isset(a)
      error('A must be a set.')
   end
   a=a(:);
   b=b(:);
   [ii,jj]=find(abs(osum(a,-b))<=tol);
end

function[]=lookup_test
 
a=[1 3 5 7 9];
b=[1 2 3 4 5 6 7];
[ii,jj]=lookup(a,b);
reporttest('LOOKUP for sample set', all(ii==(1:4)'&jj==[1 3 5 7]'))

kg=wavegrid(pi,5);       
[ii1,jj1]=lookup(kg,kg);
[ii2,jj2]=lookup(kg(:),kg(:));
reporttest('LOOKUP for 5x5 WAVEGRID, first subscripts match',all(ii1==ii2))
reporttest('LOOKUP for 5x5 WAVEGRID, second subscripts match',all(jj1==jj2))
reporttest('LOOKUP for 5x5 WAVEGRID, lookup correct',all(kg(ii2)==kg(jj2)))

kg=wavegrid(pi,4);       
[ii1,jj1]=lookup(kg,kg);
[ii2,jj2]=lookup(kg(:),kg(:));
reporttest('LOOKUP for 4x4 WAVEGRID, first subscripts match',all(ii1==ii2))
reporttest('LOOKUP for 4x4 WAVEGRID, second subscripts match',all(jj1==jj2))
reporttest('LOOKUP for 4x4 WAVEGRID, lookup correct',all(kg(ii2)==kg(jj2)))

index=randperm(length(kg(:)))';
index=index(1:floor(end/2));
temp=kg(index);
[ii1,jj1]=lookup(kg,kg(index));
[ii2,jj2]=lookup(kg(:),kg(index));

reporttest('LOOKUP for scrambled 5x5 WAVEGRID, first subscripts match',all(ii1==ii2))
reporttest('LOOKUP for scrambled 5x5 WAVEGRID, second subscripts match',all(jj1==jj2))
reporttest('LOOKUP for 5x5 WAVEGRID, lookup correct',all(kg(ii2)==temp(jj2)))
  

% k=[3.3629 + 1.4622i ; 3.3629 - 1.4622i ; 6.7257]; 
% a=[0.5900;    0.3400;    0.3400];
% phi=[  0;         0;   -1.5708];
% x1=a.*rot(phi);
% [k2,x2]=findk2(k,k,x1);
% [ii,jj]=lookup(k(:),k2(:));
% reporttest('LOOKUP for wave triad, wavenumbers match',all(k(ii)==k2(jj)))
% reporttest('LOOKUP for wave triad, coefficients match',all(x1(ii)==x2(jj)))
