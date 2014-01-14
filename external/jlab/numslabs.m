function[n]=numslabs(A,B)
% NUMSLABS  Counts the number of 'slabs' of one array relative to another.
%  
%   NUMSLABS(A,B) counts the number of size SIZE(B) elements of A, 
%   that is,  it returns LENGTH(A(:))./LENGTH(B(:)).  The size of A 
%   must be 'compatible' with that of B.
%  
%   See also ISCOMPAT.  
%
%   'numslabs --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2006 J.M. Lilly --- type 'help jlab_license' for details    

if strcmp(A, '--t')
  numslabs_test,return
end

if ~iscompat(A,B)
   error('The size of A must be compatible with that of B.  See ISCOMPAT.')
end

n=numel(A)./numel(B);

function[]=numslabs_test(A,B)
x=[1 2; 3 4];y(:,:,2)=x;

b(1)=numslabs(y,x);
b(2)=numslabs(x,x);
b(3)=numslabs(y,1);
b(4)=numslabs(x,1);

bans=[2 1 8 4];

reporttest('NUMSLABS 2D and 3D cases',aresame(b,bans))

