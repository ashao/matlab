function[bool]=islargest(i1,i2,i3)
%ISLARGEST  True for the largest-magnitude element at each index location.
%  
%   BOOL=ISLARGEST(II,X) where II is an index into row locations at 
%   which values X occur, retuns an array BOOL which equals one for 
%   those elements of X which are the largest than any others having 
%   the same value of II, and zero otherwise.
%  
%   All three arrays BOOL, II, and X have the same size.
%
%   If X is complex-valued, its absolute value is taken.
%  
%   ISLARGEST is defined to return zero for NAN elements of X.
%   ___________________________________________________________________
%
%   Two indices
%
%   BOOL=ISLARGEST(II,JJ,X) where JJ is an array into column locations 
%   at which the values X occur, leads to BOOL equal to one for those  
%   elements of X which are the largest in magnitude for each particular 
%   II,JJ location.
%   ___________________________________________________________________
%
%   Usage: bool=islargest(ii,x);
%          bool=islargest(ii,jj,x);
%
%   'islargest --t' runs a test
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2005--2008 J.M. Lilly --- type 'help jlab_license' for details


if strcmp(i1,'--t')
  islargest_test;return
end  


if nargin==2
    bool=islargest1(i1,i2);
elseif nargin==3
    bool=islargest2(i1,i2,i3);
end


function[bool]=islargest1(ir,yr)
%Column vector case

yro=yr;

ir=ir(:);
yr=yr(:);

origorder=(1:length(ir));  
[ir,index]=sort(ir);
yr=yr(index);
origorder=origorder(index);

[L,ia,ib]=blocklen(ir);
ii=[1;find(diff(ia)~=1)];
vindex(ia,ib,ii,1);

bool=true(size(ir));
for i=1:length(ia)
  index=ia(i):ib(i);
  bool(index)=0;
  [m,ii]=max(abs(yr(index)));
  bool(index(ii))=1;
end

bool(origorder)=bool;

index=find(isnan(yro));
if ~isempty(index)
  bool(index)=0;
end
  
bool=reshape(bool,size(yro));


function[bool]=islargest2(ir,jr,yr)
%Matrix case

indexr=sub2ind([maxmax(ir) maxmax(jr)],ir,jr);
bool=islargest1(indexr,yr);


function[]=islargest_test

x=       [8 5]';  
t=       [4 4]';
boolans= [1 0]';
reporttest('ISLARGEST short case', aresame(islargest(t,x),boolans))

x=       [1 10 2 8 5 6 1 nan 9]';  
t=       [1 2  2 4 4 5 6 7   5]';
boolans= [1 1  0 1 0 0 1 0   1]';
reporttest('ISLARGEST column case', aresame(islargest(t,x),boolans))
reporttest('ISLARGEST row case', aresame(islargest(t',x'),boolans'))

x=       [1 10 2 8 5 6 1 nan 9 nan]';  
j=       [1 3  1 2 2 3 4 3   3 nan]';  
t=       [1 2  2 4 4 5 6 7   5 nan]';
boolans= [1 1  1 1 0 0 1 0   1 0]';
reporttest('ISLARGEST matrix column case', aresame(islargest(t,j,x),boolans))
reporttest('ISLARGEST matrix row case', aresame(islargest(t',j',x'),boolans'))



