function[varargout]=vcellcat(varargin)
%VCELLCAT  Concatenates cell arrays of column vectors.
%
%   Y=VCELLCAT(X), where X is a cell array containing column vectors
%   of arbitrary length, concatenates all the column vectors together
%   and return the result in column vector Y.
%  
%   [Y1,Y2,...YN]=VCELLCAT(X1,X2,...XN) concatenates each of the cell
%   arrays Xi into the respective column vector Yi.
%
%   It is permissible that some of the XN be empty.
%
%   VCELLCAT(X1,X2,...XN); with no output arguments overwrites the
%   original input variables.  
%
%   See also VTOOLS, CELL2COL.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2009 J.M. Lilly --- type 'help jlab_license' for details        

if strcmp(varargin{1},'--t')  
   vcellcat_test;return
end
   
for i=1:nargin
    x=varargin{i};
    L=sum(cellsize(x,1));
    a=1;
    b=0;
    if L>0
        y=zeros(L,size(x{find(cellsize(x,1)~=0,1)},2));
        for j=1:length(x)
               a=b+1;
               b=a+size(x{j},1)-1;
               if ~isempty(x{j})
                   y(a:b,:)=x{j};
               end
        end
    else 
        y=[];
    end
    varargout{i}=y;
end

eval(to_overwrite(nargin));

function[]=vcellcat_test
x{1}=[1 1]';
x{2}=3;
x{3}=[2 2 2]';
y=x;
z=[1 1 3 2 2 2]';
vcellcat(x,y);
reporttest('VCELLCAT column vector', aresame(x,z) && aresame(y,z))

clear x
x{1}=[1 1]';
x{2}=3;
x{3}=[];
x{4}=[2 2 2]';
y=x;
z=[1 1 3 2 2 2]';
vcellcat(x,y);
reporttest('VCELLCAT column vector with empty', aresame(x,z) && aresame(y,z))

clear x
x{1}=[];
x{2}=[];
x{3}=[];
x{4}=[2 2 2]';
y=x;
z=[2 2 2]';
vcellcat(x,y);
reporttest('VCELLCAT column vector leading empties', aresame(x,z) && aresame(y,z))

clear x
x{1}=[1 1; 2 2];
x{2}=[3 3];
x{3}=[2 2; 3 3 ; 4 4];
y=x;
z=[x{1};x{2};x{3}];
vcellcat(x);

reporttest('VCELLCAT two-column array', aresame(x,z))