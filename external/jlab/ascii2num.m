function[varargout]=ascii2num(varargin)
%ASCII2NUM  Convert ASCII values for numbers into numeric values.
%
%   N=ASCII2NUM(A) converts a matrix A of ASCII charcter codes, 
%   reprenting numbers, into a vector of numeric values N. 
%
%   Thus, CHAR(A) would be a string array oriented as follows
%
%         CHAR(A) = [ '   -6.1512'; 
%                     '20103.1321'; 
%                     '- 010.1232'].  
%
%   CHAR(A) may include numbers, minus signs, and decimal points.  
%   If a decimal point exists, it must be in the same column for 
%   all rows. Minus signs may be present in any row.  
%   
%   [N1,N2,... NN]=ASCII2NUM(A1,A2,... AN) also works.
%
%   For large A, ASCII2NUM is about 1-2 orders of magnitude faster
%   than the built-in Matlab alternative, STR2NUM(CHAR(A)).
%
%   If ASCII2NUM detects that the decimal points are not lined up,
%   it switches to the slower method.
%
%   ASCII2NUM is useful in reading in large files of ASCII data.
%
%   'ascii2num --t' runs a test.
%  
%   Usage:  n=ascii2num(a);
%           [n1,n2,n3]=ascii2num(a1,a2,a3);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2010 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    ascii2num_test,return
end
 

for i=1:nargin
    x=varargin{i};
    varargout{i}=ascii2num1(x);
end

function[y]=ascii2num1(x)


%Find and remove all minus signs
index=find(x==real('-'));
sgn=ones(size(x(:,1)));
if ~isempty(index)
    x(index)=real(' ');
    [ii,jj]=ind2sub(size(x),index);
    sgn(ii)=-1;
    bool=false(size(x));
    bool(index)=1;
    if anyany(sum(bool,2)>1)
        error('More than one minus sign in at least one row.')
   end
end

vswap(x,real(' '), real('0'));

%Does there exist a decimal place?
jd=find(x(1,:)==real('.'));

usechar=0;
if isempty(jd)
    jd=size(x,2)+1;
elseif length(jd)>1
    error('More than one decimal point in first row.')
elseif length(jd)==1
    b=all(x(:,jd)==real('.'));
    if ~b
        disp('Decimal point does not exist in all rows; switching to STR2NUM(CHAR).')
        usechar=1;
    end    
    
    if any(find(x(:,[1:jd-1 jd+1:end])==real('.')))
        size(x(:,[1:jd-1 jd+1:end])==real('.'))
        (x(:,[1:jd-1 jd+1:end])==real('.'))
        disp('Decimal points exist in more than one column; switching to STR2NUM(CHAR).')
        usechar=1;
    end
end

if usechar
   % char(x)
    %size(x),size(char(x)),size(str2num(char(x)))
    y=str2num(char(x));
else
    %JD is now column of tenths place
    x=x(:,[1:jd-1 jd+1:end]);
    tens=10.^(jd-2:-1:jd-size(x,2)-1);
    y=(x-48)*tens';
end

%size(sgn)
%size(y)
y=y.*sgn;


function[]=ascii2num_test

tol=1e-10;
clear x
x(1,:)='2006.152';
x(2,:)='2007.132';
x=real(x);

reporttest('ASCII2NUM no minus signs',aresame(ascii2num(x),str2num(char(x)),tol))

clear x
x(1,:)='  - 6.152';
x(2,:)='-2007.132';
x=real(x);

reporttest('ASCII2NUM with minus signs',aresame(ascii2num(x),str2num(char(x)),tol))

N=10000;
x=round(rand(N,1)*1000*1000)/1000;
x=vnum2str(x,-3);
x=real(x);

tic;y1=str2num(char(x));et1=toc;
tic;y2=ascii2num(x);et2=toc;
reporttest('ASCII2NUM random',aresame(y1,y2,tol)) 
warning('off','MATLAB:divideByZero')
disp(['ASCII2NUM(A) was ' num2str(et1./et2,3) ' times faster than STR2NUM(CHAR(A)).'])
warning('on','MATLAB:divideByZero')


