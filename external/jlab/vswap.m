function[varargout]=vswap(varargin)
%VSWAP(X,A,B) replaces A with B in numeric array X
%
%   VSWAP(X,A,B) replaces A with B in numeric array X.  A and B may be
%   numbers, NAN, +/- INF, or NAN+SQRT(-1)*NAN.
%
%   [Y1,Y2,...YN]=VSWAP(X1,X2,...XN,A,B) also works.
%
%   VSWAP(X1,X2,...XN,A,B); with no output arguments overwrites the 
%   original input variables.    
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2001, 2004 J.M. Lilly --- type 'help jlab_license' for details  

if strcmp(varargin{1}, '--t')
  vswap_test,return
end
 
  
a=varargin{end-1};
b=varargin{end};

for i=1:length(varargin)-2
  x=varargin{i};
  varargout{i}=swapnum1(x,a,b);
end

eval(to_overwrite(nargin-2))  


function[x]=swapnum1(x,a,b)
    
    
if isfinite(a)
%   if a==0
%       if ~isreal(x)
%          a=0+sqrt(-1)*0;
%       end
%   end
  index=find(x==a);
else
  if isnan(a)
    index=find(isnan(x));
  elseif isinf(a)
    if a>0
        index=find(isinf(x)&x>0);
    else
        index=find(isinf(x)&x<0);
    end
  elseif isnan(real(a)) && isnan(imag(a))
    index=find(isnan(real(x))&isnan(imag(x)));
  end
end

if ~isempty(index)
    if allall(x==0|x==1)
       %Matlab, apparently, won't let you put NANs into a boolean array
       x=x+1;
       x(index)=b;
       x=x-1;
   else
       x(index)=b;
   end
end

function[]=vswap_test
x=(1:10);
ans1=[2 (2:10)];
reporttest('VSWAP num case', aresame(vswap(x,1,2),ans1))

x=[nan (1:10)];
ans1=(0:10);
reporttest('VSWAP nan case', aresame(vswap(x,nan,0),ans1))

x=[nan*(1+sqrt(-1)) (1:10)];
ans1=(0:10);
reporttest('VSWAP complex nan case', aresame(vswap(x,nan+sqrt(-1)*nan,0),ans1))
