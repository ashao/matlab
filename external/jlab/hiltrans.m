function[y]=hiltrans(x)
%HILTRANS  Hilbert transform.
%
%   Y=HILTRANS(X) returns the Hilbert transform of column vector X.
%
%   If SIZE(X,2)>1, HILTRANS takes the Hilbert transform along columns.
%
%   'hiltrans --f' makes a sample figure
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004 J.M. Lilly --- type 'help jlab_license' for details        
  
if strcmp(x,'--f')
  hiltrans_fig;return
end

if isreal(x)
  bool=1;
else
  bool=0;
end

N=size(x,1);
X=fft(x);
sgnom=ones(N,1);
%sgnom(1)=0;
index=(1:N);
sgnom(index-1>N/2)=-1;
d=-sqrt(-1)*(sgnom);
d=oprod(d,1+zeros(size(x(1,:)))');
y=ifft(X.*d);

%Take real part if real vector was input
if bool
  y=real(y);
end

function[]=hiltrans_fig
n=256;
t=(-n:n)';
warning('off','MATLAB:divideByZero')
x=1./t;
warning('on','MATLAB:divideByZero')
index=find(~isfinite(x));
if ~isempty(index)
  x(index)=0;
end

y=hiltrans(x);
figure,plot(cumsum(y))
disp('This is a plot of cumsum(hiltrans(1/x))')
disp('Hilbert transform of 1/x should be -pi*delta function')

