function[x]=cubeinterp(t1,t2,t3,t4,x1,x2,x3,x4,t)
%CUBEINTERP  Fast cubic interpolation for arbitrary-sized arrays.
%
%   XI=CUBEINTERP(T1,T2,T3,T4,X1,X2,X3,X4,TI) interpolates a function 
%   X within an independent variable T using a fast cubic algorithm.
%
%   The function values X1--X4 may be arrays of any size provided 
%   they are all the same size.  The "time" values T1--T4 may be 
%   scalars, or arrays of the same size as X1--X4.  
%
%   If T1--T4 and X1--X4 are non-scalar arrays, T should either be
%   a scalar or an array of the same size as the other arguments.
%
%   If T1--T4 and X1--X4 are all scalars, T may be of any size, and
%   XI will have the same size as T.
%
%   The function values X1--X4 may be real or complex.
%
%   CUBEINTERP uses the exact algebraic expressions to operate
%   looplessly on matrices, thus it is very fast.    
%   __________________________________________________________________
%
%   Algorithm details
%
%   Let X and T be related as
%
%              X = A*T.^3 + B*T.^2 + C*T + D
%
%   The coefficients A, B, C, and D are uniquely determined by four 
%   (T,X) pairs.  CUBEINTERP solves for these coefficients and uses 
%   the result to interpolate X to any other value of T.   
%   __________________________________________________________________
%
%   See also LININTERP, QUADINTERP.  
%
%   Usage:   x=cubeinterp(t1,t2,t3,t4,x1,x2,x3,x4,t);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(t1, '--t')
    cubeinterp_test,return
end

if ~aresame(size(x1),size(x2))||~aresame(size(x1),size(x3))||~aresame(size(x1),size(x4))
    error('The input arrays X1, X2, X3, X4 must all be the same size.')
end

if ~aresame(size(t1),size(t2))||~aresame(size(t1),size(t3))||~aresame(size(t1),size(t4))
    error('The input arrays T1, T2, T3, T4 must all be the same size.')
end

maxlent=max([length(t1(:)) length(t2(:)) length(t3(:)) length(t4(:))]);
maxlenx=max([length(x1(:)) length(x2(:)) length(x3(:)) length(x4(:))]);

if maxlent>1||maxlenx>1
    if ~aresame(size(t1),size(t))&&~aresame(size(x1),size(t))&&~isscalar(t)
        error('The time T must either be the same size as the other input arguments, or a scalar.')
    end
end

if 0
%To let time not be exact same size as matrices
if ~aresame(size(t1),size(x1))
    [t1mat,t2mat,t3mat,t4mat]=vzeros(size(x1));
    btime=0;
    for i=1:ndims(x1)
        if length(t1(:))==size(x1,i)
            if ~btime
                btime=1;
                for k=1:length(t1(:))
                    [t1mat,t2mat,t3mat,t4mat]=vindexinto(t1mat,t2mat,t3mat,t4mat,t1(k),t2(k),t3(k),t4(k),k,i);
                end
            end
        end
    end
    t1=t1mat;
    t2=t2mat;
    t3=t3mat;
    t4=t4mat;
end
end    
    
    

numa=frac(x1-x2,t1-t2)-frac(x2-x3,t2-t3)+frac(x3-x4,t3-t4)-frac(x4-x1,t4-t1);
dena=frac(t1.^3-t2.^3,t1-t2)-frac(t2.^3-t3.^3,t2-t3)+frac(t3.^3-t4.^3,t3-t4)-frac(t4.^3-t1.^3,t4-t1);
a=numa./dena;

numb=frac(x1-x2,t1-t2)-frac(x2-x3,t2-t3)-a.*frac(t1.^3-t2.^3,t1-t2)+a.*frac(t2.^3-t3.^3,t2-t3);
denb=(t1-t3);
b=numb./denb;

numc=(x1-x2)-a.*(t1.^3-t2.^3)-b.*(t1.^2-t2.^2);
denc=(t1-t2);
c=numc./denc;

d=x1-a.*t1.^3-b.*t1.^2-c.*t1;


x=a.*t.^3+b.*t.^2+c.*t+d;


function[]=cubeinterp_test

cubeinterp_test_real;
cubeinterp_test_complex;
cubeinterp_test_doubly_complex;

function[]=cubeinterp_test_real

a=1;
b=7;
c=-23;
d=18;

t1=2;
t2=4;
t3=8;
t4=11;

%t=(-20:.1:20)';
t=[t1 t2 t3 t4];
x=a.*t.^3+b.*t.^2+c.*t+d;

%ti=(-20:.1:20)';
%xi=cubeinterp(t(1),t(2),t(3),t(4),x(1),x(2),x(3),x(4),ti);
%plot(ti,xi),hold on,plot(t,x,'go')

N=100;
bool=false(N,1);
for i=1:100
   ti=randn(1)*10;
   xi=cubeinterp(t(1),t(2),t(3),t(4),x(1),x(2),x(3),x(4),ti);
   xi2=a.*ti.^3+b.*ti.^2+c.*ti+d;
   bool(i)=aresame(xi,xi2,1e-10);
end

reporttest('CUBEINTERP real scalar case',all(bool(i)))


function[]=cubeinterp_test_complex
a=1+sqrt(-1)*4;
b=7+sqrt(-1)*2;
c=-23+sqrt(-1)*12;
d=18-sqrt(-1)*7;

t1=2;
t2=4;
t3=8;
t4=11;

%t=(-20:.1:20)';
t=[t1 t2 t3 t4];
x=a.*t.^3+b.*t.^2+c.*t+d;

N=100;
bool=false(N,1);
for i=1:100
   ti=randn(1)*10;
   xi=cubeinterp(t(1),t(2),t(3),t(4),x(1),x(2),x(3),x(4),ti);
   xi2=a.*ti.^3+b.*ti.^2+c.*ti+d;
   bool(i)=aresame(xi,xi2,1e-10);
end


reporttest('CUBEINTERP complex scalar case',all(bool(i)))



function[]=cubeinterp_test_doubly_complex
a=1+sqrt(-1)*4;
b=7+sqrt(-1)*2;
c=-23+sqrt(-1)*12;
d=18-sqrt(-1)*7;

t1=2+sqrt(-1);
t2=4-sqrt(-1);
t3=8-3*sqrt(-1);
t4=11+17*sqrt(-1);

%t=(-20:.1:20)';
t=[t1 t2 t3 t4];
x=a.*t.^3+b.*t.^2+c.*t+d;

N=100;
bool=false(N,1);
for i=1:100
   ti=randn(1)*10;
   xi=cubeinterp(t(1),t(2),t(3),t(4),x(1),x(2),x(3),x(4),ti);
   xi2=a.*ti.^3+b.*ti.^2+c.*ti+d;
   bool(i)=aresame(xi,xi2,1e-10);
end


reporttest('CUBEINTERP doubly complex scalar case',all(bool(i)))



