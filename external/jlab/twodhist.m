function[mat,xmid,ymid,index]=twodhist(xdata,ydata,xbin,ybin,str)
%TWODHIST  Two-dimensional histogram.
%
%   MAT=TWODHIST(X,Y,XBIN,YBIN) where X and Y are arrays of the same
%   length, creates a two-dimensional histogram MAT with bin edges
%   specified by XBIN and YBIN. 
%
%   TWODHIST uses a fast (exact) algorithm which is particularly 
%   efficient for large arrays.
%  
%   The (X,Y) data points are sorted according to their X-values, which 
%   determine a column within MAT, and their Y-values, which determine 
%   a row within MAT.  
%
%   If XBIN and YBIN are length N and M, respectively, then MAT is of
%   size M-1 x N-1.
%
%   XBIN and YBIN must be monotonically increasing. 
%
%   [MAT,XMID,YMID]=TWODHIST(...) optionally returns the midpoints XMID
%   and YMID of the bins.
%
%   [MAT,XMID,YMID,INDEX]=TWODHIST(...) also returns INDEX, an array of
%   the same size as X and Y giving the index into the matrix MAT of 
%   each data point. 
%
%   TWODHIST(...,'slow') uses a slow algorithm which uses less memory.  
%   By default, TWODHIST uses a fast but memory-intensive algorithm.  
%   Use this flag if you get an out-of-memory error.  
%
%   See also TWODMED, TWODSTATS, TWODSORT.
%
%   'twodhist --t' runs a test.
%
%   Usage: mat=twodhist(x,y,xbin,ybin);
%          [mat,xmid,ymid]=twodhist(x,y,xbin,ybin);
%          [mat,xmid,ymid,index]=twodhist(x,y,xbin,ybin);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details    

  
if strcmp(xdata,'--t')
   twodhist_test;return
end

if nargin==4
  str='fast';
end

xbin=xbin(:);
ybin=ybin(:);
if any(diff(xbin)<0)
  error('XBIN must be monotonically increasing')
end
if any(diff(ybin)<0)
  error('YBIN must be monotonically increasing')
end

bool1=isfinite(xdata)&isfinite(ydata);
bool2=xdata<xbin(end)&xdata>xbin(1)&ydata<ybin(end)&ydata>ybin(1);
bool=bool1&bool2;

xdata=xdata(bool);
ydata=ydata(bool);

if ~isempty(xdata)
    if ~isempty(strfind(str,'fast'))
      [mat,index]=twodhist_fast(xdata,ydata,xbin,ybin);
    else
      mat=twodhist_slow(xdata,ydata,xbin,ybin);
    end
else
    disp('Warning: Data contains only NANs and / or INFs.')
    mat=0*osum(ybin(1:end-1),xbin(1:end-1)); 
end

if nargout>1
  xmid=(xbin+vshift(xbin,1,1))./2;
  xmid=xmid(1:end-1);
end
if nargout>2
  ymid=(ybin+vshift(ybin,1,1))./2;
  ymid=ymid(1:end-1);
end

function[mat,index]=twodhist_fast(xdata,ydata,xbin,ybin)

[xnum,xi,xmid]=bindata(xbin,xdata);
[ynum,yi,ymid]=bindata(ybin,ydata);

mat=zeros([length(ybin)-1,length(xbin)-1]);
index=nan*zeros(size(xdata));

nani=(~isnan(xnum)&~isnan(ynum));

if sum(nani(:))>0
    index(nani)=sub2ind([length(ybin)-1,length(xbin)-1],ynum(nani),xnum(nani));
    [indexsorted,sorter]=sort(index(nani));
    
    [L,ia]=blocklen(indexsorted);
    mat(indexsorted(ia))=L(ia);
end





function[mat]=twodhist_slow(xdata,ydata,xbin,ybin)
mat=0*osum(ybin,xbin); 
[xbinb,ybinb]=vshift(xbin,ybin,1,1);
for i=1:length(xbin)
   for j=1:length(ybin)
         mat(j,i)=length(find(xdata>xbin(i)&xdata<=xbinb(i)&...
			      ydata>ybin(j)&ydata<=ybinb(j)));
   end
end
mat=mat(1:end-1,:);
mat=mat(:,1:end-1); 

function[]=twodhist_test
L=10000;
xdata=3*abs(rand(L,1));
ydata=3*abs(rand(L,1));
xbin=(0:.1:2);
ybin=(0:.2:2);
tic;
mat1=twodhist(xdata,ydata,xbin,ybin,'fast');
dt1=toc;
tic
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
dt2=toc;
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm',bool)
disp(['TWODHIST hist algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])

xdata=-3*abs(rand(L,1));
ydata=-3*abs(rand(L,1));
xbin=(-2:.1:0);
ybin=(-2:.2:0);
mat1=twodhist(xdata,ydata,xbin,ybin,'fast');
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm, negative bins',bool)

xdata=randn(L,1);
ydata=randn(L,1);
xbin=(-2:.1:2);
ybin=(-2:.2:2);
mat1=twodhist(xdata,ydata,xbin,ybin,'fast');
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm, crossing zero',bool)




