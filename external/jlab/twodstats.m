function[mat,xmid,ymid,num,std]=twodstats(xdata,ydata,zdata,xbin,ybin,str)
%TWODSTATS  Mean, variance, and covariance of functions of two variables.
%
%   TWODSTATS computes the first- and second-order statistics of a function
%   of two variables in prescribed bins.  This function may either be a 
%   scalar-valued or vector-valued quantity at each point. 
% 
%   An example of a scalar-valued dataset is temperature as a function of
%   latitude and longitude. An example of a vector-valued dataset is wind
%   or current velocity as a function of latitude and longitude.
%
%   TWODSTATS uses a fast (exact) algorithm which is particularly efficient 
%   for large arrays.  TWODSTATS can be two orders of magnitude faster than
%   the obvious way of sorting the data in bins using explicit loops.   
%   __________________________________________________________________
%  
%   Mean and standard deviation of a scalar-valued function 
%
%   MZ=TWODSTATS(X,Y,Z,XBIN,YBIN) where X, Y and Z are arrays of the same
%   length, forms the mean of Z over the XY plane.  
%
%   The values of Z are sorted into bins according to the associated 
%   (X,Y) value, with bin edges specified by XBIN and YBIN, and the mean
%   of all finite values of Z in each bin is returned as MZ.
%  
%   If XBIN and YBIN are length N and M, respectively, then MZ is of 
%   size M-1 x N-1.  Bins with no data are assigned a value of NAN.
%
%   XBIN and YBIN must be monotonically increasing. 
%  
%   [MZ,XMID,YMID]=TWODSTATS(...) optionally returns the midpoints XMID
%   and YMID of the bins.
%
%   [MZ,XMID,YMID,NUMZ]=TWODSTATS(...) also returns the number of good
%   data points in each of the (X,Y) bins.  NUMZ is the same size as MZ.
%
%   [MZ,XMID,YMID,NUMZ,STDZ]=TWODSTATS(...) also returns the standard 
%   deviation of Z in the (X,Y) bins.  STDZ is the same size as MZ.
%   __________________________________________________________________
%  
%   Mean and covariance of a vector-valued function 
%   
%   TWODSTATS can also be used to analyze a function which contains more
%   than one value at each (X,Y) point.  
%
%   If Z represents a vector with K components, then Z should have the same
%   size as X and Y in all but its last dimension, which will be length K.
%
%   MZ=TWODSTATS(X,Y,Z,XBIN,YBIN) then returns MZ, containing the mean 
%   values of each component of Z in each bin.  If M and N are the lengths 
%   of XBIN and YBIN, MZ is of size M-1 x N-1 x K. 
%
%   [MZ,XMID,YMID,NUMZ,COVZ]=TWODSTATS(...) returns the full covariance
%   matrix COVZ in each of the bins.  As the covariance of Z is K x K, the 
%   size of the output matrix COVZ is M-1 x N-1 x K x K.
%   __________________________________________________________________
%
%   Additional comments
%
%   You can use TWODSTATS for fast binning of data over the plane.  For
%   the case in which Z is so sparsely distributed over X and Y, such 
%   that no bins will have more than one entry, the mean in each bin
%   is just the value of the data point in the bin.  
%
%   TWODSTATS(...,'slow') uses a slow algorithm which uses less memory.  
%   By default, TWODSTATS uses a fast but memory-intensive algorithm.  
%   Use this flag if you get an out-of-memory error.  
%   __________________________________________________________________
%   
%   See also TWODHIST, TWODMED, TWODSORT.
%
%   'twodstats --t' runs a test.
%
%   Usage: mz=twodstats(x,y,z,xbin,ybin);
%          [mz,xmid,ymid]=twodstats(x,y,z,xbin,ybin);
%          [mz,xmid,ymid,numz]=twodstats(x,y,z,xbin,ybin);
%          [mz,xmid,ymid,numz,stdz]=twodstats(x,y,z,xbin,ybin);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2011 J.M. Lilly --- type 'help jlab_license' for details    


if strcmp(xdata,'--t')
   twodstats_test;return
end

if nargin==5
  str='fast';
end

vcolon(xbin,ybin);
if any(diff(xbin)<0)
  error('XBIN must be monotonically increasing')
end
if any(diff(ybin)<0)
  error('YBIN must be monotonically increasing')
end
if ~aresame(size(xdata),size(ydata))
     error('X and Y should have the same size.')
end
    
if nargout>4
    stdflag=1;
else
    stdflag=0;
end
vcolon(xdata,ydata);
K=length(zdata(:))/length(xdata);
if ~isint(K)
     error('Size of Z is not compatible with size of X and Y.')
else
    zdata=reshape(zdata,length(xdata),K);
end   

bool1=isfinite(xdata)&isfinite(ydata)&isfinite(sum(zdata,2));
bool2=xdata<xbin(end)&xdata>xbin(1)&ydata<ybin(end)&ydata>ybin(1);
bool=bool1&bool2;

xdata=xdata(bool);
ydata=ydata(bool);
zdata=zdata(bool,:);


if ~isempty(zdata)
    if ~isempty(strfind(str,'fast'))
      [mat,num,std]=twodstats_fast(xdata,ydata,zdata,xbin,ybin,stdflag);
    elseif ~isempty(strfind(str,'slow'))
      [mat,num,std]=twodstats_slow(xdata,ydata,zdata,xbin,ybin,stdflag);
    end
else
    disp('Warning: Data contains only NANs and / or INFs.')
    mat=0*osum(ybin(1:end-1),xbin(1:end-1)); 
    num=mat;
    std=mat;
end


if nargout>1
  xmid=(xbin+vshift(xbin,1,1))./2;
  xmid=xmid(1:end-1);
end
if nargout>2
  ymid=(ybin+vshift(ybin,1,1))./2;
  ymid=ymid(1:end-1);
end

index=find(num==0);
if ~isempty(index)
    mat(index)=nan;
    if stdflag
        %for i=1:size(std,3)...
        std(index)=nan;
    end
end

function[mat,num,covmat]=twodstats_fast(xdata,ydata,zdata,xbin,ybin,stdflag)
%vsize(xdata,ydata,zdata,xbin,ybin,stdflag)
[num,mat,covmat]=vzeros((length(ybin)-1)*(length(xbin)-1),1,'nan');
mat=vrep(mat,size(zdata,2),2);
covmat=vrep(vrep(covmat,size(zdata,2),2),size(zdata,2),3);

xnum=bindata(xbin,xdata);
ynum=bindata(ybin,ydata);

index=sub2ind([length(ybin)-1,length(xbin)-1],ynum,xnum);

%figure,plot(xnum,ynum,'.')
%vsize(index,xdata,ydata,zdata)


if ~isempty(index)
    [index,sorter]=sort(index);
    zdata=double(zdata(sorter,:));
    [L,ia,ib,numblock]=blocklen(index);
    
    num(index(ia))=L(ia);
    vswap(num,0,nan);
     
    if 1
        %Cumsum has problem: Force double precision for good results
        cumsumzdata=cumsum(zdata,1);
        mat(index(ia),:)=cumsumzdata(ib,:)-cumsumzdata(ia,:)+zdata(ia,:);
    else
        %For testing purposes
        for i=1:length(ia)
            mat(index(ia(i)),:)=sum(zdata(ia(i):ib(i)),:);
        end
    end
    mat=mat./vrep(num,size(mat,2),2);

    if stdflag
        %This trickery simply puts the mean back where it belongs in the column
        blockmean=zeros(length(index),size(zdata,2)); 
        blockmean(1,:)=mat(index(ia(1)),:);
        blockmean(ia(2:end),:)=mat(index(ia(2:end)),:)-mat(index(ia(1:end-1)),:);
        blockmean=cumsum(blockmean,1);
        
        zprime=zdata-blockmean;
        if size(zdata,2)==1        
            zvar=abs(zprime).^2;
        else
            zprime2=conj(permute(zprime,[1 3 2]));
            zprime=vrep(zprime,size(zdata,2),3);
            zprime2=vrep(zprime2,size(zdata,2),2);
            zvar=zprime.*zprime2;
        end
        
        if 1
            %Cumsum has problem: Force double precision for good results
            cumsumzvar=cumsum(double(zvar),1);
            covmat(index(ia),:,:)=cumsumzvar(ib,:,:)-cumsumzvar(ia,:,:)+zvar(ia,:,:);
        else
            %Just so you know what the above line actually does
            for i=1:length(ia)
                covmat(index(ia(i)),:,:)=sum(zvar(ia(i):ib(i)),:,:);
            end
        end
        
        covmat=covmat./vrep(vrep(num,size(covmat,2),2),size(covmat,3),3);
                
        if size(zdata,2)==1 
            covmat=sqrt(covmat);
        end
            
    end
end
vswap(num,nan,0);
mat=reshape(mat,length(ybin)-1,length(xbin)-1,size(zdata,2));
num=reshape(num,length(ybin)-1,length(xbin)-1);
covmat=reshape(covmat,length(ybin)-1,length(xbin)-1,size(zdata,2),size(zdata,2));


function[mat,num,std]=twodstats_slow(xdata,ydata,zdata,xbin,ybin,stdflag)
vcolon(xdata,ydata,zdata,xbin,ybin);
index=find(isfinite(xdata)&isfinite(ydata)&isfinite(zdata));
vindex(xdata,ydata,zdata,index,1);

mat=0*osum(ybin,xbin); 
num=0*osum(ybin,xbin); 
    
if stdflag 
    std=0*osum(ybin,xbin); 
else
    std=[];
end

[xbinb,ybinb]=vshift(xbin,ybin,1,1);
for i=1:length(xbin)
   for j=1:length(ybin)
         index=find(xdata>xbin(i)&xdata<=xbinb(i)&ydata>ybin(j)&ydata<=ybinb(j));
         if ~isempty(index)
             mat(j,i)=vmean(zdata(index),1);
             num(j,i)=length(index);
             if stdflag
                 std(j,i)=vstd(zdata(index),1);
             end
         end
   end
end
mat=mat(1:end-1,:);
mat=mat(:,1:end-1); 
num=num(1:end-1,:);
num=num(:,1:end-1);
    
if stdflag
    std=std(1:end-1,:);
    std=std(:,1:end-1);
end

function[]=twodstats_test
twodstats_test1;
twodstats_test2;

function[]=twodstats_test1
load npg2006
use npg2006
lono=[-19:.05:-16.5];
lato=[42.8:.05:44.3];


tic;
[mat1,xmid,ymid,num1,std1]=twodstats(lon,lat,t,lono,lato,'fast');
dt1=toc;
tic
[mat2,xmid,ymid,num2,std2]=twodstats(lon,lat,t,lono,lato,'slow');
dt2=toc;

%index=isfinite(mat1)&isfinite(mat2);
index=1:length(mat1(:));
%size(find(index))

%figure,pcolor(mat1),shading flat
%figure,pcolor(mat2),shading flat

bool1=aresame(mat1(index),mat2(index),1e-10);
bool2=aresame(std1(index),std2(index),1e-10);
bool3=aresame(num1(index),num2(index),1e-10);
reporttest('TWODSTATS fast vs. slow algorithm for npg2006, mean',bool1)
reporttest('TWODSTATS fast vs. slow algorithm for npg2006, std',bool2)
reporttest('TWODSTATS fast vs. slow algorithm for npg2006, num',bool3)
disp(['TWODSTATS fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])


function[]=twodstats_test2
L=10000;
xdata=3*abs(rand(L,1));
ydata=3*abs(rand(L,1));
zdata=randn(L,1);
xbin=(0:.1:2);
ybin=(0:.2:2);
tic;
[mat1,xmid,ymid,num1,std1]=twodstats(xdata,ydata,zdata,xbin,ybin,'fast');
dt1=toc;
tic
[mat2,xmid,ymid,num2,std2]=twodstats(xdata,ydata,zdata,xbin,ybin,'slow');
dt2=toc;
bool1=aresame(mat1,mat2,1e-10);
bool2=aresame(std1,std2,1e-10);
bool3=aresame(num1,num2,1e-10);
reporttest('TWODSTATS fast vs. slow algorithm, mean',bool1)
reporttest('TWODSTATS fast vs. slow algorithm, std',bool2)
reporttest('TWODSTATS fast vs. slow algorithm, num',bool3)
disp(['TWODSTATS fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])

xdata=-4*rand(L,1);
ydata=-5*rand(L,1);
zdata(1:round(L/10):end)=nan;
xbin=(-2:.1:0);
ybin=(-3:.2:0);
tic;[mat1,xmid,ymid,num1,std1]=twodstats(xdata,ydata,zdata,xbin,ybin,'fast');
dt1=toc;
tic;[mat2,xmid,ymid,num2,std2]=twodstats(xdata,ydata,zdata,xbin,ybin,'slow');
dt2=toc;
bool1=aresame(mat1,mat2,1e-10);
bool2=aresame(std1,std2,1e-10);
bool3=aresame(num1,num2,1e-10);
reporttest('TWODSTATS fast vs. slow algorithm, with negative bins and NANs, mean',bool1)
reporttest('TWODSTATS fast vs. slow algorithm, with negative bins and NANs, std',bool2)
reporttest('TWODSTATS fast vs. slow algorithm, with negative bins and NANs, num',bool3)



disp(['TWODSTATS fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])


function[]=twodstats_test3
load nordicdrifters
use nordicdrifters
cell2col(lon,lat,cv);
[x,y]=latlon2xy(lat,lon,70,10);
z=rot(-pi/4)*(x+sqrt(-1)*y);
x=real(z);
y=imag(z);
xbin=[-1000:20:1000];
ybin=[-600:20:200];



 [mz,xmid,ymid,numz,stdz]=twodstats(x,y,cv,xbin,ybin);
 [mz,xmid,ymid,numz,cov]=twodstats(x,y,[real(cv) imag(cv)],xbin,ybin);




[xg,yg]=meshgrid(xmid,ymid);

ellipseplot(kappa,lambda,th,xg+sqrt(-1)*yg,'skip',[1 5])

load jetdrifters
use jetdrifters

cv=latlon2uv(num,lat,lon);
 
xbin=[-80:1/2:-50];
ybin=[30:1/2:45];
[mz,xmid,ymid,numz,covz]=twodstats(lon,lat,[real(cv) imag(cv)],xbin,ybin);


%Diagonalize the covariance matrix
[d1,d2,th]=specdiag(covz(:,:,1,1),covz(:,:,2,2),covz(:,:,1,2));

%Convert to ellipse parameters, as in notes
a=sqrt(d1);b=real(sqrt(d2));  %bn can sometimes have a very small spurious negative part

[xg,yg]=meshgrid(xmid,ymid);
spd=sqrt(vsum(squared(mz),3));

figure,
contourf(xmid,ymid,spd,40),nocontours
hold on,latratio(37.5),caxis([0 115])


[kappa,lambda]=ab2kl(a,b);
kappa(spd>30)=nan;
ellipseplot(kappa/200,lambda,th,xg+sqrt(-1)*yg,get(gca,'dataaspectratio'),'skip',[1 1], 'G')

[kappa,lambda]=ab2kl(a,b);
kappa(spd<30)=nan;
h=ellipseplot(kappa/400,lambda,th,xg+sqrt(-1)*yg,get(gca,'dataaspectratio'),'skip',[1 1], 'w');
linestyle -h h 4w
h=ellipseplot(kappa/400,lambda,th,xg+sqrt(-1)*yg,get(gca,'dataaspectratio'),'skip',[1 1], 'm');
linestyle -h h 2m


title('Variance Ellipses over the Gulf Stream')
h=colorbar;
axes(h),ylabel('Speed of Mean Current (cm/s)')
orient landscape
fontsize 14 12 12 12
cd_figures
print -dpng gulfellipses.png


