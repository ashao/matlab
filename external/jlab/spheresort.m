function[varargout]=spheresort(varargin)
%SPHERESORT  Sorted great circle distances to nearby points on a sphere.
%
%   Computing great circle distances on the sphere between two sets of
%   points --- data points and a fixed grid, say --- is computationally 
%   expensize.  SPHERESORT speeds this up substantially. 
%
%   [DS,XS,YS]=SPHERESORT(LAT,LON,LATO,LONO,CUTOFF) returns the great
%   circle distances DS between data points at locations LAT, LON and 
%   nearby grid points located at LATO, LONO, sorting in order of distance. 
%
%   XS and YS are the corresponding positions of the sorted data points 
%   in a local tangent plane about each grid point. 
%
%   SPHERESORT only computes distances for nearby points.  CUTOFF is the
%   maximum distance (in kilometers) for which we wish the great circle
%   distance to be computed.
% 
%   LAT and LON are arrays of the same size into data point locations.
%   LATO and LONO are arrays of length M and N, say, specifying the
%   latitudes and longitudes of an M x N matrix of grid points, i.e.
%
%       LATO =  [LATO_1;    LONO= [LONO_1 LONO_2 ... LONO_N]. 
%                LATO_2; 
%                  ...
%                LATO_M]
%
%   The output arrays are then each M x N x P arrays of column vectors, 
%   where P is the maximum number of points in the CUTOFF neighborhood at
%   any grid point.  Points farther away are filled with NANs.
%
%   DS gives the distances of data points less than CUTOFF kilometers 
%   from the (m,n)th grid point, sorted in order of increasing distance.  
%
%   XS and YS are also M x N x P, and give the coordinates of the 
%   corresponding data points in a Cartesian expansion about the (m,n)th
%   grid point, in kilometers.  See LATLON2XY for details.
%   _________________________________________________________________
% 
%   Additional input parameters
%
%   Let's say some additional variables Z1, Z2,...,ZK are given at
%   the data locations LAT, LON.  Then 
%   
%   [DS,XS,YS,Z1S,Z2S,...,ZKS]=
%
%                SPHERESORT(LAT,LON,Z1,Z2,...,ZK,LATO,LONO,CUTOFF);
%
%   also returns the sorted values of these variables.
%
%   Z1S, Z2S,...,ZKS are the same size as the other output arguments.  
%   Z1S then gives the value of Z1 at data points no more than CUTOFF 
%   kilometers from the (m,n)th grid point, etc.
%
%   To output an index into the data point locations, use
%  
%   [DS,XS,YS,INDEX]=SPHERESORT(LAT,LON,1:LENGTH(LAT(:)),LATO,LONO,CUTOFF).
%   _________________________________________________________________
%
%   See also TWODSORT, POLYSMOOTH.
%
%   'spheresort --t' runs a test.
%
%   Usage: ds=spheresort(lat,lon,lato,lono,cutoff);
%          [ds,xs,ys]=spheresort(lat,lon,lato,lono,cutoff);
%          [ds,xs,ys,z1s,z2s]=spheresort(lat,lon,z1,z2,lato,lono,cutoff);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008--2009 J.M. Lilly --- type 'help jlab_license' for details
 

%
%   Algorithm options
%
%   SPHERESORT can use two different algorithms, one which optimizes
%   memory and another which optimizes speed.  These give identical
%   results.
%
%   SPHERESORT(...,'memory') optimizes memory, the default.  
%
%   SPHERESORT(...,'speed') optimizes speed.  If you try this and you 
%   get an `Out of Memory' error, use the default algorithm.
%   _________________________________________________________________

%          [ds,xs,ys,index]=spheresort(lat,lon,lato,lono,cutoff);

if strcmp(varargin{1}, '--t')
    spheresort_test,return
end

lat=varargin{1};
lon=varargin{2};
if isstr(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='memory';
end

P=inf;
    
lato=varargin{end-2};
lono=varargin{end-1};
cutoff=varargin{end};
args=varargin(3:end-3);

if ~aresame(size(lat),size(lon))
    error('LAT and LON must be the same size.')
end

cutoff=abs(cutoff);

for i=1:nargin
    varargout{i}=[];
end
indexo=find(isfinite(lat)&isfinite(lon));
if ~isempty(indexo)
    vcolon(lat,lon);
    vindex(lat,lon,indexo,1);
else
    disp(['No finite data values.']), return
end

vcolon(lato,lono);
lono=lono';


[d,xd,yd,indexd]=spheresort_memory(lat,lon,lato,lono,cutoff);

varargout{1}=d;
varargout{2}=xd;
varargout{3}=yd;

K=nargout;
for k=1:length(args)
    temp=cell(length(lato),length(lono));
    for i=1:size(indexd,1)
        for j=1:size(indexd,2)
            if ~isempty(indexd{i,j})
                temp{i,j}=args{k}(indexo(indexd{i,j}));
            end
        end
    end
    varargout{3+k}=temp;
end
%I think this can be done with no cell loop, just stick into long array
%and the stick this with clever index into varargout
%easy just cell2col (so it can't be i,j to start) and then index in
cells=varargout;
xs=cells{1};
L=zeros(size(xs));
for i=1:size(xs,1)
    for j=1:size(xs,2)
        L(i,j)=length(xs{i,j});
    end
end
maxL=maxmax(L);

for k=1:K
    varargout{k}=vzeros(size(xs,1),size(xs,2),maxL,'nan');
end

for i=1:size(xs,1)
    for j=1:size(xs,2)
        if L(i,j)>=1
            for k=1:K
                varargout{k}(i,j,1:L(i,j))=cells{k}{i,j};
            end
        end
    end
end

disp('SPHERESORT finished.')
     

function [d,xd,yd,indexd]=spheresort_memory(lat,lon,lato,lono,cutoff)

dlat_cutoff=jrad2deg(cutoff./radearth);

sd=cell(length(lato),length(lono));
indexd=cell(length(lato),length(lono));
d=cell(length(lato),length(lono));
xd=cell(length(lato),length(lono));
yd=cell(length(lato),length(lono));

for j=1:length(lato)
     disp(['SPHERESORT computing latitude band ' int2str(j) ' of ' int2str(length(lato)) '.'])
     index=find(abs(lat-lato(j))<=dlat_cutoff);
     if ~isempty(index)
                  
         [latj,lonj]=vindex(lat,lon,index,1);
        
         lonomat=vrep(lono,length(index),1);
         latomat=lato(j)+zeros(size(lonomat));
                  
         [latjmat,lonjmat,indexmat]=vrep(latj,lonj,index,length(lono),2); 
         
         [dj,xjmat,yjmat]=vzeros(size(latomat),'nan');
                  
         %Form a ``wedge'' of nearby points; speeds things up quite a bit
         dlon=abs(angle(exp(sqrt(-1)*jdeg2rad(lonjmat-lonomat))));
         
         maxlat=min(abs(lato(j))+jrad2deg(frac(cutoff,radearth)),90);
         neari=find(dlon.*radearth.*cosd(maxlat)<=cutoff);
         
         if ~isempty(neari)
             [xjmat(neari),yjmat(neari),dj(neari)]=latlon2xy(latjmat(neari),lonjmat(neari),latomat(neari),lonomat(neari));
         end

         dj(dj>cutoff)=nan;
             
         %Sort and return 
         if size(dj,1)>1
             [dj,sorter]=sort(dj);
             
             sorterjj=vrep(1:size(sorter,2),size(sorter,1),1);
             indexsorter=sub2ind(size(xjmat),sorter,sorterjj);

             xjmat=xjmat(indexsorter);
             yjmat=yjmat(indexsorter);
             indexmat=indexmat(indexsorter);
         end
       
         last=find(~isnan(vsum(dj,2)),1,'last');
         if ~isempty(last)
             vindex(dj,lonjmat,latjmat,indexmat,1:last,1);

             for i=1:length(lono)
                   ii=find(isfinite(dj(:,i)));
                   if ~isempty(ii)
                       d{j,i}=dj(ii,i);
                       xd{j,i}=xjmat(ii,i);
                       yd{j,i}=yjmat(ii,i);
                       indexd{j,i}=indexmat(ii,i);
                   end
             end
         end
     end
end

function[]=spheresort_test

N=100;
lat=rand(N,1)*180-90;
lon=rand(N,1)*360;

lono=(0:5:360);
lato=(-80:5:80);

cutoff=1000;

tic;[d,xd,yd,latd,lond]=spheresort(lat,lon,lat,lon,lato,lono,cutoff);etime1=toc;

tic;
d2=zeros(size(d));

for i=1:length(lato)
    for j=1:length(lono)
        d2(i,j,:)=spheredist(lato(i),lono(j),latd(i,j,:),lond(i,j,:));
    end
end
etime2=toc;

%disp(['SPHERESORT was ' num2str(etime2./etime1) ' times faster than a simple loop.'])

bool1=zeros(length(lato),length(lono));
bool2=zeros(length(lato),length(lono));
for i=1:length(lato)
    for j=1:length(lono)
        bool1(i,j)=all(d(i,j,:)<=cutoff|isnan(d(i,j,:)));
        bool2(i,j)=aresame(d(i,j,:),d2(i,j,:),1e-6);
    end
end

reporttest('SPHERESORT all distances less than or equal to cutoff',allall(bool1))
reporttest('SPHERESORT verify distances',allall(bool2))




