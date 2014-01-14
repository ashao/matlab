function [ outx outy outz ] = procdup(x,y,z,tol,method)
% PROCDUP: Finds all points within a search radius and returns values for the duplicate points
%	Input: 
%		x,y,z: Input data (vectors Nx1)
%		method: function handle to process duplicates (e.g. @max, @mean)
%	Algorithm:
%		1) Calculate distance of one point to all the others
%		2) Find all points within the tol (xr,yr,zr)
%		3) Process duplicates: [Mean(xr), Mean(yr) method(zr)]
%		4) Delete all points that were processed
%		5) Return vector of outx, outy, outz

flag = 1;
counter=0;
while flag
	counter=counter+1;
	npts=length(x);

	diffx=ones(npts,1)*x(1)-x;
	diffy=ones(npts,1)*y(1)-y;
	dist=sqrt(diffx.^2+diffy.^2);
	
	idx=find(dist<tol);
	if length(idx)>1
		xr=x(idx);
		yr=y(idx);
		zr=z(idx);
		outx(counter)=mean(xr);
		outy(counter)=mean(yr);
		outz(counter)=method(zr);
	else
		outx(counter)=x(1);
		outy(counter)=y(1);
		outz(counter)=z(1);
	end
	x(idx)=[];
	y(idx)=[];
	z(idx)=[];
	flag=~isempty(x);
end
	
	
