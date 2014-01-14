function [hcont,hcont2,hlab]=jcontour(x,y,z,arg4,arg5,arg6,arg7,arg8)
%JCONTOUR  Contouring with labelled and unlabelled contours.
%
%   JCONTOUR draws nonlabelled contours of one style and labelled
%   contours with a second style.
%
%   JCONTOUR(X,Y,Z,V1,V2,STY1,STY2), where V1 and V2 are vectors,
%   draws nonlabelled contours with style STY1 at levels specified
%   in V1, and labelled contours with STY2 at the V2 levels. 
%
%   If V1 and V2 are scalars, they specify the contour interval for
%   the nonlabelled and labelled contours, respectively.
%
%   X and Y may both be vectors, or both matrices, or one of each.
%
%   The SYTLE strings follow the format specified in LINESTYLE.  Thus
%   STY1='2b--' draws blue dotted lines of width 2.
%
%   JCONTOUR(X,Y,Z,V1,V2,STY1) uses SYT1 for all contours.
%  
%   JCONTOUR(X,Y,Z,V1,V2) and JCONTOUR(X,Y,Z,V), like CONTOUR, uses 
%   Z values and the current color map to determine line colors.
%  
%   JCONTOUR(X,Y,Z,V,STY) draws contours with style STY at levels 
%   in V (or with an increment V if V is a scalar), with no labelling.
%
%   JCONTOUR(..., 'nolabels') supresses labelling of contours, and may
%   be used to draw unlabelled contours in two different styles.
%
%   [H1,H2,HLAB]=JCONTOUR(X,Y,Z,V1,V2,STY1,STY2) outputs vector 
%   handles to the labelled and nonlabelled contours (H1 and H2
%   respectively) as well as to the label text (HLAB).	
%	
%   Example:  [x,y,z]=peaks; jcontour(x,y,z,1,3,'D2--','2k');
%
%   Usage:  jcontour(x,y,x,v,sty);
%           jcontour(x,y,x,v1,v2,sty);
%           jcontour(x,y,x,v1,v2,sty1,sty2);
%           [h,h1,hlab]=jcontour(x,y,x,v1,v2,sty1,sty2);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2009 J.M. Lilly --- type 'help jlab_license' for details        
na=nargin;
bstyin=0;
bsty2in=0;
bv1in=0;
bv2in=0;
hcont=[];hcont2=[];hlab=[];

%if there are more than three input arguments, look for 
%style strings corresponding to the labeled and labelless contours

if na>3
    ntharg=eval(['arg',int2str(na)]);
elseif na==3
    ntharg=z;
elseif na==2
    ntharg=y;
elseif na==1
    ntharg=x;
end
blabel=1;
if strcmp(ntharg,'nolabels')
  blabel=0;
  na=na-1;
end


if na>3
	for i=4:na
		itharg=eval(['arg',int2str(i)]);
		if ischar(itharg)
			if bstyin==0
				sty=itharg;
				bstyin=1;
			elseif bstyin==1
				sty2=itharg;
				bsty2in=1;
			end
			if na>i
				argadvance(nargin,i+1);
			end
		end
	end
end
na=na-bstyin-bsty2in;
if na>=4
	minz=min(min(z));
	maxz=max(max(z));
	v=arg4;
	if length(v)==1;
		v=(floor(minz/v)*v:v:ceil(maxz/v)*v);
	end
end
if na==5
	v2=arg5;
	if length(v2)==1
		v2=(floor(minz/v2)*v2:v2:ceil(maxz/v2)*v2);
	end
	bv2in=1;
	if bsty2in
		index=true(size(v));
		for j=1:length(v)
			if any(v(j)==v2)
				index(j)=0;
			end
        end
		v=v(index);
	end
end


washold=0;
if ishold
	ax=axis;
	washold=1;
end

if length(x)~=length(y)
	if size(z,2)~=length(x)
		z=conj(z');
	end
else 
	%disp('Make sure Z is oriented correctly')
end

if ~(min(size(x))==1 && min(size(y))==1)
	if min(size(x))==1
		x=ones(length(y),1)*conj(x(:))';
	elseif min(size(y))==1
		y=row2col(y)*ones(1,length(x));
	end
end


if na<=3
	if bstyin
		[c,hcont]=contour(x,y,z,'k');
		linestyle(hcont,sty);
	else
		[c,hcont]=contour(x,y,z,'k');
    end
else
	if bsty2in
		if bstyin	
			if ~isempty(v)
			   [c,hcont]=contour(x,y,z,v,'k');
               linestyle(hcont,sty);
			end
		else
			if ~isempty(v)
		       [c,hcont]=contour(x,y,z,v);
			end
		end
		hold on
		if ~isempty(v2)
		   [c,hcont2]=contour(x,y,z,v2,'k');
           linestyle(hcont2,sty2);
		end
		if washold, axis(ax),end
		if ~isempty(v2) && blabel
		   hlab=clabel(c,hcont2,v2);
		end
	else
		if bstyin	
			if ~isempty(v)
			   [c,hcont]=contour(x,y,z,v,'k'); 
			   linestyle(hcont,sty);
			end
		else
			if ~isempty(v)
			   [c,hcont]=contour(x,y,z,v);
			end
		end
		hold on
		if bv2in
			if bstyin	
				if ~isempty(v2)
				   [c,hcont2]=contour(x,y,z,v2,'k'); 
				   linestyle(hcont2,sty);
				end
			else
				if ~isempty(v2)
				   [c,hcont2]=contour(x,y,z,v2);
				end
			end
			if ~isempty(v2) && blabel
			   hlab=clabel(c,hcont2,v2);
			end
		end
		if washold, axis(ax), end
	end

end
set(gca,'box','on')



function[evalme]=argadvance(na,n)

if nargin==1
	n=2;
end
evalme=[];
for j=n:na
	evalme=[evalme 'arg',int2str(j-1),'=arg',int2str(j),';'];
end
