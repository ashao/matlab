function[]= fixlabels(arg1,arg2)
%FIXLABELS  Specify precision of axes labels.
%
%   FIXLABELS(PREC) rewrites the X and Y axis labels of the current
%   axes to precision PREC.  PREC denotes the tenths place of the
%   highest-precision digit, so PREC=-2 includes the hundredths place.
%
%   Unlike standard matlab labelling, FIXLABELS uses zeros instead of
%   blank spaces to the right of the decimal point.
%
%   FIXLABELS([XPREC,YPREC]) uses different precision for the X and Y
%   axes (note brackets in the function call).
%
%   FIXLABELS(H,PREC) rewrites the labels of axes H.  H may also be a
%   vector of axes handles, in which case PREC may either be a scalar
%   or a vector PREC=[X1PREC Y1PREC X2PREC ...] .
%
%   FIXLABELS with no arguments rewrites all X and Y axis labels of
%   the current figure to the default (current) precision for that
%   axis, but using zeros instead of blank spaces to the right of the
%   decimal point.
%
%   Note: make sure the tickmarks are where you want them in the
%   hardcopy before calling FIXLABELS (since it sets the 'TICKMODE'
%   and 'TICKLABELMODE' properties to 'MANUAL').
%
%   See also: vnum2str, DIGIT
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 1999--2010 J.M. Lilly --- type 'help jlab_license' for details  

%	Note unexplained Matlab quirk-- setting YTICKLABELMODE
%	to MANUAL sometimes causes uppermost label to be invisible.

if strcmp(arg1,'--t')
    return
end
h1=gca;
bprec=0;
if nargin==0
	h=findallaxes(gcf);
elseif nargin==1
	if ishandle(arg1)
		h=arg1;
	else
		h=gca;
		prec=arg1;
		bprec=1;
	end
elseif nargin==2
	h=arg1;
	prec=arg2;
	bprec=1;
	if all(size(prec)>1)
		prec=prec';
		prec=prec(:);
	end
end

if bprec
	if length(prec)==1
		prec=prec*ones(length(h)*3,1);
	end
end

for i=1:length(h)
	axes(h(i));
	xtick=get(gca,'xtick')';
	ytick=get(gca,'ytick')';
	ztick=get(gca,'ztick')';

	%find precision used by current labels assuming a decimal point
	nx=size(get(gca,'xticklabel'),2);
	ny=size(get(gca,'yticklabel'),2);
	nz=size(get(gca,'zticklabel'),2);
	
	if ~isempty(xtick)&&~isempty(get(gca,'xticklabel'))
		if ~bprec
			if all(xtick==floor(xtick))
				%no decimal point
				xlab=vnum2str(xtick,'spaces');
			else	
				%minus signs are difficult
				if any(xtick<0)
					nd=floor(log10(max(abs(xtick))))+1;
					if nd<1, nd=1;end
					nx1=nd+2-nx;	
				else
					nd=floor(log10(max(xtick)))+1;
					if nd<1, nd=1;end
					nx1=nd+1-nx;	
				end
				xlab=vnum2str(xtick,nx1,'spaces');
			end
		else
			xlab=vnum2str(xtick,prec(3*i-2),'spaces');
		end
		xlab=flushleft(xlab);
		set(gca,'xtickmode','manual')
		set(gca,'xticklabelmode','manual')
		set(gca,'xticklabel',xlab)
    end
	if ~isempty(ytick)&&~isempty(get(gca,'yticklabel'))
		if ~bprec
			if all(ytick==floor(ytick))
				%no decimal point
				ylab=vnum2str(ytick,'spaces');
			else
				%minus signs are difficult
				if any(ytick<0)
					nd=floor(log10(max(abs(ytick))))+1;
					if nd<1, nd=1;end
					ny1=nd+2-ny;	
				else
					nd=floor(log10(max(ytick)))+1;
					if nd<1, nd=1;end
					ny1=nd+1-ny;	
				end
				ylab=vnum2str(ytick,ny1,'spaces');
			end
		else
			ylab=vnum2str(ytick,prec(3*i-1),'spaces');
		end
       		ylab=flushleft(ylab);
		set(gca,'ytickmode','manual')
		set(gca,'yticklabelmode','manual')
		set(gca,'yticklabel',ylab)
    end

	if ~isempty(ztick)&&~isempty(get(gca,'zticklabel'))
		if ~bprec
			if all(ztick==floor(ztick))
				%no decimal point
				ylab=vnum2str(ztick,'spaces');
			else
				%minus signs are difficult
				if any(ztick<0)
					nd=floor(log10(max(abs(ztick))))+1;
					if nd<1, nd=1;end
					nz1=nd+2-nz;	
				else
					nd=floor(log10(max(ztick)))+1;
					if nd<1, nd=1;end
					nz1=nd+1-nz;	
				end
				zlab=vnum2str(ztick,nz1,'spaces');
			end
		else
			zlab=vnum2str(ztick,prec(3*i),'spaces');
		end
       		zlab=flushleft(zlab);
		set(gca,'ztickmode','manual')
		set(gca,'zticklabelmode','manual')
		set(gca,'zticklabel',zlab)
	end
end

axes(h1)


function[h]=findallaxes(fignum)
%FINDALLAXES    Returns handles to all axes children.
%       FINDALLAXES returns a vector of handles to all axes childen
%       regardless of which figure is their parent.
%
%       FINDALLAXES(FIGNUM) returns a vector of handles to axes
%       children of Figure FIGNUM. 
%
%       See also ALLCHILD, FINDALL

%       Author: J.M. Lilly, 1/15/98


if nargin==1
        h=get(fignum,'children');
else
        h=[];
        for i=1:length(get(0,'children'))
                h=[h;get(i,'children')];
        end
end
bool=false(size(h));
for i=1:length(h)
        if strcmp(get(h(i),'type'),'axes')
                bool(i)=1;
        end
end
h=h(bool); 


