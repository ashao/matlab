function[ticks,ticklabs]=timelabel(yearf,style,n)
%TIMELABEL  Put month, day, or hour labels on a time axes.
%
%   TIMELABEL(YEARF,STYLE) puts labels of style SYTLE on a time
%   axis with YEAR.FRAC format (see NUM2YF).  Valid STYLEs are:
%  
%       'month' for month
%	'day' for day of current year
%	'cumday' for day of year relative to beginning of record
%	'hour' for the hour of the current day
%	'cumhour' for the hour relative to beginning of record
% 
%   YEARF must be a column vector.
%
%   TIMELABEL(YEARF,SYTLE,N), where N is an interger, labels every
%   Nth month (day); the default labels every twelvth month (day).
%
%   TIMELABEL defines the first day of the year, January 1, as day 0
%   for the purposese of 'cumday'.  Similarly, for 'cumhour' the first
%   hour is defined as hour 0.
%    
%   See also NUM2YF, DATENUM, DATESTR
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000, 2004 J.M. Lilly --- type 'help jlab_license' for details      
    
if size(yearf,2)>1
  error('YEARF must be a column vector')
end

dyearf=yearf(end)-yearf(1);

if nargin==1	
	if dyearf>=1,
	   style='none';
	elseif dyearf>=0.5  &&  dyearf<1
	   style='month';
	elseif (dyearf*365)>=6  &&  dyearf<0.5
           style='day';  
	else
           style='cumhour';  
	end
end

X=12;


if nargin<3	
	if strcmp(style,'month')
	   n=floor(dyearf*12/X);
	   if n<1, n=1./floor(X/dyearf/12); end
        elseif strcmp(style,'day') || strcmp(style,'cumday')
	   n=floor(dyearf*365/X);
	   if n<1, n=1./floor(X/dyearf/365); end
        elseif strcmp(style,'cumhour') || strcmp(style,'hour')
	   n=floor(dyearf*365*24/X);
	   if n<1, n=1./floor(X/dyearf/365/24); end
	end
end

num=yf2num(yearf);
[yr,mo,dy,hr]=datevec(num);


if strcmp(style,'cumday')
	yearday=floor(num)-floor(datenum(yr(1),mo(1),dy(1)));
elseif strcmp(style,'day')
	yearday=floor(num)-datenum(yr(1),1,1);
elseif strcmp(style,'cumhour') || strcmp(style,'hour')
	hrd=diff(hr);
	ii=find(hrd<0);
	if ii(end)~=length(hrd)
	   hrd(ii)=hrd(ii+1);
	else
	   hrd(ii(1:end-1))=hrd(ii(1:end-1)+1);
	   hrd(ii(end))=hrd(ii(end)-1);
	end
	hrd(end+1)=hrd(end);
	cmhr=cumsum(hrd);
	cmhr=cmhr-cmhr(1);
elseif ~strcmp(style,'month')
  error('Unsupported label type.')
end


if strcmp(style,'month')
        moi=find(diff(mo)~=0)+1;
	if n>1
		index=1:n:length(moi);
	else
		index=1:1:length(moi);
	end
	moi=moi(index);
%	molabs=['J>F>M>A>M>J>J>A>S>O>N>D>'];
%	molabs=[reshape(molabs,2,12)]';
	molabs='JFMAMJJASOND';
	molabs=reshape(molabs,1,12)';
	ticks=yearf(moi);
	ticklabs=molabs(mo(moi),:);	
elseif strcmp(style,'cumday') || strcmp(style,'day')
	dyi=[1;find(diff(yearday)~=0)+1];
	if any(yearday==0);
	   bi=min(find(yearday(dyi)==0));
	else
	  bi=1;
	end
	%Adjust for sample rate less than 1 day
	n=round(n./(yearday(dyi(2))-yearday(dyi(2)-1)));
        index=bi:n:length(dyi);
	if dyi(index(end))>length(yearf)
	  index=index(1:end-1);
	end
	dyi=dyi(index);
	yearday=yearday(dyi);
	daylabs=vnum2str(floor(yearday),'spaces');
	ticks=yearf(dyi);
	ticklabs=daylabs;
elseif strcmp(style,'hour') ||  strcmp(style,'cumhour')
	hri=[1;find(diff(cmhr)~=0)+1];
	if any(hr==0);
	   bi=min(find(hr(hri)==0));
	else
	  bi=1;
	end
	%Adjust for sample rate less than 1 hour
	n=round(n./(cmhr(hri(2))-cmhr(hri(2)-1)));
        index=bi:n:length(hri);
	if hri(index(end))>length(yearf)
	  index=index(1:end-1);
	end
	hri=hri(index);
	ticks=yearf(hri);
	if strcmp(style,'hour')
	  ticklabs=vnum2str(hr(hri),'spaces');
	else
  	  ticklabs=vnum2str(cmhr(hri),'spaces');
	end
end

set(gca,'xtick',ticks);
set(gca,'xticklabel',ticklabs)

if nargout==1
	clear ticklabs
elseif nargout==0
	clear ticklabs ticks
end


