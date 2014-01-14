function [acpot]=accelpot(sa,te,pr,ga,gvalue,plevel)

% function [acpot]=accelpot(sa,te,pr,ga,gvalue,plevel)
%
% Given matrices or vectors of sa, te, pr, and ga (density), Computes the 
% acceleration potential, or Montgomery function (Montgomery, 1937). acpot
% see page 21 of Reid (1965, Interm. Waters of the Pac. Ocean, TJHP) on the
% density surface with value gvalue, and relative to the isobaric surface with 
% value plevel;
%
% unit = geopotential meters (Gill, 1982, p. 46) [1gpm=9.8m^2s^-2=9.8Jkg^-1]
%
% 
% (See also function gamsurfval for details on interpolating routine)
%
% A.H. Orsi  August, 1995
% Modified 4/4/97 G. C. Johnson 
%

%
% Check for multiple density values
%

[m,n]=size(sa);

[nsurf,dummy] = size(gvalue);
if nsurf>1
!echo select only one value of Density!! ';
gvalue=gvalue(1);
end;
nsurf;

% integrate vertically to get specific volume and gepotential anomalies

sva = sw_svan(sa,te,pr);
gpa = sw_gpan(sa,te,pr);

% interpolate geopotential anomaly on pressure surface

dhlvl=interp1(pr,gpa,plevel);

% interpolate values on density surfaces station by station

for i=1:n

ii=find(ga(:,i)<=gvalue); if length(ii)>0, ii=ii(1); end % if
jj=find(ga(:,i)>=gvalue); if length(jj)>0, jj=jj(1); end % if

if (length(ii)>0 & length(jj)>0 & ii==jj);

%prsgn(i)=ga(ii,i); % seems wrong =>
prsgn(i)=pr(ii,i);  % changed by Sabine 10/7/98 
svagn(i)=sva(ii,i);
gpagn(i)=gpa(ii,i);

elseif (length(ii)>0 & length(jj)>0 & ii~=jj)

prsgn(i)=interp1(ga(ii:jj,i),pr(ii:jj),gvalue);
svagn(i)=interp1(ga(ii:jj,i),sva(ii:jj,i),gvalue);
gpagn(i)=interp1(ga(ii:jj,i),gpa(ii:jj,i),gvalue);

elseif (length(ii)==0 | length(jj)==0);

prsgn(i)=NaN;
svagn(i)=NaN;
gpagn(i)=NaN;

end % if

end % for

% make units of Acceleration Potential in geopotential meters
%       , e.g. see Gill (1982), 1 gpm = 9.8 J/kg=9.8m^2/s^

frst=dhlvl-gpagn;
scnd=prsgn.*svagn*10000;
acpot=(frst+scnd)/9.8;

