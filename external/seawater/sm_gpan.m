function [ga,ap] = sm_gpan(S,T,P,pref)

% SM_GPAN    Geopotential anomaly
%=========================================================================
% SM_GPAN  modification of $Revision: 1.3 $  $Date: 1994/10/10 05:01:00 $
%                           Copyright (C) CSIRO, Phil Morgan 1992.
%
% USAGE:  [gpan,acpo]= sw_gpan(S,T,P,pref)
%
% DESCRIPTION:
%   Geopotential Anomaly calculated as the integral of svan from the
%   the sea surface to the bottom.  Thus RELATIVE TO SEA SURFACE.
%
% INPUT:  (all must have same dimensions)
%   S = salinity    [psu      (PSS-78)]
%   T = temperature [degree C (ITP-68)]
%   P = Pressure    [db]
%       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
%   pref = (optional) reference pressure
%          by default geopotential is refernced to the surface assuming
%          constant svan from the lowest pressure to 0 db (cf. "top")
%
% OUTPUT:
%  gpan = Geopotential Anomaly   [m^3 kg^-1 Pa == m^2 s^-2 == J kg^-1]
%  acpo = Acceleration Potential [m^3 kg^-1 Pa == m^2 s^-2 == J kg^-1]
%
% AUTHOR:  Phil Morgan 92-11-05  (morgan@ml.csiro.au)
%          Sabine Mecking 98-08-06  (mecking@ocean.washington.edu)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCE: S. Pond & G.Pickard  2nd Edition 1986
%            Introductory Dynamical Oceanogrpahy
%            Pergamon Press Sydney.  ISBN 0-08-028728-X
%
% Note that older literature may use units of "dynamic decimeter' for above.
%
% Adapted method from Pond and Pickard (p76) to calc gpan rel to sea 
% surface whereas P&P calculated relative to the deepest common depth.
%
% CHANGED: 8/5/1998 by Sabine Mecking
%          routine returns also acceleration potentail now
%          ac = ga + p*svan  (Reid,Intermediate Water of the Pacific Ocean,
%                             John Hopkins Press, 1965)
%          
%          routine allows to specify reference pressure which is used for
%          geopotential anomaly as well as acceleration potential
%
%=========================================================================

%
% CALLER: general purpose
% CALLEE: sw_svan.m 

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin <3
   error('sm_gpan.m: Must pass at least parameters')
end %if

% CHECK S,T,P dimensions and verify consistent
[ms,ns] = size(S);
[mt,nt] = size(T);
[mp,np] = size(P);
[mpr,npr] = size(pref);

  
% CHECK THAT S & T HAVE SAME SHAPE
if (ms~=mt) | (ns~=nt)
   error('check_stp: S & T must have same dimensions')
end %if

% CHECK OPTIONAL SHAPES FOR P
if     mp==1  & np==1      % P is a scalar.  Fill to size of S
   P = P(1)*ones(ms,ns);
elseif np==ns & mp==1      % P is row vector with same cols as S
   P = P( ones(1,ms), : ); %   Copy down each column.
elseif mp==ms & np==1      % P is column vector
   P = P( :, ones(1,ns) ); %   Copy across each row
elseif mp==ms & np==ns     % PR is a matrix size(S)
   % shape ok 
else
   error('check_stp: P has wrong dimensions')
end %if
[mp,np] = size(P);

if mpr==1 & npr==1
   % shape ok
else
   error('check_stp: pref has wrong dimensions')
end

  
% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0; 
if mp == 1  % row vector 
  P       =  P(:); 
  T       =  T(:); 
  S       =  S(:);   
  Transpose = 1;
end %if
%***check_stp

%------
% BEGIN
%------
db2Pascal  = 1e4;
[m,n]      = size(P);
svan       = sw_svan(S,T,P);
mean_svan  = 0.5*(svan(2:m,:) + svan(1:m-1,:) );

if n==1
   top = svan(1,1).*P(1,1)*db2Pascal;
else
   top = svan(1,:).*P(1,:)*db2Pascal;
end %if

%press_diff = diff(P);

delta_ga   = (mean_svan.*diff(P))*db2Pascal;
ga         = cumsum([top; delta_ga]);

if nargin>3          % reference pressure given
  % find out whether columns of P vector are all the same
  P_even_spacing = 1;
  for i=1:n
    if P(i,:)~=P(1,:) 
      P_even_spacing = 0;
    end
  end
 
  if P_even_spacing        % can do all columns at once
    ga_ref=interp1(P(:,1),ga,pref);
  else
    for i=1:n
      ga_ref(i)=interp1(P(:,i),ga(:,i),pref);
    end
  end

  % reference geopotential anomaly to reference pressure
  ga_ref = ga_ref( ones(1,ms), : );   %   Copy down each column.
  ga=ga_ref-ga;

end

% calc.ulate accerlation potential
ap=ga+P*db2Pascal.*svan;

if Transpose
   ga = ga';
   ap = ap';
end %if

return
%--------------------------------------------------------------------
