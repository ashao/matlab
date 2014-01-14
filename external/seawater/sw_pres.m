
function pres = sw_pres(DEPTH,LAT)

% SW_PRES    Pressure from depth
%===========================================================================
% SW_PRES   $Revision: 1.5 $  $Date: 1994/10/11 01:23:32 $
%           Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE:  pres = sw_pres(depth,lat)
%
% DESCRIPTION:
%    Calculates pressure in dbars from depth in meters.
%    A FUNCTION OBTAINED BY INVERTING THE EMPIRICAL DEPTH
%    FUNCTION GIVEN BY FOFONOFF IN 1985
%
% INPUT:  (all must have same dimensions)
%   depth = depth [metres]  
%   lat   = Latitude in decimal degress north [-90..+90]
%           (LAT may have dimensions 1x1 or 1xn where depth(mxn) )
%
% OUTPUT:
%  pres   = Pressure    [db]
%
% AUTHOR:  Phil Morgan 93-06-25  (morgan@ml.csiro.au)
%     -- WORLEY DECEMBER 1985 ------
% adapted for matlab by A.H. Orsi, June 1995
%
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.  
%   See the file sw_copy.m for conditions of use and licence.
%
%=========================================================================

% CALLER:  general purpose
% CALLEE:  none

%-------------
% CHECK INPUTS
%-------------
[mD,nD] = size(DEPTH);
[mL,nL] = size(LAT);
if mL==1 & nL==1
  LAT = LAT*ones(size(DEPTH));
end %if  

if (mD~=mL) | (nD~=nL)              % DEPTH & LAT are not the same shape
     if (nD==nL) & (mL==1)          % LAT for each column of DEPTH
        LAT = LAT( ones(1,mD), : ); %     copy LATS down each column
                                    %     s.t. dim(DEPTH)==dim(LAT)
     else
        error('press85.m:  Inputs arguments have wrong dimensions')
     end %if
end %if

Transpose = 0;
if mD == 1  % row vector
   DEPTH   =  DEPTH(:);
   LAT     =  LAT(:);
   Transpose = 1;
end %if

%-------------
% BEGIN
%-------------

DEG2RAD = pi/180;
X       = (sin(abs(LAT)*DEG2RAD)).^2.;

%     THE GRAVITY ANOMALY USING Z INSTEAD OF P FOR AN APPROXIMATION

      GR=9.780318.*(1.0+(5.2788E-3+2.36E-5.*X).*X) + 1.092E-6.*DEPTH;
      DP=GR.*DEPTH;

%     SCALE DP TO ACCOUNT FOR THE SCALING DONE IN 4TH ORDER LINEAR
%     MODEL USED TO OBTAIN THE COEFFICIENTS

      DP=DP/1.0E3;
      pres=(((-2.641561E-7.*DP+2.524876E-5).*DP+2.267149E-2).*DP+1.028361E2).*DP;

% C1      = 5.92E-3+X.^2*5.25E-3;
% pres    = ((1-C1)-sqrt(((1-C1).^2)-(8.84E-6*DEPTH)))/4.42E-6;

if Transpose
   pres = pres';
end %if

return
%===========================================================================
