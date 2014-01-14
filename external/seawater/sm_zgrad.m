function [pp_grad,p_ave] = sm_zgrad(PP,P,LAT)
%
% SM_ZGRAD   Gradient of a Propeperty in vertical (z) direction  
%===========================================================================
% SM_PGRD Sabine Mecking, July 1996
%
% USAGE:  [pp_grad,p_ave] = sm_zgrad(PP,P,LAT) 
%
% DESCRIPTION:
%    Calculates vertical gradients of property at midpoint of pressure 
%    coordinates. Treats each column (or row respectively) of PP separately.
%
%    Note: If no latitude values are passed pressure is assumed to be equal to
%    depth (i.e. no conversion with sw_dpth). Thus, this function can also 
%    be used to calculate gradients directly in z-space by passing vertcial
%    disctance (or any other) in the pressure variable.  
%
% INPUT:  (all must have same dimensions MxN)
%   PP    = property (MxN)   
%   P     = pressure [db] 
%           (MxN) or (Mx1) or (1XN)
%
%   OPTIONAL:
%     LAT = Latitude in decimal degrees north [-90 ... +90]
%           May have dimensions 1x1 or 1xn where PP (mxn)
%           (will calc d(z) instead of d(p) in numerator)
%
% OUTPUT:
%   pp_grad  =  property gradient (M-1xN) 
%   p_ave =  mid pressure between P grid (M-1xN) [db] (or units of P)
%
% AUTHOR:  Sabine Mecking 96-07-08 (mecking@ocean.washington.edu)
%
% REFERENCES: similar to sw_bfrq.m
%=========================================================================

% CALLER:  general purpose
% CALLEE:  sw_dpth 


%-------------
% CHECK INPUTS
%-------------
if nargin<2
   error('sm_zgrad.m: Must pass at least 2 parameters.')
end %if
if nargin<3
  LAT=[];
end

% CHECK PP,P dimensions and verify consistent
[mpp,npp] = size(PP);
[mp,np] = size(P);
  

% CHECK OPTIONAL SHAPES FOR P
if np==npp & mp==1          % P is row vector with same cols as PP 
   P = P( ones(1,mp), : ); %   Copy down each column.
elseif mp==mpp & np==1      % P is column vector
   P = P( :, ones(1,np) ); %   Copy across each row
elseif mp==mpp & np==npp     % P is a matrix size(PP)
   % shape ok 
else
   error('check_stp: P has wrong dimensions')
end %if
[mp,np] = size(P);
  
% IF ALL ROW VECTORS ARE PASSED THEN LET US PRESERVE SHAPE ON RETURN.
Transpose = 0;
if mp == 1  % row vector
   P           =  P';
   PP          =  PP';

   Transpose = 1;
end %if
[mp,np] = size(P);
[mpp,npp] = size(PP);

% IF LAT PASSED THEN VERIFY DIMENSIONS
if ~isempty(LAT)
   [mL,nL] = size(LAT);
   if mL==1 & nL==1
      LAT = LAT*ones(size(PP));
   end %if

   if (mpp~=mL) | (npp~=nL)              % PP & LAT are not the same shape
       if (npp==nL) & (mL==1)           % copy LATS down each column
          LAT = LAT( ones(1,mpp), : );  % s.t. dim(PP)==dim(LAT)
       else
          error('sm_zgrad.m:  Inputs arguments have wrong dimensions')
       end %if
   end %if
end %if
%***check_input

%------
% BEGIN
%------
if ~isempty(LAT)
   Z = sw_dpth(P,LAT);
else
   Z = P;
end %if

[m,n] = size(Z);
iup   = 1:m-1;
ilo   = 2:m;
z_ave = (Z(iup,:)+Z(ilo,:) )/2;
p_ave = (P(iup,:)+P(ilo,:) )/2;
pp_up = PP(iup,:);
pp_lo = PP(ilo,:);
z_up = Z(iup,:);
z_lo = Z(ilo,:);
 
dif_pp    = pp_up - pp_lo;
%dif_z = diff(Z);    % gives wrong sign
dif_z = z_up - z_lo;
pp_grad   = dif_pp ./ dif_z;

  
if Transpose
  pp_grad   = pp_grad';
  p_ave   =   p_ave';
end %if
return
%-------------------------------------------------------------------
