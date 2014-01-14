
function [chgv] = sm_chgref(V,P,PREF)

% SM_CHGREF Change of Reference Level for Geostrophic Velocity 
%=========================================================================
% SM_CHGREF  Date: 07/18/1995
%          . Sabine Mecking
%
% USAGE:  [chgv]= sm_chgref(V,P,PREF)
%
% DESCRIPTION:
%   Geostrophic Velcoity referenced to a new pressuere level by
%   subtracting this levels velocity
%
% INPUT:  (all must have same dimensions)
%   V = Velocity Field 
%       (dim mxn: rows of V must relate to increasing pressure,
%        columns of V must relate to horizontal direction)
%   P = Pressure    
%       (P may be column or row vector with length m)
%   PREF = New Reference Pressure [same units as P]
%
% OUTPUT:
%  chgv = Velocities Relative to New Reference Pressure 
%         [same units as V]
%
% AUTHOR:  Sabine Mecking 95-07-18  (mecking@ocean.washington.edu)
%
% DISCLAIMER:
%
% REFERENCE: S. Pond & G.Pickard  2nd Edition 1986
%            Introductory Dynamical Oceanogrpahy
%            Pergamon Press Sydney.  ISBN 0-08-028728-X
%
%=========================================================================

%
% CALLER: general purpose
% CALLEE: ---

%----------------------
% CHECK INPUT ARGUMENTS
%----------------------
if nargin ~=3
   error('sm_chgref.m: Must pass 3 parameters')
end %if

% CHECK V,P,PREF dimensions and verify consistent
[mv,nv] = size(V);
[mp,np] = size(P);
[mpref,npref]=size(PREF);

% CHECK THAT PREF IS SCALAR  
if mpref~=1 | npref~=1     % PREF is not a scalar
   error('check_pref: PREF is not a scalar')
end %if

% CHECK OPTIONAL SHAPES FOR P
if np==mv & mp==1          % P is row vector with same cols as rows of V
   % shape ok
elseif mp==mv & np==1      % P is column vector with same rows as V
   % shape ok
elseif mp==mv & np==nv     % P is a matrix size(V)
   P=P(:,1)                % copy one col of P 
else
   error('check_p: P has wrong dimensions')
end %if
mp = length(P);
 
% CHECK IF PREF HAS A REASONABLE VALUE IN THE RANGE OF THE PRESSURE
% VALUES
if PREF<P(1)             % PREF too small
   error('check_pref: PREF is smaller than smallest P')
elseif PREF>P(mp)          % PREF too big
   error('check_pref: PREF is bigger than biggest P')
end %if 
%***check_p


%------
% BEGIN
%------
[m,n]=size(V);
ps=find(P<=PREF);
pb=find(P>=PREF);

is=ps(length(ps));
ib=pb(1);

vref=0.5*(V(is,:)+V(ib,:));
vref=vref(ones(1,m),:);   % copy down each column
chgv=V-vref;

return
%--------------------------------------------------------------------
