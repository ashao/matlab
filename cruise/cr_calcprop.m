function [ cruise ] = cr_calcprop( cruise )
% CR_CALCPROP Calculates inferred properties (e.g. density) from cruises
% Input:
%	cruise (struct) with fields: ctdprs, ctdsal, ctdtmp
% Output:
%	cruise (struct) with additional fields:
%		sw_dens: Seawater density

cruise.sw_dens=sw_dens(cruise.ctdsal,cruise.ctdtmp,cruise.ctdprs);
